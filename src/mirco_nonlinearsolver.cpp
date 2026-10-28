#include "mirco_nonlinearsolver.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

namespace
{
  using namespace MIRCO;
  void swapEntries(ViewVectorInt_d v, const int i, const int j)
  {
    if constexpr (std::is_same<MemorySpace_Host_t, MemorySpace_ofDefaultExec_t>())
    {
      std::swap(v(i), v(j));
    }
    else
    {
      Kokkos::parallel_for(
          1, KOKKOS_LAMBDA(const int) {
            const int tmp = v(i);
            v(i) = v(j);
            v(j) = tmp;
          });
    }
  }

  void removeEntry(ViewVector_d p, ViewVectorInt_d indices, const int i, const int last)
  {
    if constexpr (std::is_same<MemorySpace_Host_t, MemorySpace_ofDefaultExec_t>())
    {
      p(indices(i)) = 0.0;
      std::swap(indices(i), indices(last));
    }
    else
    {
      Kokkos::parallel_for(
          1, KOKKOS_LAMBDA(const int) {
            const int removed = indices(i);
            p(removed) = 0.0;
            indices(i) = indices(last);
            indices(last) = removed;
          });
    }
  }
}  // namespace

namespace MIRCO
{
  int nonlinearSolve(ViewVector_d& pf, ViewVectorInt_d& activeSetf, ViewVector_d& p,
      const ViewVectorInt_d activeSet0, const ViewMatrix_d matrix, const ViewVector_d b0,
      double nnlstol, int maxiter)
  {
    static auto timer = Teuchos::TimeMonitor::getNewCounter("nonlinearSolve()");
    FenceForTiming();
    Teuchos::TimeMonitor monitor(*timer);

    using minloc_t = Kokkos::MinLoc<double, int, MemorySpace_ofDefaultExec_t>;
    using minloc_value_t = typename minloc_t::value_type;
    using minloc_view_t = Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t>;
    using unmanaged_matrix_t = Kokkos::View<double**, Kokkos::LayoutLeft, Device_Default_t,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using matrix_policy_t = Kokkos::MDRangePolicy<ExecSpace_Default_t,
        Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;
    const std::string kokkosLabelPrefix = "nonlinearSolve(); ";

    constexpr double machine_eps = std::numeric_limits<double>::epsilon();

    const int n0 = b0.extent(0);

    ViewVector_d w(Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "w"), n0);

    ViewVectorInt_d activeInactiveSet(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "activeInactiveSet"),
        n0);

    int activeSetSize = 0;
    Kokkos::parallel_scan(
        n0,
        KOKKOS_LAMBDA(const int i, int& activeCount, const bool final) {
          const bool active = p(i) >= nnlstol;
          if (final)
          {
            // Active indices retain their order; inactive indices fill the tail in reverse.
            const int target = active ? activeCount : n0 - 1 - (i - activeCount);
            activeInactiveSet(target) = i;
          }
          if (active) ++activeCount;
        },
        activeSetSize);

    minloc_view_t minloc_w_i_d(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "minloc_w_i_d"));
    minloc_view_t minloc_w_iInactive_d(Kokkos::view_alloc(
        Kokkos::WithoutInitializing, kokkosLabelPrefix + "minloc_w_iInactive_d"));
    minloc_view_t minloc_alpha_i_d(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "minloc_alpha_i_d"));

    ViewVector_d b0s_workspace(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "b0_compact"), n0);
    ViewVectorInt_d ipiv_workspace;
    ViewVector_d H_workspace;

    bool init = false;
    if (activeSetSize == 0)
    {
      init = true;
      // Hp - b0 for p = 0
      Kokkos::parallel_for(n0, KOKKOS_LAMBDA(const int i) { w(i) = -b0(i); });
    }

    int iter = 0;
    while (true)
    {
      if ((init && (activeSetSize == n0)) || iter >= maxiter) break;

      if (init)
      {
        // A warm start first solves its initial contact set before a residual is available.
        Kokkos::parallel_reduce(
            n0,
            KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
              const double wi = w(i);
              if (wi < ml.val)
              {
                ml.val = wi;
                ml.loc = i;
              }
            },
            minloc_t(minloc_w_i_d));

        minloc_value_t minloc_w_i_h;
        Kokkos::deep_copy(minloc_w_i_h, minloc_w_i_d);

        if (minloc_w_i_h.val >= -nnlstol) break;

        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<ExecSpace_Default_t>(activeSetSize, n0),
            KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
              const double wi = w(activeInactiveSet(i));
              if (wi < ml.val)
              {
                ml.val = wi;
                ml.loc = i;
              }
            },
            minloc_t(minloc_w_iInactive_d));

        minloc_value_t minloc_w_iInactive_h;
        Kokkos::deep_copy(minloc_w_iInactive_h, minloc_w_iInactive_d);

        swapEntries(activeInactiveSet, minloc_w_iInactive_h.loc, activeSetSize);
        ++activeSetSize;
      }
      else
        init = true;

      while (true)
      {
        ++iter;

        // Compact versions of H and b0, i.e. H_I and \overbar{u}_I in line 6 of Algorithm 3,
        // (Bemporad & Paggi, 2015)
        ViewVector_d b0s_compact = Kokkos::subview(b0s_workspace, std::make_pair(0, activeSetSize));
        if (activeSetSize > 1)
        {
          unmanaged_matrix_t H_compact;
          ViewVectorInt_d ipiv;
          {
            static auto timer = Teuchos::TimeMonitor::getNewCounter("PrepareLinearSystem");
            FenceForTiming();
            Teuchos::TimeMonitor monitor(*timer);
            if (ipiv_workspace.extent(0) == 0)
              ipiv_workspace = ViewVectorInt_d(
                  Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "ipiv"), n0);
            ipiv = Kokkos::subview(ipiv_workspace, std::make_pair(0, activeSetSize));

            const std::size_t matrixEntries =
                static_cast<std::size_t>(activeSetSize) * static_cast<std::size_t>(activeSetSize);
            if (H_workspace.extent(0) < matrixEntries)
            {
              // Retain only the largest reduced system seen so far, with no extra growth reserve.
              // Release old storage first to avoid holding both dense allocations during growth.
              H_workspace = ViewVector_d();
              H_workspace = ViewVector_d(
                  Kokkos::view_alloc(Kokkos::WithoutInitializing, kokkosLabelPrefix + "H_compact"),
                  matrixEntries);
            }
            // Repack the active prefix with the current leading dimension for the LAPACK interface.
            H_compact = unmanaged_matrix_t(H_workspace.data(), activeSetSize, activeSetSize);

            Kokkos::parallel_for(
                matrix_policy_t({0, 0}, {activeSetSize, activeSetSize}),
                KOKKOS_LAMBDA(const int i, const int j) {
                  const int row = activeInactiveSet(i);
                  const int col = activeInactiveSet(j);
                  H_compact(i, j) = matrix(row, col);
                  if (j == 0) b0s_compact(i) = b0(row);
                });

            FenceForTiming();
          }

          // Solve H_I s_I = b0_I; b0s_compact becomes s_I
          {
            static auto timer = Teuchos::TimeMonitor::getNewCounter("LinearSolve");
            Teuchos::TimeMonitor monitor(*timer);
            KokkosLapack::gesv(H_compact, b0s_compact, ipiv);
            FenceForTiming();
          }
        }
        else if (activeSetSize == 1)
        {
          Kokkos::parallel_for(
              1, KOKKOS_LAMBDA(const int) {
                const int ii = activeInactiveSet(0);
                b0s_compact(0) = b0(ii) / matrix(ii, ii);
              });
        }

        bool allGreater = true;
        Kokkos::parallel_reduce(
            activeSetSize,
            KOKKOS_LAMBDA(const int i, bool& rallGreater) {
              const bool lesser = (b0s_compact(i) < -nnlstol);
              rallGreater = rallGreater && !lesser;
            },
            Kokkos::LAnd<bool>(allGreater));
        if (allGreater)
        {
          static auto timer = Teuchos::TimeMonitor::getNewCounter("UpdateSolutionAndResidual");
          Teuchos::TimeMonitor monitor(*timer);

          Kokkos::parallel_for(
              activeSetSize,
              KOKKOS_LAMBDA(const int i) { p(activeInactiveSet(i)) = b0s_compact(i); });

          Kokkos::parallel_for(
              n0, KOKKOS_LAMBDA(const int i) {
                double sum = 0.0;
                for (int j = 0; j < activeSetSize; ++j)
                  sum += matrix(i, activeInactiveSet(j)) * b0s_compact(j);
                w(i) = sum - b0(i);
              });

          FenceForTiming();
          break;
        }
        else
        {
          // Find min_i and argmin_i of alpha_i := \frac{p_i}{p_i - s_i}
          Kokkos::parallel_reduce(
              activeSetSize,
              KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
                if (b0s_compact(i) <= 0)
                {
                  const double alphai = p(activeInactiveSet(i)) /
                                        (machine_eps + p(activeInactiveSet(i)) - b0s_compact(i));
                  if (alphai < ml.val)
                  {
                    ml.val = alphai;
                    ml.loc = i;
                  }
                }
              },
              minloc_t(minloc_alpha_i_d));

          minloc_value_t minloc_alpha_i_h;
          Kokkos::deep_copy(minloc_alpha_i_h, minloc_alpha_i_d);

          Kokkos::parallel_for(
              activeSetSize, KOKKOS_LAMBDA(const int i) {
                p(activeInactiveSet(i)) =
                    p(activeInactiveSet(i)) +
                    minloc_alpha_i_d().val * (b0s_compact(i) - p(activeInactiveSet(i)));
              });

          if (minloc_alpha_i_h.loc > -1)
          {
            removeEntry(p, activeInactiveSet, minloc_alpha_i_h.loc, --activeSetSize);
          }
        }
      }
    }

    // Construct the final active set (the lower half of activeInactiveSet), as well as the compact
    // final pressure vector
    activeSetf = ViewVectorInt_d(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "activeSetf"), activeSetSize);
    pf = ViewVector_d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "pf"), activeSetSize);
    Kokkos::parallel_for(
        activeSetSize, KOKKOS_LAMBDA(const int i) {
          activeSetf(i) = activeSet0(activeInactiveSet(i));
          pf(i) = p(activeInactiveSet(i));
        });
    FenceForTiming();
    return iter;
  }

}  // namespace MIRCO
