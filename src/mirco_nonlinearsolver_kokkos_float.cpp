#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

// tmp
#include "tmpHelpers/Timer.hpp"

void swap_entries(ViewVectorInt_d v, const int i, const int j)
{
  if constexpr (HOSTONLY)
  {
    std::swap(v(i), v(j));
  }
  else
  {
    auto vi = Kokkos::subview(v, i);
    auto vj = Kokkos::subview(v, j);

    ViewScalarInt_d tmp("tmp");

    Kokkos::deep_copy(tmp, vi);
    Kokkos::deep_copy(vi, vj);
    Kokkos::deep_copy(vj, tmp);
  }
}

void MIRCO::NonLinearSolver::solve(const ViewMatrix_d matrix, const ViewVector_d b0_d,
    ViewVector_d& p_d, int& activeSetSize, double nnlstol, int maxiter)
{
  ///tmp{
  
  #ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
    std::cout<<"nonlinear line36\n";
  #endif
  #ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
    std::cout<<"nonlinear line39\n";
  #endif
  ///tmp}
  
  using minloc_t = Kokkos::MinLoc<double, int, MemorySpace_ofDefaultExec_t>;
  using minloc_value_t = typename minloc_t::value_type;
  
  ExecSpace_Default_t ExecSpace_Default;

  const std::string thisFctName = "MIRCO::NonLinearSolver::Solve()";
  ScopedTimer sc1(thisFctName);

  // 2^{-52}; double machine precision, i.e. smallest number such that 1.0 + eps > 1.0
  constexpr double eps = 1.0 / (1ULL << 52);
  
  // positive infinity for double
  ///constexpr double infty = 0x7FF0000000000000;
  const int n0 = b0_d.extent(0);

  ViewVector_d w_d("w_d", n0);

  ViewVectorInt_d activeInactiveSet("activeInactiveSet", n0);
  // These need to be int pointers for atomic_fetch_add() to work reliably //# says chatgpt. maybe
  // check this further
  Kokkos::View<int*, Device_Default_t> counterActive("counterActive", 1);
  Kokkos::deep_copy(counterActive, 0);
  Kokkos::View<int*, Device_Default_t> counterInactive("counterInactive", 1);
  Kokkos::deep_copy(counterInactive, 0);

  Kokkos::parallel_for(
      "build initial active set", n0, KOKKOS_LAMBDA(const int i) {
        if (p_d(i) >= nnlstol)
        {
          // Note: atomic_fetch_add() returns the old value
          activeInactiveSet(Kokkos::atomic_fetch_add(&counterActive(0), 1)) = i;
        }
        else
        {
          activeInactiveSet(n0 - 1 - (Kokkos::atomic_fetch_add(&counterInactive(0), 1))) = i;
        }
      });
  Kokkos::deep_copy(activeSetSize, Kokkos::subview(counterActive, 0));

  bool init = false;
  if (activeSetSize == 0)
  {
    init = true;

    // Hp - b0 for p = 0
    Kokkos::parallel_for(
        thisFctName + "__firstParallel", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int i) { w_d(i) = -b0_d(i); });
  }
  // iter is incremented in while(2)
  int iter = 0;
  while (1)
  {
    if ((init && (activeSetSize == n0)) || iter >= maxiter) break;
    ////minloc_value_t minloc_w_i;
    /// double min_w_i = infty;
    /// int argmin_w_i = -1;
    Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_w_i_d("minloc_w_i_d");

    Kokkos::parallel_reduce(//# THE CUDA ISSUE IS HERE
        "min_w", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
          const double wi = w_d(i);
          if (wi < ml.val)
          {
            ml.val = wi;
            ml.loc = i;
          }
        },
        minloc_t(minloc_w_i_d));
        
    minloc_value_t minloc_w_i;
    Kokkos::deep_copy(minloc_w_i, minloc_w_i_d);

    if (init && (minloc_w_i.val >= -nnlstol)) break;

    // #ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
    //     std::cout << " LAPACK";
    // #endif
    // #ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
    //     std::cout << " BLAS";
    // #endif
    if (init)
    {
      ////minloc_value_t minloc_w_iInactive;
      Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_w_iInactive_d("minloc_w_iInactive_d");
      Kokkos::parallel_reduce(
          "min_w", Kokkos::RangePolicy<ExecSpace_Default_t>(activeSetSize, n0),
          KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
            const double wi = w_d(activeInactiveSet(i));
            if (wi < ml.val)
            {
              ml.val = wi;
              ml.loc = i;
            }
          },
          minloc_t(minloc_w_iInactive_d));

      minloc_value_t minloc_w_iInactive;
      Kokkos::deep_copy(minloc_w_iInactive, minloc_w_iInactive_d);
      
      swap_entries(activeInactiveSet, minloc_w_iInactive.loc, activeSetSize);
      ++activeSetSize;
    }
    else
      init = true;
    
    while (2)
    {
      ++iter;

///std::cout<<"line153\n";
      // Compact versions of H and b0, i.e. H_I and \overbar{u}_I in line 6 of Algorithm 3,
      // (Bemporad & Paggi, 2015)
      /// std::cout << "activeSetSize=" << activeSetSize << ", n0=" << n0 << "\n";

      ViewVector_d b0s_compact_d("b0_compact", activeSetSize);
      
      if (activeSetSize > 1) {
#if(!GESV_USE_FP32)
        ViewMatrix_d H_compact_d("H_compact", activeSetSize, activeSetSize);

        Kokkos::parallel_for(
            "define compact views", activeSetSize, KOKKOS_LAMBDA(const int i) {
              const int row = activeInactiveSet(i);
              b0s_compact_d(i) = b0_d(row);
              for (int j = 0; j < activeSetSize; ++j)
              {
                const int col = activeInactiveSet(j);
                H_compact_d(i, j) = matrix(row, col);
              }
            });

        ViewVectorInt_d ipiv("ipiv", activeSetSize);

///ScopedTimer sc2("gesv()");
StandardTimer st2("st gesv()");
st2.start();
        // Solve H_I s_I = b0_I; b0s_compact_d becomes s_I
        KokkosLapack::gesv(ExecSpace_Default, H_compact_d, b0s_compact_d, ipiv);
st2.stop();

#else
      ViewVectorFloat_d b0s_compact_float_d("b0s_compact_float", activeSetSize);
      
      ViewMatrixFloat_d H_compact_d("H_compact_float", activeSetSize, activeSetSize);
///std::cout<<"line188\n";
        Kokkos::parallel_for(
            "define compact views float32", activeSetSize, KOKKOS_LAMBDA(const int i) {
              const int row = activeInactiveSet(i);
              b0s_compact_float_d(i) = static_cast<float>(b0_d(row));
              for (int j = 0; j < activeSetSize; ++j)
              {
                const int col = activeInactiveSet(j);
                H_compact_d(i, j) = static_cast<float>(matrix(row, col));
              }
            });

        ViewVectorInt_d ipiv("ipiv", activeSetSize);
///ScopedTimer sc2("gesv()");
StandardTimer st2("st gesv()");
st2.start();
        // Solve H_I s_I = b0_I; b0s_compact_float_d becomes s_I
        KokkosLapack::gesv(ExecSpace_Default, H_compact_d, b0s_compact_float_d, ipiv);
st2.stop();

Kokkos::parallel_for(
            "float32 to double", activeSetSize, KOKKOS_LAMBDA(const int i) {
              b0s_compact_d(i) = static_cast<double>(b0s_compact_float_d(i));
            });
#endif
          }
      else
      {
        ///std::cout<<"else if(!(activeSetSize > 1))\n";
        Kokkos::parallel_for(
            "", Kokkos::RangePolicy<ExecSpace_Default_t>(ExecSpace_Default, 0, 1), KOKKOS_LAMBDA(const int i) {
              const int ii = activeInactiveSet(0);
              b0s_compact_d(0) =
                  b0_d(ii) / matrix(ii, ii);
            });
      }

      bool allGreater = true;
      Kokkos::parallel_reduce(
          "allGreater", activeSetSize,
          KOKKOS_LAMBDA(const int i, bool& rallGreater) {
            const bool lesser = (b0s_compact_d(i) < -nnlstol);
            rallGreater = rallGreater && !lesser;
          },
          Kokkos::LAnd<bool>(allGreater));
///std::cout<<"allGreater="<<allGreater<<"\n";
      if (allGreater)
      {
        Kokkos::parallel_for(
            "p = s_compact", activeSetSize,
            KOKKOS_LAMBDA(const int i) { p_d(activeInactiveSet(i)) = b0s_compact_d(i); });

        Kokkos::parallel_for(
            "w = H p - b0", n0, KOKKOS_LAMBDA(const int i) {
              double sum = 0.0;
              for (int j = 0; j < activeSetSize; ++j)
                sum += matrix(i, activeInactiveSet(j)) * b0s_compact_d(j);
              w_d(i) = sum - b0_d(i);
            });

        break;
      }
      else
      {
///std::cout<<"line252\n";
        // Lines 8-11 of Algorithm 3, (Bemporad & Paggi, 2015)

        // Find min_i and argmin_i of alpha_i := \frac{p_i}{p_i - s_i}
        ////minloc_value_t minloc_alpha_i;
        Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t>minloc_alpha_i_d("minloc_alpha_i_d");
        Kokkos::parallel_reduce(
            "min_alpha", activeSetSize,
            KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
              if (b0s_compact_d(i) <= 0)
              {
                const double alphai = p_d(activeInactiveSet(i)) /
                                      (eps + p_d(activeInactiveSet(i)) - b0s_compact_d(i));
                if (alphai < ml.val)
                {
                  ml.val = alphai;
                  ml.loc = i;
                }
              }
            },
            minloc_t(minloc_alpha_i_d));
            
///<"line274\n";
            
        minloc_value_t minloc_alpha_i;
        Kokkos::deep_copy(minloc_alpha_i, minloc_alpha_i_d);

        Kokkos::parallel_for(
            "assemble reduced system", activeSetSize, KOKKOS_LAMBDA(const int i) {
              p_d(activeInactiveSet(i)) =
                  p_d(activeInactiveSet(i)) +
                  minloc_alpha_i_d().val * (b0s_compact_d(i) - p_d(activeInactiveSet(i)));
            });
///std::cout<<"line285\n";


        const int minloc_alpha_i_loc = minloc_alpha_i.loc;
        if (minloc_alpha_i_loc > -1)
        {
          ///Kokkos::deep_copy(ExecSpace_Default, Kokkos::subview(p_d, activeInactiveSet(minloc_alpha_i.loc)), 0.0);
          
          Kokkos::parallel_for(
            "", Kokkos::RangePolicy<ExecSpace_Default_t>(ExecSpace_Default, 0, 1), KOKKOS_LAMBDA(const int i) {
              p_d(activeInactiveSet(minloc_alpha_i_loc)) = 0.0;
            });
///std::cout<<"line290\n";
          // `swap_entries()` may not be very efficient on GPU because it does not use
          // `Kokkos::kokkos_swap()`, as that can only be used in parallel lambda. If performance is
          // a problem at this point, consider if it is possible to parallize this entire `while(2)`
          // loop--i.e. put it in a `Kokkos::parallel_for()` and place the `if(allGreater) { ...;
          // break; }` logic after the loop
          swap_entries(activeInactiveSet, minloc_alpha_i_loc, --activeSetSize);
        }
///std::cout<<"line298\n";
      }
///std::cout<<"line300\n";
    }
  }
  
  std::cout<<"Final iter="<<iter<<"\n";

  // # TODO: sort the active set region of activeInactiveSet and at the same time sort p_d and w_d
  // (if needed?) in the same way # wait but if all we need is the final total pressure and the
  // contact area, p_d does not need to be ordered, so that it would be most efficient to integrate
  // the last two functions into here. We can also make a switch or overload for whether we need the
  // (ordered)

  // ## p_d now still holds it in the same order though. It is like [0 p_1 p_2 0 0 p_5]. So, it is
  // essentially sorted in the way that the input p0_d was, and also the xv0 and such. To "sort" it,
  // we also just see where 0. But this can be unreliable, so let us instead
  // ## - same with w_d

  // p_d = p;
  // u_star = w + \overbar{u}
}
