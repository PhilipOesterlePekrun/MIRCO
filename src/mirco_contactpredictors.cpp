#include "mirco_contactpredictors.h"

#include <Teuchos_TimeMonitor.hpp>

namespace MIRCO
{
  void ContactSetPredictor(ViewVectorInt_d& activeSet0, ViewVector_d& xv0, ViewVector_d& yv0,
      ViewVector_d& b0, double zmax, double Delta, double w_el, const ViewMatrix_d topology,
      const ViewVector_d meshgrid)
  {
    static auto timer = Teuchos::TimeMonitor::getNewCounter("ContactSetPredictor()");
    FenceForTiming();
    Teuchos::TimeMonitor monitor(*timer);

    const int N = topology.extent(0);

    const double deltaContact = Delta + w_el - zmax;
    int n0 = 0;
    Kokkos::parallel_reduce(
        N * N,
        KOKKOS_LAMBDA(const int a, int& local_sum) {
          if (topology(a / N, a % N) >= -deltaContact) local_sum++;
        },
        n0);

    activeSet0 = ViewVectorInt_d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "activeSet0"), n0);
    xv0 = ViewVector_d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "xv0"), n0);
    yv0 = ViewVector_d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "yv0"), n0);
    b0 = ViewVector_d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "b0"), n0);

    // Compact candidates in global cell order without a shared counter.
    Kokkos::parallel_scan(
        N * N, KOKKOS_LAMBDA(const int a, int& offset, const bool final) {
          const int i = a / N;
          const int j = a % N;
          const double topology_a = topology(i, j);
          if (topology_a >= -deltaContact)
          {
            const int aa = offset++;
            if (final)
            {
              activeSet0(aa) = a;
              xv0(aa) = meshgrid(i);
              yv0(aa) = meshgrid(j);
              // Note: b0 = \overbar{u} + w_el = Delta + w_el - (zmax - topology_a);
              b0(aa) = topology_a + deltaContact;
            }
          }
        });
    FenceForTiming();
  }

}  // namespace MIRCO
