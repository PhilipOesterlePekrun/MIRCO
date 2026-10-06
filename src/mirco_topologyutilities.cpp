#include "mirco_topologyutilities.h"

#include <Teuchos_TimeMonitor.hpp>
#include <cmath>

namespace MIRCO
{
  ViewVector_d CreateMeshgrid(const int ngrid, const double GridSize)
  {
    static auto timer = Teuchos::TimeMonitor::getNewCounter("_CreateMeshgrid()");
    FenceForTiming();
    Teuchos::TimeMonitor monitor(*timer);

    ViewVector_d meshgrid("CreateMeshgrid(); meshgrid", ngrid);

    const double GridSize_2 = GridSize / 2;
    Kokkos::parallel_for(
        ngrid, KOKKOS_LAMBDA(const int i) { meshgrid(i) = GridSize_2 + i * GridSize; });

    FenceForTiming();
    return meshgrid;
  }

  double GetMax(const ViewMatrix_d topology, const char* timerName)
  {
    auto timer = Teuchos::TimeMonitor::getNewCounter(timerName);
    FenceForTiming();
    Teuchos::TimeMonitor monitor(*timer);

    const int n0 = topology.extent(0);
    const int n1 = topology.extent(1);

    double zmax = -std::numeric_limits<double>::infinity();
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {n0, n1}),
        KOKKOS_LAMBDA(int i, int j, double& update) {
          const double val = topology(i, j);
          if (val > update) update = val;
        },
        Kokkos::Max<double>(zmax));

    return zmax;
  }

}  // namespace MIRCO
