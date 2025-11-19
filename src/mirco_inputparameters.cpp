#include "mirco_inputparameters.h"

#include "mirco_topology.h"
#include "mirco_shapefactors.h"

namespace MIRCO
{
  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, int Resolution, double InitialTopologyStdDeviation,
      double Hurst, bool RandomSeedFlag, int RandomGeneratorSeed, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        pressure_green_funct_flag(PressureGreenFunFlag),
        N((1 << Resolution) + 1)
  {
    auto topology_h = CreateRmgSurface(
        Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed);
    topology = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);
    
    shape_factor = getShapeFactor(N, PressureGreenFunFlag);
    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    grid_size = LateralLength / N;
  }

  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        pressure_green_funct_flag(PressureGreenFunFlag)
  {
    auto topology_h = CreateSurfaceFromFile(TopologyFilePath);
    N = topology_h.extent(0);
    topology = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);

    shape_factor = getShapeFactor(N, PressureGreenFunFlag);
    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    grid_size = LateralLength / N;
  }

}  // namespace MIRCO
