#include "mirco_inputparameters.h"

#include "mirco_shapefactors.h"
#include "mirco_topology.h"

namespace MIRCO
{
  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, int Resolution, double InitialTopologyStdDeviation,
      double Hurst, int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag,
      bool RandomSeedFlag, std::optional<int> RandomGeneratorSeed,
      std::optional<std::string> ExportVisualizationPath)
      : InputParameters(E1, E2, nu1, nu2, Tolerance, Delta, LateralLength, Resolution,
            InitialTopologyStdDeviation, Hurst, MaxIteration, WarmStartingFlag,
            PressureGreenFunFlag, RandomSeedFlag, RandomGeneratorSeed, ExportVisualizationPath,
            true)
  {
  }

  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, int Resolution, double InitialTopologyStdDeviation,
      double Hurst, int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag,
      bool RandomSeedFlag, std::optional<int> RandomGeneratorSeed,
      std::optional<std::string> ExportVisualizationPath, bool ElasticComplianceCorrectionFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        elastic_compliance_correction_flag(ElasticComplianceCorrectionFlag),
        pressure_green_funct_flag(PressureGreenFunFlag),
        N((1 << Resolution) + 1),
        export_visualization_path(ExportVisualizationPath)
  {
    auto topology_h = CreateRmgSurface(
        Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed);
    topology = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);

    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    if (elastic_compliance_correction_flag)
    {
      shape_factor = getShapeFactor(N, PressureGreenFunFlag);
      elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    }
    grid_size = LateralLength / N;
  }

  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag,
      std::optional<std::string> ExportVisualizationPath)
      : InputParameters(E1, E2, nu1, nu2, Tolerance, Delta, LateralLength, TopologyFilePath,
            MaxIteration, WarmStartingFlag, PressureGreenFunFlag, ExportVisualizationPath, true)
  {
  }

  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag,
      std::optional<std::string> ExportVisualizationPath, bool ElasticComplianceCorrectionFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        elastic_compliance_correction_flag(ElasticComplianceCorrectionFlag),
        pressure_green_funct_flag(PressureGreenFunFlag),
        export_visualization_path(ExportVisualizationPath)
  {
    auto topology_h = CreateSurfaceFromFile(TopologyFilePath);
    N = topology_h.extent(0);
    topology = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);

    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    if (elastic_compliance_correction_flag)
    {
      shape_factor = getShapeFactor(N, PressureGreenFunFlag);
      elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    }
    grid_size = LateralLength / N;
  }

}  // namespace MIRCO
