#ifndef SRC_EVALUATE_H_
#define SRC_EVALUATE_H_

#include "mirco_inputparameters.h"
#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Relate the far-field displacement with pressure
   *
   * @param[out] pressure Mean pressure
   * @param[out] effectiveContactAreaFraction Effective contact area as percentage of the total area
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] LateralLength Lateral side of the surface [micrometers]
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] Tolerance Tolerance for the convergence of force
   * @param[in] MaxIteration Maximum number of iterations for the force to converge
   * @param[in] CompositeYoungs Composite Young's modulus
   * @param[in] WarmStartingFlag Warm-Starter Flag
   * @param[in] ElasticComplianceCorrection Elastic compliance correction
   * @param[in] topology Topology matrix containing heights
   * @param[in] zmax Maximum height
   * @param[in] meshgrid Meshgrid vector
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   * @param[in] ExportVisualizationPath Path to export visualization files to
   */
  void Evaluate(double& pressure, double& effectiveContactAreaFraction, const double Delta,
      const double LateralLength, const double GridSize, const double Tolerance,
      const int MaxIteration, const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection, const ViewMatrix_d topology, const double zmax,
      const ViewVector_d meshgrid, const bool PressureGreenFunFlag,
      std::optional<std::string> VisualizationExportPath = std::nullopt);

  /**
   * @brief Relate the far-field displacement with pressure, taking the parameters from an
   * InputParameters object
   *
   * @param[out] pressure Mean pressure
   * @param[out] effectiveContactAreaFraction Effective contact area as percentage of the total area
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] inputParams Object which holds the input parameters
   */
  inline void Evaluate(double& pressure, double& effectiveContactAreaFraction, const double Delta,
      const InputParameters& inputParams)
  {
    Evaluate(pressure, effectiveContactAreaFraction, Delta, inputParams.lateral_length,
        inputParams.grid_size, inputParams.tolerance, inputParams.max_iteration,
        inputParams.composite_youngs, inputParams.warm_starting_flag,
        inputParams.elastic_compliance_correction, inputParams.topology, inputParams.topology_max,
        inputParams.mesh_grid, inputParams.pressure_green_funct_flag,
        inputParams.export_visualization_path);
  }
}  // namespace MIRCO

#endif  // SRC_EVALUATE_H_
