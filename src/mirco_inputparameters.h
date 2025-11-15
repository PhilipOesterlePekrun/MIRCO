#ifndef SRC_INPUTPARAMETERS_H_
#define SRC_INPUTPARAMETERS_H_

#include <string>

#include "mirco_kokkostypes.h"

namespace MIRCO::Input
{
  using std::optional;

  struct Pptimization_Algorithm_NNLS
  {
    optional<double> Tolerance;
    optional<int> MaxIter;
  };  // -> all

  struct Solver_Parameters
  {
    bool ElasticCorrection;
    optional<int> MaxIterations;
    optional<double> Tolerance;
    bool PressureGreenFunFlag;
    optional<bool> WarmStartingFlag;
    // std::variant<optimization_algorithm_nnls, /*etc*/> optimization_algorithm;
    std::optional<Pptimization_Algorithm_NNLS> optimization_algorithm;
  };  // -> all

  struct Material_Parameters
  {
    double E1;
    double nu1;
    double E2;
    double nu2;
    inline double composite_youngs()
    {
      return 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    };
  };  // -> compositeYoungsLateralLength: 1000.0

  struct Topology_RMG
  {
    int Resolution;
    double HurstExponent;
    double InitialTopologyStdDeviation;
    std::optional<int> RandomGeneratorSeed;
  };

  /*struct Topology_File {
    std::string TopologyFilePath;
  };
  struct Topology_Flat {
    int N;
  };*/
  using Topology_File = std::string;
  using Topology_Flat = int;

  struct Geometrical_Parameters
  {
    double Delta;
    // optional<std::string> topology_Type;
    double Topology_LateralLength;

    std::variant<Topology_File, Topology_Flat, Topology_RMG> topology;
  };

  struct Result_Description
  {
    double ExpectedPressure;
    double ExpectedPressureTolerance;
    double ExpectedEffectiveContactAreaFraction;
    double ExpectedEffectiveContactAreaFractionTolerance;
  };



  /**
   * @brief This struct is constructed from a set of initial input parameters. It then stores or
   * derives and stores those parameters which are needed by the solver.
   *
   */
  struct InputParameters
  {
    /**
     * @brief Constructor which sets the necessary member variable parameters from an input (.xml)
     * file and creates the topology
     *
     * @param inputFileName Input file w.r.t. the calling directory
     */
    InputParameters(const std::string& inputFileName);

    /**
     * @brief Constructor which sets the necessary member variable parameters without an input
     * (.xml) file and creates the topology using the random midpoint generator
     *
     * @param E1 Young's modulus of body 1
     * @param E2 Young's modulus of body 2
     * @param nu1 Poisson's ratio of body 1
     * @param nu2 Poisson's ratio of body 2
     * @param Tolerance Tolerance for the convergence of force.
     * @param Delta Far-field displacement (Gap).
     * @param LateralLength Lateral side of the surface [micrometers]
     * @param Resolution Resolution parameter
     * @param InitialTopologyStdDeviation Initial Standard deviation for the random-midpoint
     * generator [micrometers]
     * @param Hurst Hurst Exponent (Used in random mid-point generator)
     * @param RandomSeedFlag Set `true` to fix the seed to generate psuedo random topology to
     * reproduce results. Set `false` to use random seed.
     * @param RandomGeneratorSeed Set the value of seed for the random mid-point generator
     * @param MaxIteration Maximum number of iterations for the force to converge.
     * @param WarmStartingFlag Set `true` for using the warm starter. It predicts the nodes coming
     * into contact in the next iteration and hence speeds up the computation.
     * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
     * point force
     */
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        double LateralLength, int Resolution, double InitialTopologyStdDeviation, double Hurst,
        std::optional<int> RandomGeneratorSeed, int MaxIteration, bool WarmStartingFlag,
        bool PressureGreenFunFlag);

    /**
     * @brief Constructor which sets the necessary member variable parameters without an input
     * (.xml) file and creates the topology from a specified topology (.dat) file
     *
     * @param E1 Young's modulus of body 1
     * @param E2 Young's modulus of body 2
     * @param nu1 Poisson's ratio of body 1
     * @param nu2 Poisson's ratio of body 2
     * @param Tolerance Tolerance for the convergence of force.
     * @param Delta Far-field displacement (Gap).
     * @param LateralLength Lateral side of the surface [micrometers]
     * @param TopologyFilePath Path of the input file containing the topology.
     * @param MaxIteration Maximum number of iterations for the force to converge.
     * @param WarmStartingFlag Set `true` for using the warm starter. It predicts the nodes coming
     * into contact in the next iteration and hence speeds up the computation.
     * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
     * point force
     */
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
        bool WarmStartingFlag, bool PressureGreenFunFlag);

    /**
     * @brief Constructor which sets the necessary member variable parameters without an input
     * (.xml) file and creates the topology from a specified topology (.dat) file
     *
     */
    InputParameters(Solver_Parameters& SolverParameters, Material_Parameters& MaterialParameters,
        Geometrical_Parameters& GeometricalParameters,
        optional<Result_Description&> ResultDescription = std::nullopt,
        std::optional<ViewMatrix_d&> Topology = std::nullopt);

    optional<double> elastic_compliance_correction, shape_factor, tolerance;
    optional<int> max_iterations;
    optional<bool> warm_starting_flag;
    double composite_youngs, delta, lateral_length, grid_size;
    bool pressure_green_funct_flag = false;
    ViewMatrix_d topology;
    std::optional<Result_Description> result_description;
  };
}  // namespace MIRCO::Input

#endif  // SRC_INPUTPARAMETERS_H_
