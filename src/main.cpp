#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"
#include "mirco_kokkostypes.h"
#include "mirco_topologyutilities.h"
#include "mirco_utils.h"

using namespace MIRCO;

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::cout << "-- Kokkos information --\n";
    std::cout << "Threads in use: " << ExecSpace_Default_t().concurrency() << "\n";
    std::cout << "Default execution space: " << typeid(ExecSpace_Default_t).name() << "\n";
    std::cout << "Default host execution space: " << typeid(ExecSpace_DefaultHost_t).name() << "\n";
    std::cout << "Default memory space: " << typeid(MemorySpace_ofDefaultExec_t).name() << "\n";
    std::cout << "Default host memory space: " << typeid(MemorySpace_Host_t).name() << "\n";
    std::cout << "\n";

    std::string inputFileName = argv[1];

    int numSamples = (argc == 3) ? std::stoi(argv[2]) : 1000;

    int firstSeed = 17283049;

    //{
    std::vector<double> meanPressureVect;
    std::vector<double> contactAreaFractionVect;
    std::vector<double> elapsedTimeVect;
    //}


    const auto startGlobal = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numSamples; ++i)
    {
      std::cout << "i=" << i << "\n";

      const auto start = std::chrono::high_resolution_clock::now();

      std::ifstream fin(inputFileName);
      if (!fin) throw std::runtime_error("Cannot open file: " + inputFileName);

      std::stringstream ss;
      ss << fin.rdbuf();
      std::string inString = ss.str();

      ryml::Tree tree = ryml::parse_in_arena(c4::to_csubstr(inString));
      ryml::ConstNodeRef root = tree["mirco_input"];
      ryml::ConstNodeRef parameters = root["parameters"];
      ryml::ConstNodeRef geoParams = parameters["geometrical_parameters"];
      ryml::ConstNodeRef matParams = parameters["material_parameters"];

      if (root.invalid()) throw std::runtime_error("Input incomplete: missing root `mirco_input`");
      if (parameters.invalid())
        throw std::runtime_error("Input incomplete: missing section `parameters`");
      if (geoParams.invalid())
        throw std::runtime_error("Input incomplete: missing section `geometrical_parameters`");
      if (matParams.invalid())
        throw std::runtime_error("Input incomplete: missing section `material_parameters`");

      auto exportVisualization = Utils::get_optional_bool(root, "ExportVisualization");
      std::optional<std::string> exportVisualizationPath;
      if (exportVisualization && exportVisualization.value())
        exportVisualizationPath = Utils::get_string(root, "ExportVisualizationPath");
      else
        exportVisualizationPath = std::nullopt;

      if (!Utils::get_bool(root, "RandomTopologyFlag"))
        throw std::runtime_error("Needs to be randomtopology");

      std::optional<int> randomSeed = firstSeed + i;

      // std::cout<<"\trandomSeed="<<*randomSeed<<"\n";

      InputParameters trueInputParams = InputParameters(Utils::get_double(matParams, "E1"),
          Utils::get_double(matParams, "E2"), Utils::get_double(matParams, "nu1"),
          Utils::get_double(matParams, "nu2"), Utils::get_double(geoParams, "Tolerance"),
          Utils::get_double(geoParams, "Delta"), Utils::get_double(geoParams, "LateralLength"),
          Utils::get_int(geoParams, "Resolution"),
          Utils::get_double(geoParams, "InitialTopologyStdDeviation"),
          Utils::get_double(geoParams, "HurstExponent"), Utils::get_int(root, "MaxIteration"),
          Utils::get_bool(root, "WarmStartingFlag"), Utils::get_bool(root, "PressureGreenFunFlag"),
          Utils::get_bool(root, "RandomSeedFlag"), randomSeed, exportVisualizationPath);

      int thisN = trueInputParams.topology.extent(0);
      if (i == 0)
      {
        std::cout << "thisN=" << thisN << "\n";
        std::cout << "//#\ntopology=\n";
        for (int ii = 0; ii < thisN; ++ii)
        {
          for (int jj = 0; jj < thisN; ++jj)
          {
            std::cout << trueInputParams.topology(ii, jj) << " ";
          }
          std::cout << "\n";
        }
      }

      ViewVector_d meshgrid = CreateMeshgrid(trueInputParams.N, trueInputParams.grid_size);
      const double topologyMax = GetMax(trueInputParams.topology);

      double meanPressure, effectiveContactAreaFraction;
      Evaluate(meanPressure, effectiveContactAreaFraction, trueInputParams, topologyMax, meshgrid);

      const auto finish = std::chrono::high_resolution_clock::now();

      const double elapsedTime =
          std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();

      meanPressureVect.push_back(meanPressure);
      contactAreaFractionVect.push_back(effectiveContactAreaFraction);
      elapsedTimeVect.push_back(elapsedTime);
    }



    const auto finishGlobal = std::chrono::high_resolution_clock::now();

    const double elapsedTimeGlobal =
        std::chrono::duration_cast<std::chrono::duration<double>>(finishGlobal - startGlobal)
            .count();

    //{
    std::cout << "meanPressureVect=[\n";
    for (int i = 0; i < numSamples; ++i)
    {
      std::cout << meanPressureVect[i] << "\n";
    }
    std::cout << "]\n\n";

    std::cout << "contactAreaFractionVect=[\n";
    for (int i = 0; i < numSamples; ++i)
    {
      std::cout << contactAreaFractionVect[i] << "\n";
    }
    std::cout << "]\n\n";

    std::cout << "elapsedTimeVect=[\n";
    for (int i = 0; i < numSamples; ++i)
    {
      std::cout << elapsedTimeVect[i] << "\n";
    }
    std::cout << "]\n\n";


    double mean_meanPressure = 0.0;
    double mean_contactAreaFraction = 0.0;
    double mean_elapsedTime = 0.0;

    for (int i = 0; i < numSamples; ++i)
    {
      mean_meanPressure += meanPressureVect[i];
      mean_contactAreaFraction += contactAreaFractionVect[i];
      mean_elapsedTime += elapsedTimeVect[i];
    }
    mean_meanPressure /= numSamples;
    mean_contactAreaFraction /= numSamples;
    mean_elapsedTime /= numSamples;

    std::cout << "mean_meanPressure = " << mean_meanPressure << "\n";
    std::cout << "mean_contactAreaFraction = " << mean_contactAreaFraction << "\n";
    std::cout << "mean_elapsedTime = " << mean_elapsedTime << "\n\n";

    double variance_meanPressure = 0.0;
    double variance_contactAreaFraction = 0.0;
    double variance_elapsedTime = 0.0;

    for (int i = 0; i < numSamples; ++i)
    {
      variance_meanPressure += pow((meanPressureVect[i] - mean_meanPressure), 2);
      variance_contactAreaFraction +=
          pow((contactAreaFractionVect[i] - mean_contactAreaFraction), 2);
      variance_elapsedTime += pow((elapsedTimeVect[i] - mean_elapsedTime), 2);
    }
    variance_meanPressure /= (numSamples - 1);
    variance_contactAreaFraction /= (numSamples - 1);
    variance_elapsedTime /= (numSamples - 1);

    std::cout << "variance_meanPressure = " << variance_meanPressure << "\n";
    std::cout << "variance_contactAreaFraction = " << variance_contactAreaFraction << "\n";
    std::cout << "variance_elapsedTime = " << variance_elapsedTime << "\n\n";

    std::cout << "stddev meanPressure = " << sqrt(variance_meanPressure) << "\n";
    std::cout << "stddev contactAreaFraction = " << sqrt(variance_contactAreaFraction) << "\n";
    std::cout << "stddev elapsedTime = " << sqrt(variance_elapsedTime) << "\n\n";

    std::cout << "coeff of variation meanPressure = "
              << sqrt(variance_meanPressure) / mean_meanPressure << "\n";
    std::cout << "coeff of variation contactAreaFraction = "
              << sqrt(variance_contactAreaFraction) / mean_contactAreaFraction << "\n";
    std::cout << "coeff of variation elapsedTime = "
              << sqrt(variance_elapsedTime) / mean_elapsedTime << "\n\n";

    std::cout << "elapsedTimeGlobal = " + std::to_string(elapsedTimeGlobal) + "s\n";
    //}
  }
  Kokkos::finalize();
}
