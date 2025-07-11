#include <omp.h>

#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include "Timer.hpp"
#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"
#include "mirco_topologyutilities.h"

inline std::string levelizeString(const std::string& s)
{
  std::string lvlStr = "\t";

  std::string tmpS = lvlStr;
  for (int i = 0; i < s.length(); ++i)
  {
    tmpS += s[i];
    if (s[i] == '\n' && i != s.length() - 1)  // -1 because we don't want to turn the last "\n" into
                                              // "\n\t", or else the next thing is affected
      tmpS += lvlStr;
  }
  return tmpS;
}



int main(int argc, char* argv[])
{
  TimerRegistry::globalInstance().start();
  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  /// std::string inputFileName = argv[1];

  const auto start = std::chrono::high_resolution_clock::now();

  // MIRCO::InputParameters inputParams(inputFileName);
  //  argv: ./mirco resolution delta tolerance
  MIRCO::InputParameters inputParams(1.0, 1.0, 0.3, 0.3, std::stod(argv[3]), std::stod(argv[2]),
      1000.0, std::stoi(argv[1]), 20.0, 0.7, false, 95, 100, true, true);

  // Identical Vectors/Matricies, therefore only created one here.
  auto meshgrid = MIRCO::CreateMeshgrid(inputParams.N_, inputParams.grid_size_);

  auto& topology = *(inputParams.topology_);
  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Initialise Pressure
  double pressure = 0.0;
  MIRCO::Evaluate(inputParams, pressure, max_and_mean.max_, meshgrid);

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsedTime = finish - start;
  std::cout << "Elapsed time is: " << elapsedTime.count() << "s." << std::endl;
  // std::cout<<TimerRegistry::globalInstance().timingReportStr(true);


  linearSystemSizeAvg = (double)linearSystemSizeSummed / numLinearCallsTotal;



  std::string timerFilePath =
      "/home/oesterle/rd/MIRCO_Base/TimingStudyBeforeKokkos/" + std::string(argv[4]);
  std::fstream outFile(timerFilePath, std::ios::app);
  outFile << "--delimiter--" << "\n\n";

  outFile << "\t-numThreads=" << omp_get_max_threads() << "\n\n";

  outFile << "\t-meanPressure=" << std::to_string(pressure) << "\n";
  outFile << "\t-totalTime=" << elapsedTime.count() << "s." << "\n\n";

  outFile << "\t-resolution=" << argv[1] << "\n";
  outFile << "\t-delta=" << argv[2] << "\n";
  outFile << "\t-tolerance=" << argv[3] << "\n\n";

  outFile << "\t-numNonlinearIters=" << numNonlinearIters << "\n";
  outFile << "\t-numLinearCallsTotal=" << numLinearCallsTotal << "\n";
  outFile << "\t-linearSystemSizeSummed=" << linearSystemSizeSummed << "\n";
  outFile << "\t-linearSystemSizeAvg=" << linearSystemSizeAvg << "\n";
  outFile << "\t-evaluateIterations=" << evaluateIterations << "\n\n";

  outFile << levelizeString(TimerRegistry::globalInstance().timingReportStr(true)) << "\n\n\n";
}
