#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <omp.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"
#include "mirco_topologyutilities.h"

#include <myUtils/Timers.hpp>

using namespace MIRCO;

int main(int argc, char* argv[])
{
    double compositeYoungs = std::stod(argv[1]);
    double tol = std::stod(argv[2]);
    double delta = std::stod(argv[3]);
    double lateralLength = std::stod(argv[4]);
    int res = std::stoi(argv[5]);
    double stdDev = std::stod(argv[6]);
    double hurst = std::stod(argv[7]);
    bool warmst = argv[8][0] == 't';
    bool greenf = argv[9][0] == 't';

    std::string outFile = argv[10];
    
    if (argc != 11)
    {
      std::cout << "E: Incorrect arg count: " << argc << "\n";
      return 1;
    }
    
    std::cout << "omp_get_max_threads() = " << omp_get_max_threads() << "\n";
    

    int RandomSeed =
        285928127;  // # I suppose we should keep this the same to not introduce randomness

    const auto start = std::chrono::high_resolution_clock::now();

    // IO up until here

    auto& globalR = MyUtils::Timers::TimerRegistry::globalInstance();
    globalR.start();  ////{

/*KOKKOS
InputParameters inputParams(compositeYoungs, tol, delta, lateralLength, res, stdDev, hurst,
    false, RandomSeed, 100, warmst, greenf);

ViewVector_d meshgrid = CreateMeshgrid(inputParams.N, inputParams.grid_size);
const auto maxAndMean = ComputeMaxAndMean(inputParams.topology);

// Main evaluation agorithm
double meanPressure, effectiveContactAreaFraction;
Evaluate(meanPressure, effectiveContactAreaFraction, inputParams, maxAndMean.max, meshgrid);
*/
    
    
    
    InputParameters inputParams(compositeYoungs, tol, delta, lateralLength, res, stdDev, hurst,
    false, RandomSeed, 1000, warmst, greenf);

  // Identical Vectors/Matricies, therefore only created one here.
  auto meshgrid = MIRCO::CreateMeshgrid(inputParams.N, inputParams.grid_size);

  auto& topology = *(inputParams.topology);
  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Initialise Pressure
  double meanPressure = 0.0;
  MIRCO::Evaluate(inputParams, meanPressure, max_and_mean.max_, meshgrid);

    ////}

    std::string sTim = globalR.timingReportStr();  // we do this before other stuff
    // Then do
    // MyUtils::Timers::ScopedTimer sTimer0("TIMERNAME()")
    // and remember to
    // #include <myUtils/Timers.hpp>

    std::string sTot = "";
    sTot += "\n__[[/]]\n";
    sTot += "Original Implementation (Teuchos, raw OpenMP)\n";
    sTot += "__[[inputs]]\n";
    sTot += "numThreads = " + std::to_string((int)omp_get_max_threads()) + "\n";
    sTot += "compositeYoungs = " + std::to_string(compositeYoungs) + "\n";
    sTot += "tol = " + std::to_string(tol) + "\n";
    sTot += "delta = " + std::to_string(delta) + "\n";
    sTot += "lateralLength = " + std::to_string(lateralLength) + "\n";
    sTot += "res = " + std::to_string(res) + "\n";
    sTot += "stdDev = " + std::to_string(stdDev) + "\n";
    sTot += "hurst = " + std::to_string(hurst) + "\n";
    sTot += "warmstartFlag = " + (warmst ? std::string("true") : std::string("false")) + "\n";
    sTot += "greenFunctionFlag = " + (greenf ? std::string("true") : std::string("false")) + "\n";
    sTot += "RandomSeed = " + std::to_string(RandomSeed) + "\n";

    sTot += "__[[timers]]\n";
    sTot += sTim;

    sTot += "__[[outputs]]\n";
    sTot += "meanPressure = " + std::to_string(meanPressure) + "\n";

    std::ofstream fOut(outFile, std::ios::app);  // we append
    fOut << sTot;




    std::cout << std::setprecision(16) << "Mean pressure is: " << meanPressure
              << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s" << std::endl;
}
