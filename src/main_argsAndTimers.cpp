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

#define myUtils_ENABLE_TIMERS true
#include <myUtils/Timers.hpp>

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

    double E1 = std::stod(argv[1]);
    double E2 = std::stod(argv[2]);
    double nu1 = std::stod(argv[3]);
    double nu2 = std::stod(argv[4]);

    double tol = std::stod(argv[5]);
    double delta = std::stod(argv[6]);
    int res = std::stoi(argv[7]);
    double stdDev = std::stod(argv[8]);
    double hurst = std::stod(argv[9]);
    double warmst = std::stoi(argv[10]) > 0;
    double greenf = std::stoi(argv[11]) > 0;

    const auto start = std::chrono::high_resolution_clock::now();

    // IO up until here

    auto& globalR = MyUtils::Timers::TimerRegistry::globalInstance();
    globalR.start();  ////{

    InputParameters inputParams(
        E1, E2, nu1, nu2, tol, delta, 1000, res, stdDev, hurst, false, -1, 100, warmst, greenf);

    ViewVector_d meshgrid = CreateMeshgrid(inputParams.N, inputParams.grid_size);
    const auto maxAndMean = ComputeMaxAndMean(inputParams.topology);

    // Main evaluation agorithm
    double meanPressure, effectiveContactAreaFraction;
    Evaluate(meanPressure, effectiveContactAreaFraction, inputParams, maxAndMean.max, meshgrid);

    ////}
    auto s = globalR.timingReportStr();

    std::cout << "\n" << s << "\n";

    std::cout << std::setprecision(16) << "Mean pressure is: " << meanPressure
              << "\nEffective contact area fraction is: " << effectiveContactAreaFraction
              << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s" << std::endl;
  }
  Kokkos::finalize();
}
