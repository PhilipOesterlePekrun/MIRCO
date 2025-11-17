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

    const auto start = std::chrono::high_resolution_clock::now();

    int res = std::stoi(argv[1]);
    double Delta = std::stod(argv[2]);
    double LateralLength = std::stod(argv[3]);
    double Tol = std::stod(argv[4]);
    int maxIter = std::stoi(argv[5]);
    // double Tol = std::stod(argv[4]);

    InputParameters inputParams(1.0, 1.0, 0.3, 0.3, Tol, Delta, LateralLength, res, 5.0, 0.7, false,
        1201, maxIter, false, true);

    ViewVector_d meshgrid = CreateMeshgrid(inputParams.N, inputParams.grid_size);
    const auto maxAndMean = ComputeMaxAndMean(inputParams.topology);

    // Main evaluation agorithm
    double meanPressure, effectiveContactAreaFraction;
    Evaluate(meanPressure, effectiveContactAreaFraction, inputParams, maxAndMean.max, meshgrid);

    const auto finish = std::chrono::high_resolution_clock::now();

    std::cout << std::setprecision(16) << "Mean pressure is: " << meanPressure
              << "\nEffective contact area fraction is: " << effectiveContactAreaFraction
              << std::endl;

    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s" << std::endl;
  }
  Kokkos::finalize();
}
