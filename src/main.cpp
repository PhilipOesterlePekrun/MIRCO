#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <iostream>
#include <string>

// #include "mirco_evaluate.h"
#include "mirco_evaluate_kokkos.h"
#include "mirco_inputparameters.h"
#include "mirco_topologyutilities.h"

// tmp
#include <omp.h>

#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

int main(int argc, char* argv[])
{
#if (kokkosElseOpenMP)
  Kokkos::initialize(argc, argv);
  {
    int threads_in_use = ExecSpace_Default_t::concurrency();
    std::cout << "\nthreads_in_use=" << threads_in_use << "\n";
    std::cout << "\nExecSpace_Default= " << typeid(ExecSpace_Default_t).name() << "\n";
    std::cout << "\nExecSpace_DefaultHost= " << typeid(ExecSpace_DefaultHost_t).name() << "\n\n";
    std::cout << "\nMemorySpace_Host_t= " << typeid(MemorySpace_Host_t).name() << "\n";
    std::cout << "\nMemorySpace_ofDefaultExec_t= " << typeid(MemorySpace_ofDefaultExec_t).name()
              << "\n\n";
#else
  int max_threads = omp_get_max_threads();
  std::cout << "OPENMP omp_get_max_threads=" << max_threads << "\n\n";
#endif
    TimerRegistry::globalInstance().start();

    TEUCHOS_TEST_FOR_EXCEPTION(
        argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
    // reading the input file name from the command line
    std::string inputFileName = argv[1];

    const auto start = std::chrono::high_resolution_clock::now();

    MIRCO::InputParameters inputParams(inputFileName);

    // Identical Vectors/Matricies, therefore only created one here.
    ViewVector_d meshgrid = MIRCO::CreateMeshgrid(inputParams.N_, inputParams.grid_size_);

    auto& topology = *(inputParams.topology_);
    auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

    // Main evaluation agorithm
    double meanPressure = MIRCO::Evaluate(inputParams, max_and_mean.max_, meshgrid);

    std::cout << "Mean pressure is: " << std::to_string(meanPressure) << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;

    std::cout << TimerRegistry::globalInstance().timingReportStr(true);
#if (kokkosElseOpenMP)
  }
  Kokkos::finalize();
#endif
}
