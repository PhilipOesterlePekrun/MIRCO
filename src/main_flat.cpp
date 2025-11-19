#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#include "mirco_kokkostypes.h"
#include "mirco_topologyutilities.h"
#include "mirco_contactstatus.h"
#include "mirco_matrixsetup.h"

#include <KokkosLapack_gesv.hpp>

using namespace MIRCO;

// flat and args
int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::string argv1 = std::string(argv[1]);
    if (argv1 == "-help" || argv1 == "-h" || argv1 == "--help" || argv1 == "--h")
    {
      std::cout
          << "TODO description\n"
          << "\tInputs:"
          << "\t\tN: Lateral mesh grid count (positive integer); if prefixed with a (i.e. aN, where "
            "N is a positive integer), all N from 0 to N (inclusive) will be sequentially executed"
          << "\t\tDelta (far-field displacement), "
            "Composite Young's Modulus, LateralLength, Pressure Green Function Flag\n"
          << "\tOutput: TODO\n"
          << std::endl;
      return 0;
    }
    if (argc != 6) throw std::runtime_error("The code expects 5 arguments. Use --help for detials.");

    std::cout << "-- Kokkos information --\n";
    std::cout << "Threads in use: " << ExecSpace_Default_t().concurrency() << "\n";
    std::cout << "Default execution space: " << typeid(ExecSpace_Default_t).name() << "\n";
    std::cout << "Default host execution space: " << typeid(ExecSpace_DefaultHost_t).name() << "\n";
    std::cout << "Default memory space: " << typeid(MemorySpace_ofDefaultExec_t).name() << "\n";
    std::cout << "Default host memory space: " << typeid(MemorySpace_Host_t).name() << "\n";
    std::cout << "\n";

    bool computeAllUpToN = (argv[1][0] == 'a');
    if (computeAllUpToN) argv1 = argv1.substr(1, argv1.length());
    int nominalN = std::stoi(argv1);
    double Delta = std::stod(argv[2]);
    double CompositeYoungs = std::stod(argv[3]);
    double LateralLength = std::stod(argv[4]);
    std::string argv5 = std::string(argv[5]);
    bool PressureGreenFunFlag =
        (argv5 == "t" || argv5 == "T" || argv5 == "true" || argv5 == "True" || argv5 == "1");

    const auto start = std::chrono::high_resolution_clock::now();

    if (computeAllUpToN)
      std::cout << "Shape factors (PressureGreenFunFlag=" << (PressureGreenFunFlag ? "true" : "false")
                << "):\n{\n";
    int N;
    for (N = (computeAllUpToN ? 0 : nominalN); N <= nominalN; ++N)
    {
      try
      {
        if (N <= 0)
        {
          if (N != 0) std::cout << "," << std::endl;
          std::cout << std::fixed << std::setprecision(16) << -1;
          continue;
        }
        const double gridSize = LateralLength / N;

        ViewVector_d meshgrid = CreateMeshgrid(N, gridSize);

        const int N2 = N * N;

        ViewVector_d xv0 = ViewVector_d("xv0", N2);
        ViewVector_d yv0 = ViewVector_d("yv0", N2);
        ViewScalarInt_d counter("counter");
        Kokkos::deep_copy(counter, 0);
        Kokkos::parallel_for(
            N2, KOKKOS_LAMBDA(const int a) {
              const int i = a / N;
              const int j = a % N;
              const int aa = Kokkos::atomic_fetch_add(&counter(), 1);
              xv0(aa) = meshgrid(i);
              yv0(aa) = meshgrid(j);
            });

        // For a flat indentor, the following hold: displacement field = const, active set = entire
        // domain
        ViewMatrix_d H = SetupMatrix(
            xv0, yv0, gridSize, CompositeYoungs, N2, pressureGreenFunFlag);
        ViewVector_d b0p("b0p", N2);
        Kokkos::deep_copy(b0p, Delta);
        ViewVectorInt_d ipiv("ipiv", N2);
        // Solve H s = b0; b0p becomes s
        KokkosLapack::gesv(H, b0p, ipiv);

        double totalForce;
        double contactArea;
        ComputeContactForceAndArea(
            totalForce, contactArea, b0p, gridSize, LateralLength, PressureGreenFunFlag);

        double domainArea = LateralLength * LateralLength;

        double meanPressure = totalForce / domainArea;
        double effectiveContactAreaFraction = contactArea / domainArea;

        // w_{el} = Delta = meanPressure * l * \alpha / CompositeYoungs
        double shapeFactor_alpha = Delta * CompositeYoungs * LateralLength / (totalForce);

        if (computeAllUpToN)
        {
          if (N != 0) std::cout << "," << std::endl;
          std::cout << std::fixed << std::setprecision(16) << shapeFactor_alpha;
        }
        else
        {
          std::cout << std::fixed << std::setprecision(16) << "Mean pressure is: " << meanPressure
                    << "\nEffective contact area fraction (should be 1 for a flat topology) is: "
                    << effectiveContactAreaFraction << std::endl;

          std::cout << std::fixed << std::setprecision(16)
                    << "Calculated shape factor is: " << shapeFactor_alpha;
        }
      }
      // catch (const Kokkos::Experimental::RawMemoryAllocationFailure&) {break;}
      catch (const std::bad_alloc&)
      {
        break;
      }
      catch (const std::exception&)
      {
        break;
      }
    }
    
    const auto finish = std::chrono::high_resolution_clock::now();

    if (computeAllUpToN) std::cout << "\n}\n\nComputed shape factors up to N=" << N - 1 << " ";
    if (N - 1 != nominalN) std::cerr << "(ran out of memory)";
    
    std::cout << "\n\nElapsed time is: " + std::to_string(std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count()) + "s"<<std::endl;
    
    if (N - 1 != nominalN) return 1;
  }
  Kokkos::finalize();
}
