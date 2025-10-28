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

#include <myUtils/Timers.hpp>

using namespace MIRCO;

// free helpers
#include <limits>
double maxOfMat(const ViewMatrix_h v)
{
  double max = std::numeric_limits<double>::lowest();
  FOR(i, v.extent(0))
  FOR(j, v.extent(1))
  {
    if (v(i, j) > max) max = v(i, j);
  }
  return max;
}
double minOfMat(const ViewMatrix_h v)
{
  double min = std::numeric_limits<double>::max();
  FOR(i, v.extent(0))
  FOR(j, v.extent(1))
  {
    if (v(i, j) < min) min = v(i, j);
  }
  return min;
}
bool writeBitmapGrayscale(
    const ViewMatrix_h a, double vmin, double vmax, const std::string& outFile)
{
  using u8 = uint8_t;
  using u32 = uint32_t;
  const u32 h = a.extent(0), w = a.extent(1);
  if (!w || !h) return false;
  const u32 row = 3 * w, pad = (4 - (row % 4)) % 4;
  const u32 pix = (row + pad) * h;
  const u32 fsz = 54 + pix;

  auto w16 = [](std::ofstream& os, u32 v)
  {
    u8 b[2]{u8(v), u8(v >> 8)};
    os.write((char*)b, 2);
  };
  auto w32 = [](std::ofstream& os, u32 v)
  {
    u8 b[4]{u8(v), u8(v >> 8), u8(v >> 16), u8(v >> 24)};
    os.write((char*)b, 4);
  };

  std::ofstream os(outFile, std::ios::binary);
  if (!os) return false;

  // BITMAPFILEHEADER (14)
  os.put('B');
  os.put('M');
  w32(os, fsz);
  w16(os, 0);
  w16(os, 0);
  w32(os, 54);
  // BITMAPINFOHEADER (40)
  w32(os, 40);
  w32(os, w);
  w32(os, h);
  w16(os, 1);
  w16(os, 24);
  w32(os, 0);
  w32(os, pix);
  w32(os, 2835);
  w32(os, 2835);
  w32(os, 0);
  w32(os, 0);

  const double den = (vmax > vmin) ? (vmax - vmin) : 1.0;
  char zpad[3] = {0, 0, 0};
  for (u32 y = 0; y < h; ++y)
  {
    u32 r = h - 1 - y;  // bottom-up
    for (u32 x = 0; x < w; ++x)
    {
      double t = (a(r, x) - vmin) / den;
      if (t < 0) t = 0;
      if (t > 1) t = 1;
      u8 g = (u8)std::lround(t * 255.0);
      os.put((char)g);
      os.put((char)g);
      os.put((char)g);  // B,G,R
    }
    if (pad) os.write(zpad, pad);
  }
  return (bool)os;
}

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

    if (argc > 11)
    {
      std::cout << "E " << argc << "\n";
      return 1;
    }

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

    int RandomSeed =
        285928127;  // # I suppose we should keep this the same to not introduce randomness

    const auto start = std::chrono::high_resolution_clock::now();

    // IO up until here

    auto& globalR = MyUtils::Timers::TimerRegistry::globalInstance();
    globalR.start();  ////{

    InputParameters inputParams(compositeYoungs, tol, delta, lateralLength, res, stdDev, hurst,
        false, RandomSeed, 100, warmst, greenf);

    ViewVector_d meshgrid = CreateMeshgrid(inputParams.N, inputParams.grid_size);
    const auto maxAndMean = ComputeMaxAndMean(inputParams.topology);

    // Main evaluation agorithm
    double meanPressure, effectiveContactAreaFraction;
    Evaluate(meanPressure, effectiveContactAreaFraction, inputParams, maxAndMean.max, meshgrid);

    // write deformed topology and pressure //#{
    // #}

    ////}

    std::string sTim = globalR.timingReportStr();  // we do this before other stuff
    // Then do
    // MyUtils::Timers::ScopedTimer sTimer0("TIMERNAME()")
    // and remember to
    // #include <myUtils/Timers.hpp>

    std::string sTot = "";
    sTot += "\n__[[/]]\n";

    sTot += "__[[inputs]]\n";
    sTot += "compositeYoungs = " + std::to_string(compositeYoungs) + "\n";
    sTot += "tol = " + std::to_string(tol) + "\n";
    sTot += "delta = " + std::to_string(delta) + "\n";
    sTot += "lateralLength = " + std::to_string(lateralLength) + "\n";
    sTot += "res = " + std::to_string(res) + "\n";
    sTot += "stdDev = " + std::to_string(stdDev) + "\n";
    sTot += "hurst = " + std::to_string(hurst) + "\n";
    sTot += "warmstartFlag = " + (warmst ? std::string("true") : std::string("false")) + "\n";
    sTot += "greenFunctionFlag = " + (greenf ? std::string("true") : std::string("false")) + "\n";

    sTot += "__[[alwaysConstInputs]]\n";
    sTot += "RandomSeed = " + std::to_string(RandomSeed) + "\n";

    sTot += "__[[timers]]\n";
    sTot += sTim;

    sTot += "__[[outputs]]\n";
    sTot += "meanPressure = " + std::to_string(meanPressure) + "\n";
    sTot += "effectiveContactAreaFraction = " + std::to_string(effectiveContactAreaFraction) + "\n";

    std::ofstream fOut(outFile, std::ios::app);  // we append
    fOut << sTot;



    std::cout << "\n";
    if (false)
    {
      double meanPressure2, effectiveContactAreaFraction2;
      ViewVectorInt_d activeSetf2;
      ViewVector_d pf2;
      ViewVector_d w2;
      /*
            EvaluateRet(meanPressure2, effectiveContactAreaFraction2, inputParams.delta,
         inputParams.lateral_length, inputParams.grid_size, inputParams.tolerance,
         inputParams.max_iteration, inputParams.composite_youngs, inputParams.warm_starting_flag,
              inputParams.elastic_compliance_correction, inputParams.topology, maxAndMean.max,
         meshgrid, inputParams.pressure_green_funct_flag,

              activeSetf2, pf2);//, w2);
      */

      std::cout << "Writing initial topology to file\n";

      ViewMatrix_h top_h =
          Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), inputParams.topology);
      double maxx = maxOfMat(top_h);
      double minn = minOfMat(top_h);
      std::cout << "maxx=" << maxx << ", minn=" << minn << "\n";
      writeBitmapGrayscale(top_h, minn, maxx, "outFile_top_h.bmp");


      std::cout << "Writing deformed topology and pressure to files\n";
      // auto w_h = Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), w2);
      auto pf2_h = Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), pf2);
      auto activeSetf2_h =
          Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), activeSetf2);
    }
    std::cout << "\n";


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
