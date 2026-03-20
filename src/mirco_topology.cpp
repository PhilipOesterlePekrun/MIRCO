#include "mirco_topology.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>

namespace MIRCO
{
  ViewMatrix_h CreateSurfaceFromFile(const std::string& filepath)
  {
    int N = 0;

    std::ifstream reader(filepath);
    std::string line;

    while (getline(reader, line))
    {
      ++N;
    }
    reader.clear();
    // Reuse reader
    reader.seekg(0, std::ios::beg);

    ViewMatrix_h z("CreateSurfaceFromFile(); z", N, N);
    std::ifstream stream(filepath);
    int lineCounter = 0;
    while (getline(reader, line))
    {
      ++lineCounter;
      for (int i = 0; i < N; i++)
      {
        const int separatorPosition = line.find_first_of(';');
        z(lineCounter - 1, i) = stod(line.substr(0, separatorPosition));
        line = line.substr(separatorPosition + 1, line.length());
      }
    }
    reader.close();

    return z;
  }



#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

  // Assumes ViewMatrix_h supports:
  //   ViewMatrix_h(name, rows, cols);
  //   z(i,j) get/set
  //   z.extent(0), z.extent(1) if you need them (not used below)

  static bool IsPowerOfTwo(int x) { return x > 0 && (x & (x - 1)) == 0; }

  static int ILog2Pow2(int x)
  {
    // x must be power of two
    int r = 0;
    while ((1 << r) < x) ++r;
    return r;
  }

  static std::uint32_t MixSeed(std::uint32_t base, std::uint32_t idx)
  {
    // Simple, decent mixing for seeds (deterministic).
    std::uint32_t v = base ^ (0x9E3779B9u + idx * 0x85EBCA6Bu);
    v ^= v >> 16;
    v *= 0x7FEB352Du;
    v ^= v >> 15;
    v *= 0x846CA68Bu;
    v ^= v >> 16;
    return v;
  }

  static double MeanOfBlock(const ViewMatrix_h& z, int i0, int j0, int ni, int nj)
  {
    double sum = 0.0;
    for (int i = 0; i < ni; ++i)
    {
      for (int j = 0; j < nj; ++j)
      {
        sum += z(i0 + i, j0 + j);
      }
    }
    return sum / static_cast<double>(ni * nj);
  }

  static double MinOfMatrix(const ViewMatrix_h& z, int N)
  {
    double mn = std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        mn = std::min(mn, z(i, j));
      }
    }
    return mn;
  }

  // Fits z ~= a*x + b*y + c on a cell-center grid over [0, L]x[0, L]
  // using the top-left ni x nj block of patchZ (starting at 0,0).
  // Returns (a,b,c) in plane[0..2].
  static void FitPlaneLeastSquares(
      const ViewMatrix_h& patchZ, int ni, int nj, double L, double plane[3])
  {
    const double dx = L / static_cast<double>(nj);
    const double dy = L / static_cast<double>(ni);

    // Build normal equations:
    // [Sxx Sxy Sx] [a] = [Sxz]
    // [Sxy Syy Sy] [b]   [Syz]
    // [Sx  Sy  n ] [c]   [Sz ]
    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0;
    double Sx = 0.0, Sy = 0.0;
    double Sxz = 0.0, Syz = 0.0, Sz = 0.0;
    const double n = static_cast<double>(ni * nj);

    for (int i = 0; i < ni; ++i)
    {
      const double y = (static_cast<double>(i) + 0.5) * dy;
      for (int j = 0; j < nj; ++j)
      {
        const double x = (static_cast<double>(j) + 0.5) * dx;
        const double z = patchZ(i, j);

        Sxx += x * x;
        Sxy += x * y;
        Syy += y * y;
        Sx += x;
        Sy += y;
        Sxz += x * z;
        Syz += y * z;
        Sz += z;
      }
    }

    // Solve 3x3 with Gaussian elimination (small, stable enough here).
    double A[3][4] = {{Sxx, Sxy, Sx, Sxz}, {Sxy, Syy, Sy, Syz}, {Sx, Sy, n, Sz}};

    // Forward elimination
    for (int k = 0; k < 3; ++k)
    {
      // Pivot
      int piv = k;
      double best = std::abs(A[k][k]);
      for (int r = k + 1; r < 3; ++r)
      {
        const double val = std::abs(A[r][k]);
        if (val > best)
        {
          best = val;
          piv = r;
        }
      }
      if (best == 0.0) throw std::runtime_error("Plane fit failed (singular system).");

      if (piv != k)
      {
        for (int c = k; c < 4; ++c) std::swap(A[k][c], A[piv][c]);
      }

      // Eliminate
      const double akk = A[k][k];
      for (int r = k + 1; r < 3; ++r)
      {
        const double f = A[r][k] / akk;
        for (int c = k; c < 4; ++c) A[r][c] -= f * A[k][c];
      }
    }

    // Back substitution
    double xsol[3];
    for (int r = 2; r >= 0; --r)
    {
      double rhs = A[r][3];
      for (int c = r + 1; c < 3; ++c) rhs -= A[r][c] * xsol[c];
      xsol[r] = rhs / A[r][r];
    }

    plane[0] = xsol[0];  // a
    plane[1] = xsol[1];  // b
    plane[2] = xsol[2];  // c
  }

  static void SubtractPlane(ViewMatrix_h& patchZ, int ni, int nj, double L, const double plane[3])
  {
    const double dx = L / static_cast<double>(nj);
    const double dy = L / static_cast<double>(ni);

    const double a = plane[0];
    const double b = plane[1];
    const double c = plane[2];

    for (int i = 0; i < ni; ++i)
    {
      const double y = (static_cast<double>(i) + 0.5) * dy;
      for (int j = 0; j < nj; ++j)
      {
        const double x = (static_cast<double>(j) + 0.5) * dx;
        patchZ(i, j) -= (a * x + b * y + c);
      }
    }
  }



  ViewMatrix_h CreateRmgSurface(int Resolution, double InitialTopologyStdDeviation, double Hurst,
      bool RandomSeedFlag, std::optional<int> RandomGeneratorSeed)
  {
    srand(time(NULL));

    int seed;

    if (RandomSeedFlag)
      seed = rand();
    else if (RandomGeneratorSeed)
      seed = *RandomGeneratorSeed;
    else
      throw std::runtime_error(
          "Please provide 'RandomGeneratorSeed' when 'RandomSeedFlag' is false.");

    std::default_random_engine generate(seed);
    std::normal_distribution<double> distribution(
        0.0, 1.0);  // normal distribution: mean = 0.0, standard deviation = 1.0

    int N = (1 << Resolution) + 1;
    ViewMatrix_h z("CreateRmgSurface(); z", N, N);

    const double scaling_factor = pow(2.0, 0.5 * Hurst);
    double alpha = InitialTopologyStdDeviation * scaling_factor;

    const int D_0 = N - 1;
    int D = D_0;
    int d = D_0 / 2;

    for (int i = 0; i < Resolution; i++)
    {
      alpha = alpha / scaling_factor;

      for (int j = d; j < D_0 - d + 1; j = j + D)
      {
        for (int k = d; k < D_0 - d + 1; k = k + D)
        {
          z(j, k) = (z(j + d, k + d) + z(j + d, k - d) + z(j - d, k + d) + z(j - d, k - d)) / 4 +
                    alpha * distribution(generate);
        }
      }

      alpha = alpha / scaling_factor;

      for (int j = d; j < D_0 - d + 1; j = j + D)
      {
        z(j, 0) = (z(j + d, 0) + z(j - d, 0) + z(j, d)) / 3 + alpha * distribution(generate);
        z(j, D_0) =
            (z(j + d, D_0) + z(j - d, D_0) + z(j, D_0 - d)) / 3 + alpha * distribution(generate);
        z(0, j) = (z(0, j + d) + z(0, j - d) + z(d, j)) / 3 + alpha * distribution(generate);
        z(D_0, j) =
            (z(D_0, j + d) + z(D_0, j - d) + z(D_0 - d, j)) / 3 + alpha * distribution(generate);
      }

      for (int j = d; j < D_0 - d + 1; j = j + D)
      {
        for (int k = D; k < D_0 - d + 1; k = k + D)
        {
          z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                    alpha * distribution(generate);
        }
      }

      for (int j = D; j < D_0 - d + 1; j = j + D)
      {
        for (int k = d; k < D_0 - d + 1; k = k + D)
        {
          z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                    alpha * distribution(generate);
        }
      }

      D = D / 2;
      d = d / 2;
    }

    // Finding minimum of topology
    double zmin = std::numeric_limits<double>::max();
    for (int i = 0; i < D_0 + 1; i++)
    {
      for (int j = 0; j < D_0 + 1; j++)
      {
        zmin = std::min(zmin, z(i, j));
      }
    }

    // Setting the minimum of topology to zero
    for (int i = 0; i < D_0 + 1; i++)
    {
      for (int j = 0; j < D_0 + 1; j++)
      {
#if (REGULARMIRCO_ELSEDODELTARELATIVETOOTHER)
        z(i, j) = z(i, j) - zmin;
#else
        z(i, j) = z(i, j);
#endif
      }
    }

    return z;
  }



  ViewMatrix_h CreatePatchBasedRmgSurface(int Resolution, double InitialTopologyStdDeviation,
      double Hurst, bool RandomSeedFlag, std::optional<int> RandomGeneratorSeed, int NumPatches,
      double L_patch, bool RemoveSlopePerPatch, bool RemoveMeanPerPatch, double SeamBlend)
  {
    const int patches_per_side =
        static_cast<int>(std::lround(std::sqrt(static_cast<double>(NumPatches))));
    if (patches_per_side * patches_per_side != NumPatches)
      throw std::runtime_error("NumPatches must be a perfect square.");

    if (!IsPowerOfTwo(patches_per_side))
      throw std::runtime_error(
          "patches_per_side must be a power of two for exact tiling with (2^R + 1) grids.");

    const int log2_pps = ILog2Pow2(patches_per_side);
    const int localResolution = Resolution - log2_pps;
    if (localResolution < 1)
      throw std::runtime_error("Resolution too small for the requested number of patches.");

    const int N_global = (1 << Resolution) + 1;
    std::cout << "N_global=" << N_global << "\n";
    const int N_local = (1 << localResolution) + 1;
    const int n_local = N_local - 1;  // interior span (without the overlapping last row/col)

    ViewMatrix_h z_raw("CreatePatchBasedRmgSurface(); z_raw", N_global, N_global);

    // Base seed
    std::uint32_t baseSeed = 0u;
    if (RandomSeedFlag)
    {
      std::random_device rd;
      baseSeed = static_cast<std::uint32_t>(rd());
    }
    else if (RandomGeneratorSeed)
    {
      baseSeed = static_cast<std::uint32_t>(*RandomGeneratorSeed);
    }
    else
    {
      throw std::runtime_error("Provide RandomGeneratorSeed when RandomSeedFlag is false.");
    }

    // Generate and stitch patches.
    for (int pr = 0; pr < patches_per_side; ++pr)
    {
      for (int pc = 0; pc < patches_per_side; ++pc)
      {
        const int pidx = pr * patches_per_side + pc;
        const std::uint32_t seed = MixSeed(baseSeed, static_cast<std::uint32_t>(pidx));
        std::optional<int> seedO = seed;
        // Force deterministic per-patch seed; avoid RandomSeedFlag inside loop.
        ViewMatrix_h patch =
            CreateRmgSurface(localResolution, InitialTopologyStdDeviation, Hurst, false, seedO);

        // Optionally remove best-fit plane (slope) on the (N_local-1)x(N_local-1) region.
        if (RemoveSlopePerPatch)
        {
          double plane[3];
          FitPlaneLeastSquares(patch, n_local, n_local, L_patch, plane);
          SubtractPlane(patch, n_local, n_local, L_patch, plane);
        }

        // Optionally remove mean on the same region.
        if (RemoveMeanPerPatch)
        {
          const double m = MeanOfBlock(patch, 0, 0, n_local, n_local);
          for (int i = 0; i < n_local; ++i)
            for (int j = 0; j < n_local; ++j) patch(i, j) -= m;
        }

        // Copy into global, skipping overlaps except on the outer boundary.
        const int i0 = pr * n_local;
        const int j0 = pc * n_local;

        const int ni = (pr == patches_per_side - 1) ? N_local : n_local;
        const int nj = (pc == patches_per_side - 1) ? N_local : n_local;

        for (int i = 0; i < ni; ++i)
        {
          for (int j = 0; j < nj; ++j)
          {
            z_raw(i0 + i, j0 + j) = patch(i, j);
          }
        }
      }
    }

    // Smooth seams (blend rows/cols at patch borders), like the Python code.
    ViewMatrix_h z_sm("CreatePatchBasedRmgSurface(); z_sm", N_global, N_global);
    for (int i = 0; i < N_global; ++i)
      for (int j = 0; j < N_global; ++j) z_sm(i, j) = z_raw(i, j);

    const double th = SeamBlend;

    // Horizontal seams (between patch rows)
    for (int pr = 0; pr < patches_per_side; ++pr)
    {
      const int i_top = n_local * pr;
      const int i_bottom = n_local * (pr + 1) - 1;

      if (i_top > 0)
      {
        for (int j = 0; j < N_global; ++j)
          z_sm(i_top, j) = th * z_raw(i_top, j) + (1.0 - th) * z_raw(i_top - 1, j);
      }

      if (i_bottom < N_global - 1)
      {
        for (int j = 0; j < N_global; ++j)
          z_sm(i_bottom, j) = th * z_raw(i_bottom, j) + (1.0 - th) * z_raw(i_bottom + 1, j);
      }
    }

    // Vertical seams (between patch columns)
    for (int pc = 0; pc < patches_per_side; ++pc)
    {
      const int j_left = n_local * pc;
      const int j_right = n_local * (pc + 1) - 1;

      if (j_left > 0)
      {
        for (int i = 0; i < N_global; ++i)
          z_sm(i, j_left) = th * z_raw(i, j_left) + (1.0 - th) * z_sm(i, j_left - 1);
      }

      if (j_right < N_global - 1)
      {
        for (int i = 0; i < N_global; ++i)
          z_sm(i, j_right) = th * z_raw(i, j_right) + (1.0 - th) * z_sm(i, j_right + 1);
      }
    }

    // Shift global minimum to 0.
    const double zmin = MinOfMatrix(z_sm, N_global);
    std::cout << "zmin=" << zmin << "\n";
    for (int i = 0; i < N_global; ++i)
      for (int j = 0; j < N_global; ++j)
      {
#if (REGULARMIRCO_ELSEDODELTARELATIVETOOTHER)
        z_sm(i, j) -= zmin;
#else
        z_sm(i, j);
#endif
      }

    return z_sm;
  }


}  // namespace MIRCO
