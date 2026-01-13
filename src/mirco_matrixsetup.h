#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Create the influence coefficient matrix (Discrete version of Green's function)
   *
   * @param[in] xv0 x-coordinates of the points in contact in the previous iteration.
   * @param[in] yv0 y-coordinates of the points in contact in the previous iteration.
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] CompositeYoungs The composite Young's modulus
   * @param[in] systemsize Number of nodes predicted to be in contact
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead
   * of point force
   *
   * @return Influence coefficient matrix (Discrete version of Green Function) (usually denoted H)
   */
  ViewMatrix_d SetupMatrix(const ViewVector_d xv0, const ViewVector_d yv0, const double GridSize,
      const double CompositeYoungs, const int systemsize, const bool PressureGreenFunFlag);

  /**
   * @brief Compute one entry of the full influence coefficient matriix. Use when memory is too
   * constrained to store the full matrix.
   *
   * @param[in] ix x index of first point
   * @param[in] iy y index of first point
   * @param[in] jx x index of second point
   * @param[in] jy y index of second point
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] CompositeYoungs The composite Young's modulus
   * @param[in] N Element count along one direction
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead
   * of point force
   *
   * @return Matrix entry of H
   */
  double SetupMatrixOneEntry(const int ix, const int iy, const int jx, const int jy,
      const double GridSize, const double CompositeYoungs, const int N,
      const bool PressureGreenFunFlag);
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_H_
