#ifndef SRC_MATRIXSETUP_KOKKOS_H_
#define SRC_MATRIXSETUP_KOKKOS_H_

namespace MIRCO
{
  namespace MatrixGeneration
  {
    #include "mirco_kokkostypes_kokkos.h"
    /**
     * @brief The aim of this function is to create the influence coefficient matrix (Discrete
     * version of Green function)
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
    ViewMatrix_d SetupMatrix(const ViewVector_d xv0, const ViewVector_d yv0,
    const double GridSize, const double CompositeYoungs, const int systemsize,
    const bool PressureGreenFunFlag);
  };  // namespace MatrixGeneration
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_KOKKOS_H_
