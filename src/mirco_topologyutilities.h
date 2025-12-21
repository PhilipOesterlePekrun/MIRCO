#ifndef SRC_TOPOLOGYUTILITIES_H_
#define SRC_TOPOLOGYUTILITIES_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Create a Meshgrid vector of the surface
   *
   * Creates a meshgrid with coordinates. Since the vector is identical in both directions,
   * therefore only one is created.
   *
   * @param meshgrid Meshgrid vector
   * @param ngrid Number of grid points in one direction
   * @param GridSize Grid size (length of each cell)
   */
  ViewVector_d CreateMeshgrid(const int ngrid, const double GridSize);

  /**
   * @brief Compute the maximum value of a ViewMatrix_d v.
   */
  double GetMax(ViewMatrix_d v);
}  // namespace MIRCO

#endif  // SRC_TOPOLOGYUTILITIES_H_
