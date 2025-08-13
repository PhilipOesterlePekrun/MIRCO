#ifndef SRC_WARMSTART_KOKKOS_H_
#define SRC_WARMSTART_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief This function is used to determine the nodes which were in contact in the last iteration
   * from the current contact set. This helps in making an initial guess of the nodes in contact in
   * the current iteration and speeds up the computation.
   *
   * @param[in] xv0 x-coordinates of the points in contact in the previous iteration.
   * @param[in] yv0 y-coordinates of the points in contact in the previous iteration.
   * @param[in] xvf x-coordinates of the points in contact in the previous iteration.
   * @param[in] yvf y-coordinates of the points in contact in the previous iteration.
   * @param[in] pf Contact force at (xvf,yvf) predicted in the previous iteration.
   *
   * @return p0 vector of contact forces predicted in the previous iteration but are a part of
   * currect predicted contact set.
   */
  ViewVector_h Warmstart(const ViewVector_h& xv0, const ViewVector_h& yv0, const ViewVector_h& xvf,
      const ViewVector_h& yvf, const ViewVector_h& pf);
}  // namespace MIRCO

#endif  // SRC_WARMSTART_KOKKOS_H_
