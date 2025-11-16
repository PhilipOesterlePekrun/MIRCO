#ifndef SRC_VISUALIZATIONEXPORT_H_
#define SRC_VISUALIZATIONEXPORT_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  struct VisualizationData
  {
    ViewMatrix_d topology;
    ViewVector_d meshgrid;

    ViewVector_d xv;
    ViewVector_d yv;
    ViewVector_d activeSet;
    ViewVector_d p;

    std::optional<double> delta;

    void computeFullDisplacementField();
  };
  struct SolverData
  {
    // iteration counts and such. Maybe timers here? #ifdef timers_on?
    // also if you choose, you can store the forces or whatever at each iteration of whatever
  };

  void ExportVisualizations(VisualizationData& visData);
}  // namespace MIRCO

#endif  // SRC_VISUALIZATIONEXPORT_H_
