#ifndef SRC_TOPOLOGY_H_
#define SRC_TOPOLOGY_H_

#include <optional>
#include <string>

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Construct a topology by reading topology from an input (.dat) file.
   *
   * @param[in] filepath Path of the input file containing the topology, relative to the calling
   * directory or absolute
   *
   * @return Topology heightfield matrix
   */
  ViewMatrix_h CreateSurfaceFromFile(const std::string& filepath);

  /**
   * @brief Construct a topology using Random Midpoint Generator.
   *
   * @param[in] resolution Resolution parameter
   * @param[in] initialTopologyStdDeviation Initial Standard deviation for the random-midpoint
   * generator [micrometers]
   * @param[in] hurst Hurst exponent
   * @param[in] randomGeneratorSeed Seed for the random mid-point generator
   *
   * @return Topology heightfield matrix
   */
  ViewMatrix_h CreateRmgSurface(int Resolution, double InitialTopologyStdDeviation, double Hurst,
      bool RandomSeedFlag, std::optional<int> RandomGeneratorSeed);

}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
