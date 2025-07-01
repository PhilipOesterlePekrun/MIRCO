#include "mirco_topologyutilities.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "mirco_topology.h"

void MIRCO::CreateMeshgrid(std::vector<double>& meshgrid, const int ngrid, const double GridSize)
{
#pragma omp parallel for schedule(static, 16)  // Same amount of work -> static
  for (int i = 0; i < ngrid; i++)
  {
    meshgrid[i] = (GridSize / 2) + i * GridSize;
  }
}

void MIRCO::CreateSurfaceObject(const int Resolution, const double InitialTopologyStdDeviation,
    const double Hurst, const bool RandomSeedFlag, const std::string TopologyFilePath,
    const bool RandomTopologyFlag, const int RandomGeneratorSeed,
    Teuchos::RCP<MIRCO::TopologyGeneration>& surfacegenerator)  // # dont need this maybe
{
  if (RandomTopologyFlag)
  {
    surfacegenerator = Teuchos::rcp(new MIRCO::Rmg(
        Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed));
  }
  else
  {
    surfacegenerator = Teuchos::rcp(new MIRCO::ReadFile(Resolution, TopologyFilePath));
  }
}

MIRCO::TopologyMaxAndMean MIRCO::ComputeMaxAndMean(
    const Teuchos::SerialDenseMatrix<int, double>& topology)
{
  double zmax = -1.0;
  double zmean = 0.0;
#pragma omp parallel for schedule(guided, 16) reduction(+ : zmean) reduction(max : zmax)
  // Static and Guided seem even but Guided makes more sense
  for (int i = 0; i < topology.numRows(); i++)
  {
    for (int j = 0; j < topology.numCols(); j++)
    {
      zmean += topology(i, j);
      if (topology(i, j) > zmax)
      {
        zmax = topology(i, j);
      }
    }
  }

  zmean = zmean / (topology.numCols() * topology.numRows());
  return TopologyMaxAndMean{zmax, zmean};
}
