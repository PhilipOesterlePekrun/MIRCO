#include <fstream>
#include <sstream>

#include "mirco_inputparameters.h"
#include "mirco_utilsIO.h"

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  std::ifstream fin(inputFileName);
  if (!fin) throw std::runtime_error("Cannot open file: " + inputFileName);

  std::stringstream ss;
  ss << fin.rdbuf();
  std::string inString = ss.str();

  ryml::Tree tree = ryml::parse_in_arena(c4::to_csubstr(inString));
  ryml::ConstNodeRef root = tree["mirco_input"];
  ryml::ConstNodeRef parameters = root["parameters"];
  ryml::ConstNodeRef geoParams = parameters["geometrical_parameters"];
  ryml::ConstNodeRef matParams = parameters["material_parameters"];

  if (root.invalid()) throw std::runtime_error("Input incomplete: missing root `mirco_input`");
  if (parameters.invalid())
    throw std::runtime_error("Input incomplete: missing section `parameters`");
  if (geoParams.invalid())
    throw std::runtime_error("Input incomplete: missing section `geometrical_parameters`");
  if (matParams.invalid())
    throw std::runtime_error("Input incomplete: missing section `material_parameters`");

  // Set the surface generator based on RandomTopologyFlag
  if (UtilsIO::get_bool(root, "RandomTopologyFlag"))
  {
    *this = InputParameters(UtilsIO::get_double(matParams, "E1"),
        UtilsIO::get_double(matParams, "E2"), UtilsIO::get_double(matParams, "nu1"),
        UtilsIO::get_double(matParams, "nu2"), UtilsIO::get_double(geoParams, "Tolerance"),
        UtilsIO::get_double(geoParams, "Delta"), UtilsIO::get_double(geoParams, "LateralLength"),
        UtilsIO::get_int(geoParams, "Resolution"),
        UtilsIO::get_double(geoParams, "InitialTopologyStdDeviation"),
        UtilsIO::get_double(geoParams, "HurstExponent"), UtilsIO::get_bool(root, "RandomSeedFlag"),
        UtilsIO::get_int(root, "RandomGeneratorSeed"), UtilsIO::get_int(root, "MaxIteration"),
        UtilsIO::get_bool(root, "WarmStartingFlag"),
        UtilsIO::get_bool(root, "PressureGreenFunFlag"));
  }
  else
  {
    std::string topology_file_path = UtilsIO::get_string(root, "TopologyFilePath");
    // The following function generates the actual path of the topology file
    MIRCO::UtilsIO::changeRelativePath(topology_file_path, inputFileName);

    *this =
        InputParameters(UtilsIO::get_double(matParams, "E1"), UtilsIO::get_double(matParams, "E2"),
            UtilsIO::get_double(matParams, "nu1"), UtilsIO::get_double(matParams, "nu2"),
            UtilsIO::get_double(geoParams, "Tolerance"), UtilsIO::get_double(geoParams, "Delta"),
            UtilsIO::get_double(geoParams, "LateralLength"), topology_file_path,
            UtilsIO::get_int(root, "MaxIteration"), UtilsIO::get_bool(root, "WarmStartingFlag"),
            UtilsIO::get_bool(root, "PressureGreenFunFlag"));
  }
}
