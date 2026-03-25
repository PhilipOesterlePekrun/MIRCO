#include <filesystem>
#include <fstream>
#include <sstream>

#include "mirco_inputparameters.h"
#include "mirco_inpututilities.h"

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
  ryml::ConstNodeRef geometricalParameters = parameters["geometrical_parameters"];
  ryml::ConstNodeRef materialParameters = parameters["material_parameters"];

  if (root.invalid()) throw std::runtime_error("Input incomplete: missing root `mirco_input`");
  if (parameters.invalid())
    throw std::runtime_error("Input incomplete: missing section `parameters`");
  if (geometricalParameters.invalid())
    throw std::runtime_error("Input incomplete: missing section `geometrical_parameters`");
  if (materialParameters.invalid())
    throw std::runtime_error("Input incomplete: missing section `material_parameters`");

  auto exportVisualization = rget_optional<bool>(root, "ExportVisualization");
  std::optional<std::string> exportVisualizationPath;
  if (exportVisualization && exportVisualization.value())
    exportVisualizationPath = rget<std::string>(root, "ExportVisualizationPath");
  else
    exportVisualizationPath = std::nullopt;

  // Set the surface generator based on RandomTopologyFlag
  if (rget<bool>(root, "RandomTopologyFlag"))
  {
    *this = InputParameters(rget<double>(materialParameters, "E1"),
        rget<double>(materialParameters, "E2"), rget<double>(materialParameters, "nu1"),
        rget<double>(materialParameters, "nu2"), rget<double>(geometricalParameters, "Tolerance"),
        rget<double>(geometricalParameters, "LateralLength"),
        rget<int>(geometricalParameters, "Resolution"),
        rget<double>(geometricalParameters, "InitialTopologyStdDeviation"),
        rget<double>(geometricalParameters, "HurstExponent"), rget<int>(root, "MaxIteration"),
        rget<bool>(root, "WarmStartingFlag"), rget<bool>(root, "PressureGreenFunFlag"),
        rget<bool>(root, "RandomSeedFlag"), rget_optional<int>(root, "RandomGeneratorSeed"),
        exportVisualizationPath);
  }
  else
  {
    std::string topology_file_path = rget<std::string>(root, "TopologyFilePath");
    // If the path is relative, it is relative to the input (.yaml) file
    std::filesystem::path new_path = topology_file_path;
    if (new_path.is_relative())
      new_path = std::filesystem::path(inputFileName).parent_path() / new_path;
    topology_file_path = new_path.string();

    *this = InputParameters(rget<double>(materialParameters, "E1"),
        rget<double>(materialParameters, "E2"), rget<double>(materialParameters, "nu1"),
        rget<double>(materialParameters, "nu2"), rget<double>(geometricalParameters, "Tolerance"),
        rget<double>(geometricalParameters, "LateralLength"), topology_file_path,
        rget<int>(root, "MaxIteration"), rget<bool>(root, "WarmStartingFlag"),
        rget<bool>(root, "PressureGreenFunFlag"), exportVisualizationPath);
  }
}
