#include "mirco_utils.h"

#include <filesystem>

namespace
{
  inline bool has_key(ryml::ConstNodeRef node, const std::string& key)
  {
    return node.has_child(ryml::to_csubstr(key));
  }
}  // namespace

namespace MIRCO::Utils
{
  void changeRelativePath(std::string& targetfilename, const std::string& sourcefilename)
  {
    std::filesystem::path targetfilepath = targetfilename;

    if (targetfilepath.is_relative())
    {
      std::filesystem::path sourcefilepath = sourcefilename;
      std::filesystem::path root_dir = sourcefilepath.parent_path();
      root_dir += "/";
      root_dir += targetfilepath;
      targetfilename = root_dir.c_str();
    }
  }

  std::string get_string(ryml::ConstNodeRef node, const std::string& key)
  {
    auto child = node[ryml::to_csubstr(key)];
    if (child.invalid()) throw std::runtime_error("Parameter \"" + key + "\" not found");
    ryml::csubstr v = child.val();
    return std::string(v.str, v.len);
  }
  bool get_bool(ryml::ConstNodeRef node, const std::string& key)
  {
    std::string s = get_string(node, key);
    if (s == "true" || s == "True" || s == "1") return true;
    if (s == "false" || s == "False" || s == "0") return false;

    throw std::runtime_error("Parameter \"" + key + "\" has invalid bool value");
  }
  double get_double(ryml::ConstNodeRef node, const std::string& key)
  {
    std::string s = get_string(node, key);
    return std::stod(s);
  }
  int get_int(ryml::ConstNodeRef node, const std::string& key)
  {
    std::string s = get_string(node, key);
    return std::stoi(s);
  }

  std::optional<std::string> get_optional_string(ryml::ConstNodeRef node, const std::string& key)
  {
    if (has_key(node, key)) return get_string(node, key);
    return std::nullopt;
  }
  std::optional<bool> get_optional_bool(ryml::ConstNodeRef node, const std::string& key)
  {
    if (has_key(node, key)) return get_bool(node, key);
    return std::nullopt;
  }
  std::optional<double> get_optional_double(ryml::ConstNodeRef node, const std::string& key)
  {
    if (has_key(node, key)) return get_double(node, key);
    return std::nullopt;
  }
  std::optional<int> get_optional_int(ryml::ConstNodeRef node, const std::string& key)
  {
    if (has_key(node, key)) return get_int(node, key);
    return std::nullopt;
  }

}  // namespace MIRCO::Utils
