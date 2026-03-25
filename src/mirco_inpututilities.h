#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <optional>
#include <ryml.hpp>
#include <ryml_std.hpp>
#include <string>

namespace MIRCO
{
  /*
   * \brief The following function is used to get the value of a parameter in a list (node) in
   * a yaml tree, using ryml
   */
  template <class T>
  T rget(ryml::ConstNodeRef node, c4::csubstr key)
  {
    auto child = node[key];
    if (child.invalid())
      throw std::runtime_error("Parameter \"" + std::string(key.str, key.len) + "\" not found");

    try
    {
      T value{};
      child >> value;
      return value;
    }
    catch (...)
    {
      throw std::runtime_error(
          "Parameter \"" + std::string(key.str, key.len) + "\" has invalid value");
    }
  }

  /*
   * \brief The following function is for getting an optional parameter
   */
  template <class T>
  std::optional<T> rget_optional(ryml::ConstNodeRef node, c4::csubstr key)
  {
    auto child = node[key];
    if (child.invalid()) return std::nullopt;

    try
    {
      T value{};
      child >> value;
      return value;
    }
    catch (...)
    {
      throw std::runtime_error(
          "Parameter \"" + std::string(key.str, key.len) + "\" has invalid value");
    }
  }
}  // namespace MIRCO

#endif  // SRC_UTILS_H_
