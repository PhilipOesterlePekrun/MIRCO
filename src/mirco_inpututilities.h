#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <optional>
#include <ryml.hpp>
#include <ryml_std.hpp>
#include <string>

namespace MIRCO
{
  /*
   * @brief The following function is used to get the value of a parameter in a list (node) in
   * a yaml tree
   */
  template <class T>
  T rget(ryml::ConstNodeRef node, c4::csubstr key)
  {
    auto child = node.find_child(key);
    if (!child.readable())
      throw std::runtime_error("Parameter \"" + std::string(key.str, key.len) + "\" not found");

    T value{};
    if (!ryml::read(child, &value))
      throw std::runtime_error(
          "Parameter \"" + std::string(key.str, key.len) + "\" has invalid value");

    return value;
  }

  /*
   * \brief The following function is for getting an optional parameter
   */
  template <class T>
  std::optional<T> rget_optional(ryml::ConstNodeRef node, c4::csubstr key)
  {
    auto child = node.find_child(key);
    if (!child.readable()) return std::nullopt;

    T value{};
    if (!ryml::read(child, &value))
    {
      throw std::runtime_error(
          "Parameter \"" + std::string(key.str, key.len) + "\" has invalid value");
    }

    return value;
  }
}  // namespace MIRCO

#endif  // SRC_UTILS_H_
