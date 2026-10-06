#ifndef SRC_STACKEDTIMER_H_
#define SRC_STACKEDTIMER_H_

#include <Teuchos_StackedTimer.hpp>
#include <algorithm>
#include <iomanip>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>

namespace MIRCO
{
  // Teuchos' tabular report rounds seconds in internal streams. Read the stored
  // timings directly to retain seven significant figures in a serial report.
  class StackedTimer : public Teuchos::StackedTimer
  {
   public:
    using Teuchos::StackedTimer::StackedTimer;

    void report(std::ostream& out)
    {
      flatten();
      int nameWidth = 0;
      for (const auto& fullName : flat_names_)
      {
        const int level = std::count(fullName.begin(), fullName.end(), '@');
        const auto separator = fullName.rfind('@');
        const auto label = fullName.substr(separator == std::string::npos ? 0 : separator + 1);
        nameWidth = std::max(nameWidth, 4 * level + static_cast<int>(label.size()) + 2);
        nameWidth = std::max(nameWidth, 4 * (level + 1) + 11);  // Remainder row
      }

      std::ostringstream report;
      report << std::left;
      printLevel(report, name(), 0, 0.0, nameWidth);
      out << report.str();
    }

   private:
    double printLevel(std::ostream& out, const std::string& fullName, int level, double parentTime,
        int nameWidth) const
    {
      const auto* timer = findBaseTimer(fullName);
      const double seconds = timer->accumulatedTime();
      const auto separator = fullName.rfind('@');
      const auto label = fullName.substr(separator == std::string::npos ? 0 : separator + 1);
      printRow(out, label, level, seconds, parentTime, nameWidth, timer->numCalls());

      double childTime = 0.0;
      const auto prefix = fullName + '@';
      for (const auto& child : flat_names_)
        if (child.compare(0, prefix.size(), prefix) == 0 &&
            child.find('@', prefix.size()) == std::string::npos)
          childTime += printLevel(out, child, level + 1, seconds, nameWidth);

      if (childTime > 0.0)
        printRow(out, "Remainder", level + 1, seconds - childTime, seconds, nameWidth);
      return seconds;
    }

    static void printRow(std::ostream& out, const std::string& name, int level, double seconds,
        double parentTime, int nameWidth, std::optional<unsigned long> calls = std::nullopt)
    {
      std::string label;
      for (int i = 0; i < level; ++i) label += "|   ";
      label += name + ": ";

      std::ostringstream percentage;
      if (parentTime > 0.0) percentage << " - " << seconds / parentTime * 100 << '%';

      // Scientific notation has one leading digit and six digits after the point.
      out << std::setw(nameWidth) << label << std::scientific << std::setprecision(6)
          << std::setw(14) << seconds << std::setw(14) << percentage.str();
      if (calls) out << " [" << *calls << ']';
      out << '\n';
    }
  };
}  // namespace MIRCO

#endif  // SRC_STACKEDTIMER_H_
