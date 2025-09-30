#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <type_traits>

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  namespace Utils
  {
    /**
     * @brief Necessary helper function for kokkos_createMirrorViewAndDeepCopyAnyLayout()
     */
    template <class ViewDst, class ViewSrc, std::size_t... Extents>
    ViewDst kokkos_allocateFromExtents(
        const std::string& label, const ViewSrc& vSrc, std::index_sequence<Extents...>)
    {
      if constexpr (sizeof...(Extents) == 0)
      {
        return ViewDst(Kokkos::view_alloc(Kokkos::WithoutInitializing, label));
      }
      else
      {
        return ViewDst(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, label), vSrc.extent(Extents)...);
      }
    }

    /**
     * @brief Helper to create a mirror view and copy across spaces, even if Layout differs
     */
    template <class ExecDst, class ViewSrc, class ViewDst>
    ViewDst kokkos_createMirrorViewAndDeepCopyAnyLayout(
        const ViewSrc vSrc, const std::string& label)
    {
      auto vInter = Kokkos::create_mirror_view_and_copy(ExecDst{}, vSrc);

      if constexpr (std::is_same_v<typename ViewSrc::array_layout, typename ViewDst::array_layout>)
      {
        // Layouts match -> shallow copy
        return vInter;
      }
      else
      {
        // Transform to Dst Layout
        auto vDst = kokkos_allocateFromExtents<ViewDst>(
            label, vSrc, std::make_index_sequence<ViewSrc::rank>{});
        Kokkos::deep_copy(ExecDst{}, vDst, vInter);
        return vDst;
      }
    }
  }  // namespace Utils
}  // namespace MIRCO

#endif  // SRC_UTILS_H_
