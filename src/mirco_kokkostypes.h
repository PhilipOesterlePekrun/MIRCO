#ifndef SRC_KOKKOSTYPES_H_
#define SRC_KOKKOSTYPES_H_

#include <Kokkos_Core.hpp>
#include <type_traits>
#include <utility>  //#

// This file defines some commonly used Kokkos type aliases. More aliases, as well as Kokkos-related
// macros, utilities, etc., can be added as needed
namespace MIRCO
{
  using ExecSpace_DefaultHost_t = Kokkos::DefaultHostExecutionSpace;
  using ExecSpace_Default_t = Kokkos::DefaultExecutionSpace;

  using MemorySpace_Host_t = Kokkos::HostSpace;
  using MemorySpace_ofDefaultExec_t = ExecSpace_Default_t::memory_space;

  using Device_Host_t = Kokkos::Device<ExecSpace_DefaultHost_t, MemorySpace_Host_t>;
  using Device_Default_t = Kokkos::Device<ExecSpace_Default_t, MemorySpace_ofDefaultExec_t>;

  using ViewVector_h = Kokkos::View<double*, Device_Host_t>;
  using ViewMatrix_h = Kokkos::View<double**, Device_Host_t>;

  using ViewVector_d = Kokkos::View<double*, Device_Default_t>;
  using ViewMatrix_d = Kokkos::View<double**, Device_Default_t>;

  using ViewScalarInt_d = Kokkos::View<int, Device_Default_t>;
  using ViewVectorInt_d = Kokkos::View<int*, Device_Default_t>;


  /*template <class ViewDst, class ViewSrc, std::size_t... I>
  inline ViewDst allocate_like_impl(
      std::string const& label, ViewSrc const& vSrc, std::index_sequence<I...>)
  {
    return ViewDst{Kokkos::view_alloc(label, Kokkos::WithoutInitializing{}), vSrc.extent(I)...};
  }

  /**
   * @brief Helper to create a mirror view and copy across spaces, even if Layout differs
   */
  /*template <class ExecDst, class ViewSrc, class ViewDst>
  ViewDst Kokkos_CreateMirrorViewAndDeepCopyAnyLayout(const ViewSrc vSrc, const std::string& label)
  {
    using SrcLayout = typename ViewSrc::array_layout;
    using DstLayout = typename ViewDst::array_layout;

    auto vInter = Kokkos::create_mirror_view_and_copy(ExecDst{}, vSrc);

    if constexpr (std::is_same_v<SrcLayout, DstLayout>)
    {
      // Layouts match -> shallow copy
      return vInter;
    }
    else
    {
      // Transform to Dst layout
      ViewDst vDst =
          allocate_like_impl<ViewDst>(label, vSrc, std::make_index_sequence<ViewSrc::rank>{});
      Kokkos::deep_copy(ExecDst{}, vDst, vInter);
      return vDst;
    }
  }*/

}  // namespace MIRCO

#endif  // SRC_KOKKOSTYPES_H_
