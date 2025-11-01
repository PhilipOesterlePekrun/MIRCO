cmake_minimum_required(VERSION 3.28)

# try to prevent modification of source directory
# note: some files may still be written before CMake can abort and need to be removed manually
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source build not allowed. "
      "Please create a new directory, preferably next to the source directory, and run CMake from there. "
      "You may want to remove CMakeCache.txt and CMakeFiles/ which were created in the source directory."
    )
endif()

project(mirco VERSION 0.1.0)

# Print CMake version to screen
message(STATUS "Using CMake ${CMAKE_VERSION}")

# Enforce the C++ standard we are using and turn off compiler-specific extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# We do not use C++ modules (yet). Turn off scanning to avoid issues with clang-tidy.
# If you want to add module support, this problem needs to be revisited. Our hope
# is that CMake and/or clang-tidy will be updated to handle modules better.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Ensure cmake setup the correct runtime path when installing
# see here for more information: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
include(GNUInstallDirs)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

option(myUtils_ENABLE_TIMERS "Enable timers" ON)
add_compile_definitions(
  myUtils_ENABLE_TIMERS=$<IF:$<BOOL:${myUtils_ENABLE_TIMERS}>,1,0>
)

# Add libraries
if(GTEST_IN_MIRCO)
  add_subdirectory("extern/googletest")
endif()

# Kokkos
option(KOKKOS_IN_MIRCO "Find Kokkos in the MIRCO project. If set to OFF, use the Kokkos installation from an upstream project" ON)
if(KOKKOS_IN_MIRCO)
  set(CMAKE_PREFIX_PATH ${KOKKOS_PATH} ${CMAKE_PREFIX_PATH})
  find_package(Kokkos REQUIRED)
else()
  # Use a Kokkos from an upstream project
  if(TARGET Kokkos::kokkos)
    message(STATUS "Kokkos included from upstream project.")
  else()
    message(FATAL_ERROR "KOKKOS_IN_MIRCO=OFF but no Kokkos target was available in the upstream project.")
  endif()
endif()

# Print details on the Kokkos configuration
if(KOKKOS_IN_MIRCO)
  if(Kokkos_FOUND)
    message("\nFound Kokkos! Here are the details:")
    get_target_property(Kokkos_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
    message("   Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
    get_target_property(Kokkos_COMPILE_DEFINITIONS Kokkos::kokkos INTERFACE_COMPILE_DEFINITIONS)
    message("   Kokkos_COMPILE_DEFINITIONS = ${Kokkos_COMPILE_DEFINITIONS}")
    get_target_property(Kokkos_COMPILE_OPTIONS Kokkos::kokkos INTERFACE_COMPILE_OPTIONS)
    message("   Kokkos_COMPILE_OPTIONS = ${Kokkos_COMPILE_OPTIONS}")
    message("End of Kokkos details\n")
  else()
    message(FATAL_ERROR "Could not find Kokkos!")
  endif()
endif()

# Kokkos-Kernels
option(KOKKOS_KERNELS_IN_MIRCO "Find Kokkos-Kernels in the MIRCO project. If set to OFF, use the Kokkos-Kernels installation from an upstream project" ON)
if(KOKKOS_KERNELS_IN_MIRCO)
  set(CMAKE_PREFIX_PATH ${KOKKOS_KERNELS_PATH} ${CMAKE_PREFIX_PATH})
  find_package(KokkosKernels REQUIRED)
else()
  # Use a KokkosKernels from an upstream project
  if(TARGET KokkosKernels::kokkoskernels)
    message(STATUS "KokkosKernels included from upstream project.")
  else()
    message(FATAL_ERROR "KOKKOS_KERNELS_IN_MIRCO=OFF but no KokkosKernels target was available in the upstream project.")
  endif()
endif()

# Print details on the Kokkos-Kernels configuration
if(KOKKOS_KERNELS_IN_MIRCO)
  if(KokkosKernels_FOUND)
    message("\nFound KokkosKernels! Here are the details:")
    get_target_property(KokkosKernels_INCLUDE_DIRS KokkosKernels::kokkoskernels INTERFACE_INCLUDE_DIRECTORIES)
    message("   KokkosKernels_INCLUDE_DIRS = ${KokkosKernels_INCLUDE_DIRS}")
    get_target_property(KokkosKernels_COMPILE_DEFINITIONS KokkosKernels::kokkoskernels INTERFACE_COMPILE_DEFINITIONS)
    message("   KokkosKernels_COMPILE_DEFINITIONS = ${KokkosKernels_COMPILE_DEFINITIONS}")
    get_target_property(KokkosKernels_COMPILE_OPTIONS KokkosKernels::kokkoskernels INTERFACE_COMPILE_OPTIONS)
    message("   KokkosKernels_COMPILE_OPTIONS = ${KokkosKernels_COMPILE_OPTIONS}")
    message("End of KokkosKernels details\n")
  else()
    message(FATAL_ERROR "Could not find KokkosKernels!")
  endif()
endif()

# --- myUtils ---
set(MYUTILS_DIR "$ENV{HOME}/misc/myUtils_Base/install")
set(CMAKE_PREFIX_PATH ${MYUTILS_DIR} ${CMAKE_PREFIX_PATH})
find_package(myUtils REQUIRED)

# File(s) which need(s) Kokkos-Kernels
add_library(mirco_needKK
  src/mirco_nonlinearsolver.cpp
  )
target_link_libraries(mirco_needKK PRIVATE Kokkos::kokkos KokkosKernels::kokkoskernels myUtils::myUtils)

# Compile mirco library
add_library(mirco_core
  src/mirco_evaluate.cpp
  src/mirco_matrixsetup.cpp
  src/mirco_contactpredictors.cpp
  src/mirco_contactstatus.cpp
  src/mirco_warmstart.cpp
  )
target_link_libraries(mirco_core PRIVATE mirco_needKK Kokkos::kokkos myUtils::myUtils)

add_library(mirco_topology
  src/mirco_topology.cpp
  src/mirco_topologyutilities.cpp
  )
target_link_libraries(mirco_topology PRIVATE Kokkos::kokkos myUtils::myUtils)

add_library(mirco_inputparameters
  src/mirco_inputparameters.cpp
  )
target_link_libraries(mirco_inputparameters PRIVATE mirco_topology Kokkos::kokkos myUtils::myUtils)

add_library(mirco_lib INTERFACE)
target_link_libraries(mirco_lib INTERFACE mirco_needKK mirco_core mirco_topology mirco_inputparameters myUtils::myUtils)


# arg exe with timers
add_executable(timedMirco src/main_argsAndTimers.cpp)

target_include_directories(mirco_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
  $<INSTALL_INTERFACE:include/mirco/>
  )

add_library(mirco::mirco_lib ALIAS mirco_lib)


target_link_libraries(timedMirco PUBLIC mirco::mirco_lib Kokkos::kokkos myUtils::myUtils)

# Install mirco (to be used as a library by other codes)
install(TARGETS mirco_lib mirco_needKK mirco_core mirco_topology mirco_inputparameters
  EXPORT mirco_libTargets
  ARCHIVE LIBRARY PUBLIC_HEADER
  )

install(EXPORT mirco_libTargets
  NAMESPACE mirco::
  DESTINATION lib/cmake/mirco
  )

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/mirco_libConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/mirco_libConfig.cmake
  INSTALL_DESTINATION lib/cmake/mirco
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/mirco_libConfig.cmake
  DESTINATION lib/cmake/mirco
  )

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include/mirco FILES_MATCHING PATTERN "*.h")
