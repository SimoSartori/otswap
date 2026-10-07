cmake_minimum_required(VERSION 3.16)

project(otswap VERSION 0.1.0 LANGUAGES C CXX)

# PROJECT_IS_TOP_LEVEL needs CMake 3.21; this is the 3.16 equivalent, true
# for a standalone build and false under add_subdirectory or FetchContent.
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(OTSWAP_IS_TOP_LEVEL ON)
else()
  set(OTSWAP_IS_TOP_LEVEL OFF)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Only a top-level build may choose these for the whole build tree; as a
# subproject, the parent's choices stand.
if (OTSWAP_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

# scikit-build-core defines SKBUILD. The wheel carries the Python module only,
# so the C++ install rules are off there.
if (DEFINED SKBUILD)
  set(OTSWAP_INSTALL_DEFAULT OFF)
else()
  set(OTSWAP_INSTALL_DEFAULT ON)
endif()

option(OTSWAP_BUILD_TESTS "Build the C++ test suite" ${OTSWAP_IS_TOP_LEVEL})
option(OTSWAP_BUILD_EXAMPLES "Build the examples" ${OTSWAP_IS_TOP_LEVEL})
option(OTSWAP_BUILD_PYTHON "Build the Python bindings" OFF)
option(OTSWAP_INSTALL "Generate the C++ install rules" ${OTSWAP_INSTALL_DEFAULT})

include(GNUInstallDirs)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ---------------------------------------------------------------------------
# Floating-point contraction
#
# No a*b+c is fused into one FMA anywhere: compilers differ in what they
# contract by default, and so do the instruction sets, so contraction would
# make the results depend on the platform. meshsearch sets it on its own
# targets; this covers otswap's, fdlibm, a fetched meshsearch, the tests'
# Healpix and the extension module. It is not exported: all of otswap's
# arithmetic is in its own archive.
# ---------------------------------------------------------------------------

if (NOT MSVC)
  add_compile_options(-ffp-contract=off)
endif()

# ---------------------------------------------------------------------------
# meshsearch
#
# An installed meshsearch of the same minor version is used when found;
# otherwise the pinned tag is fetched and built from source. A parent
# project that already defines meshsearch::meshsearch provides its own.
# ---------------------------------------------------------------------------

set(OTSWAP_MESHSEARCH_VERSION 0.2)
set(OTSWAP_MESHSEARCH_TAG "v0.2.0" CACHE STRING
    "meshsearch tag to fetch when no installed meshsearch is found")
set(OTSWAP_MESHSEARCH_REPOSITORY "https://github.com/SimoSartori/meshsearch.git"
    CACHE STRING "meshsearch repository to fetch from; may be a local path")

if (NOT TARGET meshsearch::meshsearch)
  find_package(meshsearch ${OTSWAP_MESHSEARCH_VERSION} CONFIG QUIET)
endif()

if (TARGET meshsearch::meshsearch)
  if (meshsearch_FOUND)
    message(STATUS "otswap: using meshsearch ${meshsearch_VERSION} from ${meshsearch_DIR}")
  else()
    message(STATUS "otswap: using the parent project's meshsearch::meshsearch")
  endif()
else()
  message(STATUS "otswap: no installed meshsearch ${OTSWAP_MESHSEARCH_VERSION}; "
                 "fetching ${OTSWAP_MESHSEARCH_TAG}")

  # meshsearch's install rules follow otswap's, so that an installed otswap
  # comes with the meshsearch it was built against, and the wheel installs
  # neither. A normal variable overrides meshsearch's option (CMP0077).
  set(MESHSEARCH_INSTALL ${OTSWAP_INSTALL})

  include(FetchContent)
  FetchContent_Declare(meshsearch
    GIT_REPOSITORY ${OTSWAP_MESHSEARCH_REPOSITORY}
    GIT_TAG        ${OTSWAP_MESHSEARCH_TAG}
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(meshsearch)
endif()

# ---------------------------------------------------------------------------
# External dependencies
# ---------------------------------------------------------------------------

# AppleClang needs to be pointed at a libomp; Homebrew's is the usual one,
# used unless OpenMP_ROOT is set, as a CMake variable or in the environment.
# otswapConfig.cmake repeats this for consumers.
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND
    NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
  foreach (candidate /opt/homebrew/opt/libomp /usr/local/opt/libomp)
    if (EXISTS ${candidate}/include/omp.h)
      set(OpenMP_ROOT ${candidate})
      break()
    endif()
  endforeach()
endif()

find_package(OpenMP COMPONENTS CXX)

# cmake/FindCFITSIO.cmake, installed with the package for consumers.
find_package(CFITSIO)
if (NOT CFITSIO_FOUND)
  message(FATAL_ERROR
    "cfitsio was not found. Set CFITSIO_ROOT, or pass -DCFITSIO_INCLUDE_DIR "
    "and -DCFITSIO_LIBRARY.")
endif()

# ---------------------------------------------------------------------------
# Deterministic math, vendored
#
# sin, cos, exp, log1p, pow, atan2 and asin from FreeBSD's msun, copied
# unmodified into extern/fdlibm/msun/src (see its README), so that every
# platform computes them with the same code rather than with its own libm.
# Third-party code that must not be edited: its warnings are suppressed, and
# extern/fdlibm/otswap_fdlibm.h is force-included to rename its symbols.
# Only src/detmath.h calls it. An object library, whose objects go into
# libotswap.a: nothing of it is installed or exported on its own.
# ---------------------------------------------------------------------------

set(OTSWAP_FDLIBM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/fdlibm)

add_library(otswap_fdlibm OBJECT
  ${OTSWAP_FDLIBM_DIR}/msun/src/s_sin.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/s_cos.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/k_sin.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/k_cos.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/k_rem_pio2.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/e_exp.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/s_log1p.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/e_pow.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/e_atan2.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/s_atan.c
  ${OTSWAP_FDLIBM_DIR}/msun/src/e_asin.c)

# GNU C: math_private.h needs the BSD integer types and endianness macros,
# which glibc declares only outside strict ISO mode.
set_target_properties(otswap_fdlibm PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  C_STANDARD 99
  C_STANDARD_REQUIRED ON
  C_EXTENSIONS ON)

if (NOT APPLE)
  target_include_directories(otswap_fdlibm PRIVATE ${OTSWAP_FDLIBM_DIR}/shim)
endif()

if (MSVC)
  message(FATAL_ERROR "the vendored msun sources are built with GCC or clang only")
else()
  target_compile_options(otswap_fdlibm PRIVATE -w -include ${OTSWAP_FDLIBM_DIR}/otswap_fdlibm.h)
endif()

# ---------------------------------------------------------------------------
# The library
# ---------------------------------------------------------------------------

add_library(otswap STATIC
  src/DistanceTable.cpp
  src/Lightcone.cpp
  src/Mask.cpp
  src/Pixel.cpp
  src/Reconstruct.cpp
  src/RSD.cpp
  src/Support.cpp
  src/Table.cpp
  $<TARGET_OBJECTS:otswap_fdlibm>)
add_library(otswap::otswap ALIAS otswap)

set_target_properties(otswap PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(otswap
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_features(otswap PUBLIC cxx_std_17)

# meshsearch is public: the documented meshsearch::Error can reach callers.
# OpenMP is private: neither public header needs it, and a consumer of the
# static archive gets its link requirement only.
target_link_libraries(otswap
  PUBLIC  meshsearch::meshsearch
  PRIVATE CFITSIO::CFITSIO)

set(OTSWAP_WITH_OPENMP OFF)
if (OpenMP_CXX_FOUND)
  target_link_libraries(otswap PRIVATE OpenMP::OpenMP_CXX)
  set(OTSWAP_WITH_OPENMP ON)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(otswap PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Tests and examples
# ---------------------------------------------------------------------------

if (OTSWAP_BUILD_TESTS)
  include(CTest)
  enable_testing()

  # Healpix, vendored: HEALPix 3.82, copied unmodified into extern/healpix
  # (see its README). The library has its own pixel code; the tests compare
  # it with Healpix's. Third-party code that must not be edited, so its
  # warnings are suppressed and its include directories are SYSTEM.
  set(OTSWAP_HEALPIX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/healpix)

  add_library(otswap_healpix STATIC
    ${OTSWAP_HEALPIX_DIR}/Healpix_cxx/healpix_base.cc
    ${OTSWAP_HEALPIX_DIR}/Healpix_cxx/healpix_tables.cc
    ${OTSWAP_HEALPIX_DIR}/cxxsupport/error_handling.cc
    ${OTSWAP_HEALPIX_DIR}/cxxsupport/geom_utils.cc
    ${OTSWAP_HEALPIX_DIR}/cxxsupport/pointing.cc
    ${OTSWAP_HEALPIX_DIR}/cxxsupport/string_utils.cc)

  target_compile_features(otswap_healpix PUBLIC cxx_std_17)
  target_include_directories(otswap_healpix SYSTEM PUBLIC
    ${OTSWAP_HEALPIX_DIR}/Healpix_cxx
    ${OTSWAP_HEALPIX_DIR}/cxxsupport)

  if (MSVC)
    target_compile_options(otswap_healpix PRIVATE /w)
  else()
    target_compile_options(otswap_healpix PRIVATE -w)
  endif()

  foreach (test_name test_support test_distance_table test_reconstruct test_table test_mask
                   test_detmath test_pixel test_rsd)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${test_name} PRIVATE otswap CFITSIO::CFITSIO otswap_healpix)
    if (OpenMP_CXX_FOUND)
      target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(${test_name} PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()

  # Runs the library on raw input files for the Python tests, which compare
  # its outputs with those of the bindings; not a test by itself.
  add_executable(cpp_reference tests/cpp_reference.cpp)
  target_link_libraries(cpp_reference PRIVATE otswap)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(cpp_reference PRIVATE -Wall -Wextra)
  endif()
endif()

if (OTSWAP_BUILD_EXAMPLES)
  add_executable(example_box examples/example_box.cpp)
  target_link_libraries(example_box PRIVATE otswap)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(example_box PRIVATE -Wall -Wextra)
  endif()
endif()

# ---------------------------------------------------------------------------
# Python bindings
#
# The extension module otswap._otswap, installed into the otswap package of
# python/otswap. Built by pip through scikit-build-core (see pyproject.toml).
# ---------------------------------------------------------------------------

if (OTSWAP_BUILD_PYTHON)
  find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
  find_package(nanobind CONFIG REQUIRED)

  nanobind_add_module(otswap_python NB_STATIC python/otswap_python.cpp)
  set_target_properties(otswap_python PROPERTIES OUTPUT_NAME _otswap)
  target_link_libraries(otswap_python PRIVATE otswap)
  install(TARGETS otswap_python LIBRARY DESTINATION otswap)
endif()

# ---------------------------------------------------------------------------
# Install and export
#
# find_package(otswap 0.1 CONFIG) then otswap::otswap. The archive holds
# fdlibm's objects; its other dependencies are found again on the consumer's
# machine by otswapConfig.cmake: meshsearch of the same minor version,
# cfitsio through the FindCFITSIO.cmake installed beside it, and OpenMP when
# otswap was built with it.
# ---------------------------------------------------------------------------

if (OTSWAP_INSTALL)
  include(CMakePackageConfigHelpers)

  set(OTSWAP_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/otswap)

  install(TARGETS otswap EXPORT otswapTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/otswap DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(EXPORT otswapTargets
          NAMESPACE otswap::
          DESTINATION ${OTSWAP_CMAKE_DIR})

  configure_package_config_file(cmake/otswapConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/otswapConfig.cmake
    INSTALL_DESTINATION ${OTSWAP_CMAKE_DIR})

  # Before 1.0 a minor release may change the interface, so a request for 0.1
  # accepts any 0.1.x and nothing else.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/otswapConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/otswapConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/otswapConfigVersion.cmake
                cmake/FindCFITSIO.cmake
          DESTINATION ${OTSWAP_CMAKE_DIR})
endif()
