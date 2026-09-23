cmake_minimum_required(VERSION 3.16)

project(otswap VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(OTSWAP_BUILD_TESTS "Build the C++ test suite" ON)
option(OTSWAP_BUILD_EXAMPLES "Build the examples" ON)

# ---------------------------------------------------------------------------
# meshsearch
#
# It ships no CMake package configuration, so find_package(meshsearch CONFIG)
# is not available: its install step writes nothing at all. It is therefore
# consumed from source, against a pinned tag.
# ---------------------------------------------------------------------------

set(OTSWAP_MESHSEARCH_TAG "v0.1.0" CACHE STRING
    "meshsearch tag to build against")
set(OTSWAP_MESHSEARCH_REPOSITORY "https://github.com/SimoSartori/meshsearch.git"
    CACHE STRING "meshsearch repository; may be a local path")

include(FetchContent)

set(MESHSEARCH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(MESHSEARCH_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(MESHSEARCH_BUILD_PYTHON OFF CACHE BOOL "" FORCE)

FetchContent_Declare(meshsearch
  GIT_REPOSITORY ${OTSWAP_MESHSEARCH_REPOSITORY}
  GIT_TAG        ${OTSWAP_MESHSEARCH_TAG}
  GIT_SHALLOW    TRUE)

FetchContent_MakeAvailable(meshsearch)

# ---------------------------------------------------------------------------
# External dependencies
# ---------------------------------------------------------------------------

# AppleClang needs to be pointed at a libomp; Homebrew's is the usual one.
if (APPLE AND NOT DEFINED OpenMP_ROOT)
  foreach (candidate /opt/homebrew/opt/libomp /usr/local/opt/libomp)
    if (EXISTS ${candidate}/include/omp.h)
      set(OpenMP_ROOT ${candidate})
      break()
    endif()
  endforeach()
endif()

find_package(OpenMP COMPONENTS CXX)

find_path(CFITSIO_INCLUDE_DIR fitsio.h
  HINTS ${CFITSIO_ROOT} ENV CFITSIO_ROOT
  PATH_SUFFIXES include include/cfitsio cfitsio)
find_library(CFITSIO_LIBRARY NAMES cfitsio
  HINTS ${CFITSIO_ROOT} ENV CFITSIO_ROOT PATH_SUFFIXES lib)

if (NOT CFITSIO_INCLUDE_DIR OR NOT CFITSIO_LIBRARY)
  message(FATAL_ERROR
    "cfitsio was not found. Set CFITSIO_ROOT, or pass -DCFITSIO_INCLUDE_DIR "
    "and -DCFITSIO_LIBRARY.")
endif()

find_path(HEALPIX_INCLUDE_DIR healpix_base.h
  HINTS ${HEALPIX_ROOT} ENV HEALPIX_ROOT
  PATH_SUFFIXES include/healpix_cxx include healpix_cxx)
find_library(HEALPIX_LIBRARY NAMES healpix_cxx
  HINTS ${HEALPIX_ROOT} ENV HEALPIX_ROOT PATH_SUFFIXES lib)

if (NOT HEALPIX_INCLUDE_DIR OR NOT HEALPIX_LIBRARY)
  message(FATAL_ERROR
    "healpix_cxx was not found. Set HEALPIX_ROOT to a Healpix C++ installation "
    "(the directory holding include/healpix_cxx and lib), or pass "
    "-DHEALPIX_INCLUDE_DIR and -DHEALPIX_LIBRARY.")
endif()

# ---------------------------------------------------------------------------
# The library
# ---------------------------------------------------------------------------

add_library(otswap STATIC
  src/DistanceTable.cpp
  src/Lightcone.cpp
  src/Mask.cpp
  src/Reconstruct.cpp
  src/Support.cpp
  src/Table.cpp)

set_target_properties(otswap PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(otswap
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CFITSIO_INCLUDE_DIR}
    ${HEALPIX_INCLUDE_DIR})

target_compile_features(otswap PUBLIC cxx_std_17)

target_link_libraries(otswap
  PUBLIC  meshsearch
  PRIVATE ${CFITSIO_LIBRARY} ${HEALPIX_LIBRARY})

if (OpenMP_CXX_FOUND)
  target_link_libraries(otswap PUBLIC OpenMP::OpenMP_CXX)
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

  foreach (test_name test_support test_distance_table test_reconstruct test_table test_mask)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE otswap)
    target_include_directories(${test_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/src ${CFITSIO_INCLUDE_DIR} ${HEALPIX_INCLUDE_DIR})
    target_link_libraries(${test_name} PRIVATE ${CFITSIO_LIBRARY} ${HEALPIX_LIBRARY})
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(${test_name} PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()
endif()

if (OTSWAP_BUILD_EXAMPLES)
  add_executable(example_box examples/example_box.cpp)
  target_link_libraries(example_box PRIVATE otswap)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(example_box PRIVATE -Wall -Wextra)
  endif()
endif()
