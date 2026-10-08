# Includes the otswap source tree as a subproject. Configure with
# -DOTSWAP_SOURCE_DIR=<the otswap checkout>, CFITSIO_ROOT where cfitsio is
# not in a system path, and no CMAKE_BUILD_TYPE: as a subproject otswap must
# build neither its tests nor the tests' Healpix, and must leave the build
# type to the parent, so it has to stay empty here. The meshsearch it fetches
# must not build its own tests or example either.

cmake_minimum_required(VERSION 3.16)

project(otswap_subproject_consumer LANGUAGES C CXX)

add_subdirectory(${OTSWAP_SOURCE_DIR} otswap)

foreach (target test_support test_distance_table test_reconstruct test_table test_mask
                test_detmath test_pixel test_rsd cpp_reference otswap_healpix test_meshgrid
                example)
  if (TARGET ${target})
    message(FATAL_ERROR "otswap built ${target} as a subproject")
  endif()
endforeach()

if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
  message(FATAL_ERROR "otswap set CMAKE_BUILD_TYPE to '${CMAKE_BUILD_TYPE}' as a subproject")
endif()

add_executable(consumer ../consumer.cpp)
target_link_libraries(consumer PRIVATE otswap::otswap)
