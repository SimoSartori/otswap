# Builds the example against an installed otswap, from this folder alone.
# Configure with -DCMAKE_PREFIX_PATH=<the prefix otswap was installed to>,
# and CFITSIO_ROOT (and ZLIB_ROOT for a static cfitsio) where cfitsio is not
# in a system path:
#
#   cmake -S examples/periodic_box -B build/examples/periodic_box -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build build/examples/periodic_box
#   build/examples/periodic_box/periodic_box [data_dir]
#
# The program reads the examples' data folder, unless data_dir is given,
# and writes to the folder output/ next to this file.

cmake_minimum_required(VERSION 3.16)

project(otswap_example_periodic_box LANGUAGES CXX)

find_package(otswap 0.4 CONFIG REQUIRED)

add_executable(periodic_box periodic_box.cpp)
target_link_libraries(periodic_box PRIVATE otswap::otswap)
target_compile_definitions(periodic_box PRIVATE
  OTSWAP_EXAMPLES_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data"
  OTSWAP_EXAMPLE_OUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/output")
