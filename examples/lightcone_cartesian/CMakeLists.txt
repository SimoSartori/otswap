# Builds the example against an installed otswap, from this folder alone.
# Configure with -DCMAKE_PREFIX_PATH=<the prefix otswap was installed to>,
# and CFITSIO_ROOT (and ZLIB_ROOT for a static cfitsio) where cfitsio is not
# in a system path:
#
#   cmake -S examples/lightcone_cartesian -B build/examples/lightcone_cartesian -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build build/examples/lightcone_cartesian
#   build/examples/lightcone_cartesian/lightcone_cartesian [data_dir]
#
# The program reads the examples' data folder, unless data_dir is given,
# and writes to the folder output/ next to this file.

cmake_minimum_required(VERSION 3.16)

project(otswap_example_lightcone_cartesian LANGUAGES CXX)

find_package(otswap 0.4 CONFIG REQUIRED)

add_executable(lightcone_cartesian lightcone_cartesian.cpp)
target_link_libraries(lightcone_cartesian PRIVATE otswap::otswap)
target_compile_definitions(lightcone_cartesian PRIVATE
  OTSWAP_EXAMPLES_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data"
  OTSWAP_EXAMPLE_OUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/output")
