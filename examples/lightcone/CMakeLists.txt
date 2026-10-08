# Builds the example against an installed otswap, from this folder alone.
# Configure with -DCMAKE_PREFIX_PATH=<the prefix otswap was installed to>,
# and CFITSIO_ROOT (and ZLIB_ROOT for a static cfitsio) where cfitsio is not
# in a system path:
#
#   cmake -S examples/lightcone -B build/examples/lightcone -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build build/examples/lightcone
#   build/examples/lightcone/lightcone [data_dir]
#
# The program reads the examples' data folder, unless data_dir is given,
# and writes to the folder output/ next to this file.

cmake_minimum_required(VERSION 3.16)

project(otswap_example_lightcone LANGUAGES CXX)

find_package(otswap 0.1 CONFIG REQUIRED)

add_executable(lightcone lightcone.cpp)
target_link_libraries(lightcone PRIVATE otswap::otswap)
target_compile_definitions(lightcone PRIVATE
  OTSWAP_EXAMPLES_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data"
  OTSWAP_EXAMPLE_OUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/output")
