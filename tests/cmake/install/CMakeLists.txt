# Finds otswap as an installed package. Configure with
# -DCMAKE_PREFIX_PATH=<the prefix otswap was installed to>, and CFITSIO_ROOT
# (and ZLIB_ROOT for a static cfitsio) where cfitsio is not in a system path.

cmake_minimum_required(VERSION 3.16)

project(otswap_install_consumer LANGUAGES CXX)

find_package(otswap 0.1 CONFIG REQUIRED)

add_executable(consumer ../consumer.cpp)
target_link_libraries(consumer PRIVATE otswap::otswap)
