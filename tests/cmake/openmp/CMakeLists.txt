# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# A package beside otswap that links the OpenMP runtime of otswap's macOS
# wheel, as "Linking the same OpenMP runtime" (docs/python.rst) tells one to:
# OpenMP_ROOT from otswap.openmp_root(), and the rpath of a package installed
# next to otswap. Built by check.py, which gives Python_EXECUTABLE, the
# architecture and OUTPUT_DIRECTORY, a directory beside otswap's.

cmake_minimum_required(VERSION 3.16)

project(otswap_openmp_consumer LANGUAGES CXX)

execute_process(COMMAND "${Python_EXECUTABLE}" -c "import otswap; print(otswap.openmp_root() or '')"
  OUTPUT_VARIABLE OpenMP_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE status)
if (NOT status EQUAL 0 OR NOT OpenMP_ROOT)
  message(FATAL_ERROR "otswap.openmp_root() gave no directory (exit ${status}).")
endif()
message(STATUS "OpenMP_ROOT from otswap.openmp_root(): ${OpenMP_ROOT}")

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(openmp_consumer MODULE consumer.cpp)
target_link_libraries(openmp_consumer PRIVATE OpenMP::OpenMP_CXX)
set_target_properties(openmp_consumer PROPERTIES
  PREFIX ""
  BUILD_WITH_INSTALL_RPATH ON
  INSTALL_RPATH "@loader_path/../otswap/openmp/lib"
  LIBRARY_OUTPUT_DIRECTORY "${OUTPUT_DIRECTORY}")
