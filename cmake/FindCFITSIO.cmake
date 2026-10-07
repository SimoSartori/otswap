# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# FindCFITSIO
# -----------
#
# Finds the cfitsio library and defines the imported target CFITSIO::CFITSIO.
#
# Hints: CFITSIO_ROOT, as a CMake variable or in the environment, or
# CFITSIO_INCLUDE_DIR and CFITSIO_LIBRARY set directly.
#
# A static cfitsio needs zlib at every link, so the target then carries
# ZLIB::ZLIB, found with find_package(ZLIB) (hint: ZLIB_ROOT). A shared one
# resolves it itself.
#
# Result variables: CFITSIO_FOUND, CFITSIO_INCLUDE_DIRS, CFITSIO_LIBRARIES.

find_path(CFITSIO_INCLUDE_DIR fitsio.h
  HINTS ${CFITSIO_ROOT} ENV CFITSIO_ROOT
  PATH_SUFFIXES include include/cfitsio cfitsio)
find_library(CFITSIO_LIBRARY NAMES cfitsio
  HINTS ${CFITSIO_ROOT} ENV CFITSIO_ROOT
  PATH_SUFFIXES lib)

set(_cfitsio_static OFF)
if (CFITSIO_LIBRARY MATCHES "\\${CMAKE_STATIC_LIBRARY_SUFFIX}$")
  set(_cfitsio_static ON)
  find_package(ZLIB QUIET)
endif()

include(FindPackageHandleStandardArgs)
if (_cfitsio_static)
  find_package_handle_standard_args(CFITSIO
    REQUIRED_VARS CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR ZLIB_FOUND)
else()
  find_package_handle_standard_args(CFITSIO
    REQUIRED_VARS CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR)
endif()

if (CFITSIO_FOUND)
  set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})
  if (_cfitsio_static)
    list(APPEND CFITSIO_LIBRARIES ZLIB::ZLIB)
  endif()

  if (NOT TARGET CFITSIO::CFITSIO)
    add_library(CFITSIO::CFITSIO UNKNOWN IMPORTED)
    set_target_properties(CFITSIO::CFITSIO PROPERTIES
      IMPORTED_LOCATION "${CFITSIO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INCLUDE_DIR}")
    if (_cfitsio_static)
      set_target_properties(CFITSIO::CFITSIO PROPERTIES INTERFACE_LINK_LIBRARIES ZLIB::ZLIB)
    endif()
  endif()
endif()

mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)
unset(_cfitsio_static)
