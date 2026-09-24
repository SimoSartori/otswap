#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# macOS, before each wheel build: points PREFIX/omp, the OpenMP_ROOT of the
# build, at the libomp built by build_deps.sh for the architecture named in
# ARCHFLAGS, which cibuildwheel sets per build.
#
# Usage: select_omp.sh PREFIX

set -euo pipefail

PREFIX=$1

case "${ARCHFLAGS:-}" in
  "-arch arm64")  arch=arm64 ;;
  "-arch x86_64") arch=x86_64 ;;
  *) echo "select_omp: ARCHFLAGS is '${ARCHFLAGS:-}'; one of arm64 or x86_64 is required" >&2
     exit 1 ;;
esac

if [ ! -f "$PREFIX/omp-$arch/lib/libomp.dylib" ]; then
  echo "select_omp: $PREFIX/omp-$arch/lib/libomp.dylib is missing; run build_deps.sh" >&2
  exit 1
fi

ln -sfn "omp-$arch" "$PREFIX/omp"
echo "select_omp: $PREFIX/omp -> omp-$arch"
