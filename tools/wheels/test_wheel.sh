#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# Runs the full Python test suite against an installed wheel, including the
# bit-identity tests: tests/cpp_reference is built for the wheel's
# architecture, from the same sources and against the same dependencies, and
# its absence is an error rather than a skip. Then the reproducibility test,
# tests/test_determinism.py, which compares the wheel's outputs with the
# hashes recorded in tests/determinism_hashes.json.
#
# Usage: test_wheel.sh PROJECT WHEEL
#
# The reference driver is built once per architecture, in
# PROJECT/build/wheel-deps/cpp-reference-ARCH.

set -euo pipefail

PROJECT=$1
WHEEL=$2
DEPS="$PROJECT/build/wheel-deps"

case "$WHEEL" in
  *_x86_64.whl)  arch=x86_64 ;;
  *_arm64.whl)   arch=arm64 ;;
  *_aarch64.whl) arch=aarch64 ;;
  *) echo "test_wheel: no architecture in $WHEEL" >&2; exit 1 ;;
esac

build="$DEPS/cpp-reference-$arch"
reference="$build/cpp_reference"

if [ ! -x "$reference" ]; then
  args=(-DCMAKE_BUILD_TYPE=Release -DOTSWAP_BUILD_TESTS=ON -DOTSWAP_BUILD_EXAMPLES=OFF
        -DCFITSIO_ROOT="$DEPS")
  if [ "$(uname)" = Darwin ]; then
    # the wheel's own deployment target, from its tag: macosx_10_13_x86_64
    target=$(basename "$WHEEL" | sed -E 's/.*-macosx_([0-9]+)_([0-9]+)_.*/\1.\2/')
    args+=(-DCMAKE_OSX_ARCHITECTURES="$arch" -DCMAKE_OSX_DEPLOYMENT_TARGET="$target"
           -DOpenMP_ROOT="$DEPS/omp-$arch")
  fi
  ZLIB_ROOT="$DEPS" cmake -S "$PROJECT" -B "$build" "${args[@]}"
  cmake --build "$build" --target cpp_reference --parallel
fi

OTSWAP_CPP_REFERENCE="$reference" python -m pytest "$PROJECT/tests/test_otswap.py" \
  "$PROJECT/tests/test_rsd.py" "$PROJECT/tests/test_determinism.py" -p no:cacheprovider -rs
