#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# Builds the native dependencies of the wheels from pinned releases, each
# checked against its SHA-256 before use, into one prefix:
#
#   zlib     static, position-independent: cfitsio's only required dependency
#   cfitsio  static, position-independent, without curl, bzip2, pthreads,
#            SSE, the network drivers, the utilities and the test programs
#   libomp   macOS only: a shared library per architecture, in omp-arm64/
#            and omp-x86_64/, which delocate copies into each wheel
#
# On macOS zlib and cfitsio are universal (arm64 and x86_64), and everything
# is built for MACOSX_DEPLOYMENT_TARGET, or 11.0 for arm64 if that is later.
# On Linux OpenMP is the compiler's libgomp, which auditwheel copies into
# each wheel.
#
# Usage: build_deps.sh PREFIX
#
# A prefix already built with the same versions is left alone.

set -euo pipefail

PREFIX=$1

ZLIB_VERSION=1.3.2
ZLIB_URL=https://zlib.net/zlib-${ZLIB_VERSION}.tar.gz
ZLIB_SHA256=bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16

CFITSIO_VERSION=4.7.0
CFITSIO_URL=https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${CFITSIO_VERSION}.tar.gz
CFITSIO_SHA256=ce573bbea8e75b429f8c3d3e86498741ba3dc9628a1530d2f65268397ad059e8

LLVM_VERSION=23.1.2
LLVM_URL=https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VERSION}/llvm-project-${LLVM_VERSION}.src.tar.xz
LLVM_SHA256=c98bbef08a2b4c2613cd50e9aa9ae7b69b1fe6c16b2c40373bc0ab6116fdf78a

STAMP="zlib ${ZLIB_VERSION}, cfitsio ${CFITSIO_VERSION}, libomp ${LLVM_VERSION}, macOS ${MACOSX_DEPLOYMENT_TARGET:-none}"

if [ -f "$PREFIX/.built" ] && [ "$(cat "$PREFIX/.built")" = "$STAMP" ]; then
  echo "build_deps: $PREFIX already holds $STAMP"
  exit 0
fi

mkdir -p "$PREFIX"
PREFIX=$(cd "$PREFIX" && pwd)
WORK="$PREFIX/src"
rm -rf "$WORK"
mkdir -p "$WORK"

sha256_of () {
  if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

# fetch URL SHA256: downloads into $WORK and refuses a mismatching file
fetch () {
  local file="$WORK/$(basename "$1")"
  curl --fail --silent --show-error --location --retry 3 --output "$file" "$1"
  local actual
  actual=$(sha256_of "$file")
  if [ "$actual" != "$2" ]; then
    echo "build_deps: $file has SHA-256 $actual, expected $2" >&2
    exit 1
  fi
  echo "$file"
}

CMAKE_COMMON=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
              -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib)
if [ "$(uname)" = Darwin ]; then
  CMAKE_COMMON+=(-DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:?}"
                 "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64")
fi

jobs=$( (nproc || sysctl -n hw.ncpu) 2> /dev/null )

# --------------------------------------------------------------------- zlib

tar -xzf "$(fetch "$ZLIB_URL" "$ZLIB_SHA256")" -C "$WORK"
cmake -S "$WORK/zlib-$ZLIB_VERSION" -B "$WORK/zlib-build" "${CMAKE_COMMON[@]}" \
  -DZLIB_BUILD_SHARED=OFF -DZLIB_BUILD_STATIC=ON -DZLIB_BUILD_TESTING=OFF
cmake --build "$WORK/zlib-build" --parallel "$jobs"
cmake --install "$WORK/zlib-build"

# ------------------------------------------------------------------ cfitsio

# The network drivers are compiled whenever gethostbyname and connect are
# found; presetting the results of those checks leaves them out.
tar -xzf "$(fetch "$CFITSIO_URL" "$CFITSIO_SHA256")" -C "$WORK"
cmake -S "$WORK/cfitsio-$CFITSIO_VERSION" -B "$WORK/cfitsio-build" "${CMAKE_COMMON[@]}" \
  -DBUILD_SHARED_LIBS=OFF -DUSE_CURL=OFF -DUSE_BZIP2=OFF -DUSE_PTHREADS=OFF \
  -DUSE_SSE2=OFF -DUSE_SSSE3=OFF -DTESTS=OFF -DUTILS=OFF -DITERPROGS=OFF \
  -DCMAKE_HAVE_GETHOSTBYNAME=OFF -DCMAKE_HAVE_CONNECT=OFF \
  -DZLIB_INCLUDE_DIR="$PREFIX/include" -DZLIB_LIBRARY="$PREFIX/lib/libz.a"
cmake --build "$WORK/cfitsio-build" --parallel "$jobs"
cmake --install "$WORK/cfitsio-build"

# ------------------------------------------------------------------- libomp

if [ "$(uname)" = Darwin ]; then
  llvm="llvm-project-$LLVM_VERSION.src"
  # Only the parts of the source tree the runtimes build of openmp reads.
  tar -xJf "$(fetch "$LLVM_URL" "$LLVM_SHA256")" -C "$WORK" \
    "$llvm/runtimes" "$llvm/openmp" "$llvm/cmake" "$llvm/third-party" \
    "$llvm/llvm/cmake" "$llvm/llvm/utils/lit" "$llvm/llvm/utils/llvm-lit"
  for arch in arm64 x86_64; do
    omp="$PREFIX/omp-$arch"
    # arm64 macOS starts at 11.0, the target cibuildwheel raises arm64
    # builds to, and the one libomp's availability checks need to see.
    target=$MACOSX_DEPLOYMENT_TARGET
    if [ "$arch" = arm64 ] && [ "${target%%.*}" -lt 11 ]; then target=11.0; fi
    cmake -S "$WORK/$llvm/runtimes" -B "$WORK/omp-build-$arch" -DLLVM_ENABLE_RUNTIMES=openmp \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$omp" -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_OSX_DEPLOYMENT_TARGET="$target" -DCMAKE_OSX_ARCHITECTURES="$arch" \
      -DLLVM_INCLUDE_TESTS=OFF -DLIBOMP_INSTALL_ALIASES=OFF -DLIBOMP_OMPT_SUPPORT=OFF \
      -DLIBOMP_OMPD_SUPPORT=OFF -DLIBOMP_USE_HWLOC=OFF -DOPENMP_ENABLE_OMPT_TOOLS=OFF
    cmake --build "$WORK/omp-build-$arch" --parallel "$jobs"
    cmake --install "$WORK/omp-build-$arch"
    # An absolute install name, so that delocate finds this library from the
    # extension that links it.
    install_name_tool -id "$omp/lib/libomp.dylib" "$omp/lib/libomp.dylib"
  done
fi

rm -rf "$WORK"
echo "$STAMP" > "$PREFIX/.built"
echo "build_deps: built $STAMP in $PREFIX"
