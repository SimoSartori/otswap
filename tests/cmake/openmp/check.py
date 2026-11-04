# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Checks the OpenMP runtime of an installed macOS wheel of otswap, for
tools/wheels/test_wheel.sh:

1. otswap/openmp/lib/libomp.dylib and otswap/openmp/include/omp.h exist;
   otswap/.dylibs/ does not.
2. _otswap loads libomp once, as @rpath/libomp.dylib, and its only rpath is
   @loader_path/openmp/lib.
3. The install name of the shipped libomp is @rpath/libomp.dylib.
4. The shipped libomp and _otswap are of the wheel's architecture only.
5. otswap.openmp_root() is the directory of 1.
6. The project beside this script builds a library that links that runtime
   through OpenMP_ROOT = otswap.openmp_root(), into a directory beside
   otswap's in site-packages, as a package installed next to otswap.
7. One process runs a parallel otswap call, then loads that library and runs
   a parallel region in it: it exits 0, without OMP: Error #15, which two
   runtimes in one process give.
8. In that process, DYLD_PRINT_LIBRARIES lists one libomp.dylib image, the
   shipped one. When dyld prints nothing (the variable removed on the way),
   this says the check was not made, and does not count it as passed.

Usage: python check.py ARCH BUILD_DIRECTORY
"""

import os
import subprocess
import sys

import otswap

RPATH_ID = "@rpath/libomp.dylib"
HERE = os.path.dirname(os.path.abspath(__file__))


def run(*args, **kwargs):
    return subprocess.run(args, check=True, capture_output=True, text=True, **kwargs).stdout


def load_commands(binary):
    lines = run("otool", "-L", binary).splitlines()[1:]
    return [line.split(" (compatibility")[0].strip() for line in lines]


def rpaths(binary):
    lines = run("otool", "-l", binary).splitlines()
    found = []
    for i, line in enumerate(lines):
        if line.strip() == "cmd LC_RPATH":
            path = lines[i + 2].strip()
            found.append(path.split("path ", 1)[1].split(" (offset")[0])
    return found


def dyld_images(stderr):
    images = []
    for line in stderr.splitlines():
        if not line.startswith("dyld[") or ": " not in line:
            continue
        rest = line.split(": ", 1)[1]
        if rest.startswith("<") and "> " in rest:
            rest = rest.split("> ", 1)[1]
        images.append(rest.strip())
    return images


def main():
    arch, build = sys.argv[1], sys.argv[2]
    failed = []

    def check(number, ok, what):
        print(f"check {number}: {'ok' if ok else 'FAILED'}: {what}")
        if not ok:
            failed.append(number)

    package = os.path.dirname(os.path.realpath(otswap.__file__))
    root = os.path.join(package, "openmp")
    library = os.path.join(root, "lib", "libomp.dylib")
    module = os.path.realpath(otswap._otswap.__file__)

    check(1, os.path.isfile(library) and os.path.isfile(os.path.join(root, "include", "omp.h"))
          and not os.path.exists(os.path.join(package, ".dylibs")),
          f"{root}/lib/libomp.dylib and include/omp.h exist, {package}/.dylibs does not")
    if not os.path.isfile(library):
        sys.exit(f"check.py: no {library}; checks 2-8 not made")

    omp = [dep for dep in load_commands(module) if "libomp" in dep]
    paths = rpaths(module)
    check(2, omp == [RPATH_ID] and paths == ["@loader_path/openmp/lib"],
          f"_otswap loads {omp}, rpaths {paths}")

    ident = run("otool", "-D", library).splitlines()[-1].strip()
    check(3, ident == RPATH_ID, f"install name of the shipped libomp: {ident}")

    archs = {name: run("lipo", "-archs", path).split() for name, path in
             (("libomp", library), ("_otswap", module))}
    check(4, all(a == [arch] for a in archs.values()), f"architectures {archs}, wheel {arch}")

    given = otswap.openmp_root()
    check(5, given is not None and os.path.realpath(given) == os.path.realpath(root),
          f"otswap.openmp_root() = {given}")

    output = os.path.join(os.path.dirname(package), "otswap_openmp_check")
    consumer = os.path.join(output, "openmp_consumer.so")
    configure = ["cmake", "-S", HERE, "-B", build, f"-DPython_EXECUTABLE={sys.executable}",
                 f"-DCMAKE_OSX_ARCHITECTURES={arch}", "-DCMAKE_BUILD_TYPE=Release",
                 f"-DOUTPUT_DIRECTORY={output}"]
    built = (subprocess.run(configure).returncode == 0
             and subprocess.run(["cmake", "--build", build]).returncode == 0
             and os.path.isfile(consumer))
    if built:
        consumer_omp = [dep for dep in load_commands(consumer) if "libomp" in dep]
        consumer_paths = rpaths(consumer)
        built = consumer_omp == [RPATH_ID] and consumer_paths == ["@loader_path/../otswap/openmp/lib"]
        check(6, built, f"{consumer} loads {consumer_omp}, rpaths {consumer_paths}")
    else:
        check(6, False, f"{consumer} not built")
    if not built:
        sys.exit(f"check.py: failed: {failed}; checks 7-8 not made")

    child = (
        "import ctypes, numpy as np, otswap\n"
        "tracers = np.random.default_rng(1).uniform(0.0, 100.0, (300, 3))\n"
        "result = otswap.reconstruct_box(tracers, mps=10.0,"
        " config=otswap.Config(n_realizations=4, seed=7, verbosity='silent'))\n"
        f"consumer = ctypes.CDLL({consumer!r})\n"
        "consumer.openmp_consumer_threads.restype = ctypes.c_int\n"
        "print('otswap realizations', result.n_realizations,"
        " 'consumer threads', consumer.openmp_consumer_threads())\n")
    env = dict(os.environ, OMP_NUM_THREADS="2", DYLD_PRINT_LIBRARIES="1")
    proc = subprocess.run([sys.executable, "-c", child], env=env, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    errors = [line for line in proc.stderr.splitlines() if not line.startswith("dyld[")]
    if errors:
        print("stderr, without dyld's lines:\n    " + "\n    ".join(errors))
    check(7, proc.returncode == 0 and "OMP: Error #15" not in proc.stderr,
          f"otswap and the consumer in one process: exit {proc.returncode}")

    images = dyld_images(proc.stderr)
    if not images:
        print("check 8: NOT MADE: DYLD_PRINT_LIBRARIES printed nothing; the runtimes loaded are "
              "unknown, only check 7 stands")
    else:
        loaded = [p for p in images if os.path.basename(p) == "libomp.dylib"]
        check(8, len(loaded) == 1 and os.path.realpath(loaded[0]) == os.path.realpath(library),
              f"{len(images)} images loaded; libomp.dylib: {loaded}")

    if failed:
        sys.exit(f"check.py: failed: {failed}")


if __name__ == "__main__":
    main()
