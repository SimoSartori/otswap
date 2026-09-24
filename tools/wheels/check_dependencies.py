# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Fails unless every repaired wheel in a directory depends on nothing
outside itself but the allowed system libraries, and prints what each
binary in it links.

macOS: every load command of every binary must name a library under
/usr/lib or /System/Library, or one shipped in the wheel.
Linux: auditwheel show must find the wheel consistent with a manylinux tag
among its own, which admits no external library beyond that policy's.

Usage: python check_dependencies.py WHEEL_DIRECTORY
"""

import json
import os
import subprocess
import sys
import tempfile
import zipfile

SYSTEM_PREFIXES = ("/usr/lib/", "/System/Library/")


def check_macos(wheel):
    problems = []
    with tempfile.TemporaryDirectory() as root:
        with zipfile.ZipFile(wheel) as z:
            z.extractall(root)
        shipped = {}
        binaries = []
        for folder, _, files in os.walk(root):
            for name in files:
                path = os.path.join(folder, name)
                if name.endswith((".so", ".dylib")):
                    binaries.append(path)
                    shipped[name] = path
        for binary in binaries:
            out = subprocess.run(["otool", "-L", binary], check=True,
                                 capture_output=True, text=True).stdout.splitlines()
            print(os.path.relpath(binary, root))
            # the first line names the file; a dylib's second is its own id
            entries = [line.split(" (compatibility")[0].strip() for line in out[1:]]
            if binary.endswith(".dylib"):
                entries = entries[1:]
            for dep in entries:
                print("    " + dep)
                if dep.startswith(SYSTEM_PREFIXES):
                    continue
                if dep.startswith(("@loader_path/", "@rpath/")) and os.path.basename(dep) in shipped:
                    continue
                problems.append(f"{os.path.relpath(binary, root)} links {dep}")
    return problems


def check_linux(wheel):
    shown = subprocess.run(["auditwheel", "show", wheel], check=True,
                           capture_output=True, text=True).stdout
    print(shown)
    report = json.loads(subprocess.run(["auditwheel", "show", "--json", wheel], check=True,
                                       capture_output=True, text=True).stdout)
    tags = os.path.basename(wheel)[:-len(".whl")].split("-")[-1].split(".")
    overall = report["overall_tag"]
    if not overall.startswith("manylinux_") or overall not in tags:
        return [f"auditwheel finds {os.path.basename(wheel)} consistent with {overall}, "
                "which is not a manylinux tag of the wheel"]
    return []


def main():
    directory = sys.argv[1]
    wheels = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".whl")]
    if not wheels:
        sys.exit(f"check_dependencies: no wheel in {directory}")
    problems = []
    for wheel in wheels:
        print(f"== {os.path.basename(wheel)}")
        problems += check_macos(wheel) if sys.platform == "darwin" else check_linux(wheel)
    if problems:
        sys.exit("check_dependencies: " + "; ".join(problems))


if __name__ == "__main__":
    main()
