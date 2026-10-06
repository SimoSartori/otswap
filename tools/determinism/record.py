# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Record tests/determinism_hashes.json, the reference of
tests/test_determinism.py, from the otswap installed in this interpreter.

The cases are run with one thread and with four, each in a fresh
interpreter, and nothing is written unless the two agree. Run by hand, once,
when an intended change moves the outputs:

    python -B tools/determinism/record.py
"""

import datetime
import json
import platform
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parents[2] / "tests"
sys.path.insert(0, str(TESTS))

import test_determinism  # noqa: E402


def main():
    runs = {t: test_determinism.compute_with_threads(t) for t in test_determinism.THREADS}
    first = runs[test_determinism.THREADS[0]]
    for threads, hashes in runs.items():
        if hashes != first:
            differ = sorted(k for k in first if hashes.get(k) != first[k])
            sys.exit(f"not recorded: with {threads} threads these outputs differ: {differ}")

    import numpy

    record = {
        "recorded": {
            "date": datetime.date.today().isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": numpy.__version__,
            "threads": list(test_determinism.THREADS),
        },
        "hashes": first,
    }
    test_determinism.REFERENCE.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n")
    print(f"recorded {len(first)} hashes in {test_determinism.REFERENCE}")


if __name__ == "__main__":
    main()
