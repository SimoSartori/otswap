# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Bit-for-bit reproducibility across platforms.

Every case below hashes the bytes of one output array, and the hashes must
equal the ones recorded in tests/determinism_hashes.json: the same inputs
give the same bits on every platform the wheels are built for, with one
thread and with four.

The inputs are identical in bytes everywhere: they come from a fixed 64-bit
linear congruential generator in pure Python, mapped to doubles by exact
division by 2^53 and then by single IEEE operations, with no libm call, no
NumPy random generator and no file in the repository. The mask is written
to a temporary FITS file from integer rules.

Each thread count runs in its own interpreter, since OpenMP reads
OMP_NUM_THREADS once. Run as a script, this file prints the hashes as JSON;
tools/determinism/record.py records them.
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "determinism_hashes.json"
THREADS = (1, 4)


# --------------------------------------------------------------------------
# inputs

class Lcg:
    """Knuth's MMIX generator, x <- 6364136223846793005 x + 1442695040888963407
    mod 2^64; a double from the top 53 bits."""

    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFFFFFFFFFF

    def next(self):
        self.state = (6364136223846793005 * self.state + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        return self.state

    def unit(self):
        return (self.next() >> 11) / 9007199254740992.0

    def uniform(self, lo, hi):
        return lo + self.unit() * (hi - lo)


def box_points(n, size, seed):
    g = Lcg(seed)
    return np.array([[g.uniform(0.0, size) for _ in range(3)] for _ in range(n)])


def sky_points(n, seed):
    """Right ascension and declination in degrees, uniform in each, and
    redshift."""
    g = Lcg(seed)
    return np.array([[g.uniform(10.0, 40.0), g.uniform(-15.0, 15.0), g.uniform(0.3, 0.6)]
                     for _ in range(n)])


def write_mask(path, nside):
    """NSIDE nside, RING. Pixels are unobserved by an integer rule, about one
    in nine; observed ones hold one of 0.25, 0.5, 0.75, 1."""
    npix = 12 * nside * nside
    values = []
    for p in range(npix):
        h = (p * 2654435761) & 0xFFFFFFFF
        values.append(0.0 if (h >> 7) % 9 == 0 else 0.25 * (1 + (h >> 3) % 4))
    pixels = np.array(values, dtype=">f4")

    def card(key, value):
        if isinstance(value, bool):
            text = f"{'T' if value else 'F':>20}"
        elif isinstance(value, int):
            text = f"{value:>20}"
        else:
            text = f"'{value:<8}'"
        return f"{key:<8}= {text}".ljust(80)

    def block(cards):
        header = "".join(cards) + "END".ljust(80)
        return header.ljust(-(-len(header) // 2880) * 2880).encode("ascii")

    primary = block([card("SIMPLE", True), card("BITPIX", 8), card("NAXIS", 0),
                     card("EXTEND", True)])
    table = block([card("XTENSION", "BINTABLE"), card("BITPIX", 8), card("NAXIS", 2),
                   card("NAXIS1", 4096), card("NAXIS2", npix // 1024), card("PCOUNT", 0),
                   card("GCOUNT", 1), card("TFIELDS", 1), card("TTYPE1", "SIGNAL"),
                   card("TFORM1", "1024E"), card("PIXTYPE", "HEALPIX"),
                   card("ORDERING", "RING"), card("NSIDE", nside)])
    data = pixels.tobytes()
    data += b"\0" * (-len(data) % 2880)
    with open(path, "wb") as f:
        f.write(primary + table + data)
    return str(path)


# --------------------------------------------------------------------------
# cases

def digest(a):
    a = np.asarray(a)
    if a.dtype == np.bool_:
        a = a.astype("|u1")
    elif a.dtype.kind == "f":
        a = a.astype("<f8")
    else:
        a = a.astype("<u4")
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def compute():
    """The hash of every case, by name."""
    import otswap

    out = {}

    table = otswap.DistanceTable.flat(0.31, 0.68, w0=-0.9, wa=0.2, z_min=0.0, z_max=3.0,
                                      n_samples=20000)
    z = np.array([3.0 * i / 997 for i in range(998)])
    out["distance_table.distance_at"] = digest(table.distance_at(z))
    out["distance_table.growth_rate_at"] = digest(table.growth_rate_at(z))
    d = np.array([10.0 + 4000.0 * i / 499 for i in range(500)])
    out["distance_table.redshift_at"] = digest(table.redshift_at(d))

    lc_table = otswap.DistanceTable.flat(0.3, 0.7, z_min=0.0, z_max=1.5, n_samples=4000)
    tracers_sky = sky_points(1200, 11)
    randoms_sky = sky_points(4 * 1200, 12)
    out["to_cartesian"] = digest(otswap.to_cartesian(tracers_sky, lc_table, angle_unit="deg"))

    tracers = box_points(600, 100.0, 21)
    randoms = box_points(4 * 600, 100.0, 22)
    mps = 11.85                         # (100^3 / 600)^(1/3), written out: no libm
    for label, given in (("generated", None), ("given", randoms)):
        r = otswap.reconstruct_box(tracers, given, mps=mps, n_realizations=3, seed=12345)
        out[f"reconstruct_box.{label}.displacement"] = digest(r.displacement)
        out[f"reconstruct_box.{label}.mean_displacement"] = digest(r.mean_displacement)

    r = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, sky_area_deg2=900.0, n_bins=2,
                                     distances=lc_table, angle_unit="deg", n_realizations=3,
                                     seed=6789)
    out["reconstruct_lightcone.displacement"] = digest(r.displacement)
    out["reconstruct_lightcone.mean_displacement"] = digest(r.mean_displacement)

    with tempfile.TemporaryDirectory() as tmp:
        mask = otswap.Mask(write_mask(Path(tmp) / "mask.fits", 64))
        out["mask.sky_area_deg2"] = digest([mask.sky_area_deg2])

        ra = np.array([360.0 * i / 720 for i in range(721)])
        dec = np.array([-90.0 + 180.0 * j / 360 for j in range(361)])
        grid_ra, grid_dec = np.meshgrid(ra, dec)
        out["mask.allows"] = digest(mask.allows(grid_ra, grid_dec, angle_unit="deg"))

        for limit in (2, 0):
            otswap.reject_mask_crossings(r, mask, limit)
            out[f"reject_mask_crossings.{limit}.valid"] = digest(r.valid)
            out[f"reject_mask_crossings.{limit}.mean_displacement"] = digest(r.mean_displacement)

    return out


def compute_with_threads(threads):
    """compute() in a fresh interpreter with OMP_NUM_THREADS = threads."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    run = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve())], env=env,
                         capture_output=True, text=True, check=False)
    if run.returncode != 0:
        raise RuntimeError(f"the cases failed with {threads} threads:\n{run.stderr}")
    return json.loads(run.stdout)


# --------------------------------------------------------------------------
# tests

def test_hashes_do_not_depend_on_the_thread_count():
    one, four = (compute_with_threads(t) for t in THREADS)
    assert one == four


def test_hashes_match_the_reference():
    import pytest

    if not REFERENCE.exists():
        pytest.skip("no reference recorded; run tools/determinism/record.py")
    expected = json.loads(REFERENCE.read_text())["hashes"]
    for threads in THREADS:
        got = compute_with_threads(threads)
        differ = sorted(k for k in expected if got.get(k) != expected[k])
        assert set(got) == set(expected), "the cases differ from the recorded ones"
        assert not differ, f"with {threads} threads these outputs changed: {differ}"


if __name__ == "__main__":
    print(json.dumps(compute(), indent=1, sort_keys=True))
