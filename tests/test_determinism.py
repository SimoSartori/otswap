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
import warnings
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


def mask_values(nside):
    """The pixel values of an NSIDE nside map, RING, as big-endian float32.
    Pixels are unobserved by an integer rule, about one in nine; observed
    ones hold one of 0.25, 0.5, 0.75, 1."""
    npix = 12 * nside * nside
    values = []
    for p in range(npix):
        h = (p * 2654435761) & 0xFFFFFFFF
        values.append(0.0 if (h >> 7) % 9 == 0 else 0.25 * (1 + (h >> 3) % 4))
    return np.array(values, dtype=">f4")


def write_mask(path, nside):
    """mask_values(nside) written as a FITS map."""
    pixels = mask_values(nside)
    npix = len(pixels)

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
    # The inverse, with a NaN row, the origin, and a point whose right
    # ascension folds from 2 pi to 0.
    xyz = np.vstack([otswap.to_cartesian(tracers_sky, lc_table, angle_unit="deg"),
                     [[np.nan, 1.0, 1.0], [0.0, 0.0, 0.0], [1000.0, -1e-17, 0.0]]])
    out["to_sky.deg"] = digest(otswap.to_sky(xyz, lc_table, angle_unit="deg"))
    out["to_sky.rad"] = digest(otswap.to_sky(xyz, lc_table, angle_unit="rad"))

    tracers = box_points(600, 100.0, 21)
    randoms = box_points(4 * 600, 100.0, 22)
    mps = 11.85                         # (100^3 / 600)^(1/3), written out: no libm
    for label, given in (("generated", None), ("given", randoms)):
        box = otswap.reconstruct_box(tracers, given, mps=mps, n_realizations=3, seed=12345,
                                     verbosity="silent")
        out[f"reconstruct_box.{label}.displacement"] = digest(box.displacement)
        out[f"reconstruct_box.{label}.mean_displacement"] = digest(box.mean_displacement)

    # The mps of the tracers' bounding box, computed by the library.
    computed = otswap.reconstruct_box(tracers, randoms, n_realizations=3, seed=12345,
                                      verbosity="silent")
    out["reconstruct_box.computed_mps.mps"] = digest([computed.mps])
    out["reconstruct_box.computed_mps.displacement"] = digest(computed.displacement)

    def lightcone():
        return otswap.reconstruct_lightcone(tracers_sky, randoms_sky, sky_area_deg2=900.0,
                                            n_bins=2, distances=lc_table, angle_unit="deg",
                                            n_realizations=3, seed=6789, verbosity="silent")

    r = lightcone()
    out["reconstruct_lightcone.displacement"] = digest(r.displacement)
    out["reconstruct_lightcone.mean_displacement"] = digest(r.mean_displacement)
    out["reconstruct_lightcone.lagrangian_sky"] = digest(r.lagrangian_sky)

    with tempfile.TemporaryDirectory() as tmp:
        mask = otswap.Mask(write_mask(Path(tmp) / "mask.fits", 64))
        out["mask.sky_area_deg2"] = digest([mask.sky_area_deg2])

        ra = np.array([360.0 * i / 720 for i in range(721)])
        dec = np.array([-90.0 + 180.0 * j / 360 for j in range(361)])
        grid_ra, grid_dec = np.meshgrid(ra, dec)
        out["mask.allows"] = digest(mask.allows(grid_ra, grid_dec, angle_unit="deg"))

        # The same map given in memory, converted from big-endian float32.
        in_memory = otswap.Mask.from_array(mask_values(64))
        area = digest([in_memory.sky_area_deg2])
        allows = digest(in_memory.allows(grid_ra, grid_dec, angle_unit="deg"))
        assert area == out["mask.sky_area_deg2"] and allows == out["mask.allows"], \
            "the mask from an array must equal the one read from FITS"
        out["mask.from_array.sky_area_deg2"] = area
        out["mask.from_array.allows"] = allows

        otswap.reject_mask_crossings(r, mask, 0)
        assert np.array_equal(r.lagrangian_sky, otswap.to_sky(r.lagrangian, lc_table, angle_unit="deg"),
                              equal_nan=True), "a filter must recompute lagrangian_sky"
        out["reject_mask_crossings.0.lagrangian_sky"] = digest(r.lagrangian_sky)
        out["reject_mask_crossings.0.valid"] = digest(r.valid)
        out["reject_mask_crossings.0.mean_displacement"] = digest(r.mean_displacement)

        # A limit of 2 on a finer mask, where arcs cross enough pixels for
        # it to reject some displacements and keep others.
        fine = otswap.Mask(write_mask(Path(tmp) / "fine.fits", 256))
        r2 = lightcone()
        otswap.reject_mask_crossings(r2, fine, 2)
        assert r2.valid.any() and not r2.valid.all(), \
            "the limit-2 case must reject some displacements and keep others"
        out["reject_mask_crossings.2.valid"] = digest(r2.valid)
        out["reject_mask_crossings.2.mean_displacement"] = digest(r2.mean_displacement)

        # A result rebuilt from its arrays has the reconstruction's means.
        rebuilt = otswap.Result.from_arrays(r2.displacement, r2.matched_random, r2.valid, r2.tracers,
                                            tracers_sky=tracers_sky, angle_unit="deg",
                                            distances=lc_table)
        assert digest(rebuilt.mean_displacement) == out["reject_mask_crossings.2.mean_displacement"] \
            and digest(rebuilt.lagrangian_sky) == digest(r2.lagrangian_sky), \
            "from_arrays must give the reconstruction's means"

        # The mask given to the reconstruction: tracers and randoms on
        # unobserved pixels left out before it, crossings rejected after it;
        # alone, and with a redshift cut.
        for label, cut in (("mask", None), ("mask_cut", (0.35, 0.55))):
            masked = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=2,
                                                  distances=lc_table, angle_unit="deg",
                                                  n_realizations=3, seed=6789, redshift_cut=cut,
                                                  verbosity="silent")
            out[f"reconstruct_lightcone.{label}.displacement"] = digest(masked.displacement)
            out[f"reconstruct_lightcone.{label}.valid"] = digest(masked.valid)
            out[f"reconstruct_lightcone.{label}.outside_mask"] = digest(masked.outside_mask)
            out[f"reconstruct_lightcone.{label}.lagrangian_sky"] = digest(masked.lagrangian_sky)

    # The redshift-space correction, on the lightcone filtered at limit 0,
    # whose tracers without a valid realization exercise the NaN paths, and
    # on the box reconstructed from given randoms. The bias table is
    # narrower than the tracers' redshifts, so b(z) is extrapolated.
    bias_z, bias = [0.35, 0.55], [1.3, 1.7]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", otswap.ExtrapolationWarning)

        positions = otswap.to_cartesian(tracers_sky, lc_table, angle_unit="deg")
        radial = otswap.radial_projection(r.mean_displacement, positions)
        along = otswap.axis_projection(box.mean_displacement, 2)
        out["rsd.projection.radial"] = digest(radial)
        out["rsd.projection.axis"] = digest(along)

        average, n_neighbours, n_realizations = otswap.neighbour_average(
            positions, radial, r.valid_realizations, sigma=10.0)
        out["rsd.neighbour_average"] = digest(average)
        out["rsd.neighbour_average.n_neighbours"] = digest(n_neighbours)
        out["rsd.neighbour_average.n_realizations"] = digest(n_realizations)
        out["rsd.neighbour_average.weighted"] = digest(otswap.neighbour_average(
            positions, radial, r.valid_realizations, sigma=10.0, weight_by_realizations=True).values)

        bias_table = otswap.BiasTable(bias_z, bias)
        factor = otswap.rsd_factor(tracers_sky[:, 2], lc_table, bias=bias_table)
        out["rsd.factor"] = digest(factor)
        out["rsd.factor_box"] = digest([otswap.rsd_factor_box(0.5, lc_table, bias=1.5)])

        out["rsd.shift.radial"] = digest(otswap.shift_radially(positions, factor * average))
        out["rsd.shift.axis"] = digest(otswap.shift_along_axis(tracers, 0.4 * along, 2))

        for sigma, weighted in ((0.0, False), (10.0, False), (10.0, True)):
            c = otswap.real_space_lightcone(r, distances=lc_table, bias=bias_table, sigma=sigma,
                                            weight_by_realizations=weighted, verbosity="silent")
            key = f"rsd.real_space_lightcone.{sigma:g}.{'weighted' if weighted else 'plain'}"
            out[key + ".positions"] = digest(c.sky)
            out[key + ".n_neighbours"] = digest(c.n_neighbours)
            out[key + ".uncorrected"] = digest(np.flatnonzero(c.status >= 2))
            out[key + ".cartesian"] = digest(c.cartesian)
            out[key + ".status"] = digest(c.status)

        cut = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, sky_area_deg2=900.0,
                                           n_bins=2, distances=lc_table, angle_unit="deg",
                                           n_realizations=3, seed=6789, redshift_cut=(0.35, 0.55),
                                           verbosity="silent")
        out["reconstruct_lightcone.cut.displacement"] = digest(cut.displacement)
        out["reconstruct_lightcone.cut.outside_redshift_cut"] = digest(cut.outside_redshift_cut)
        out["reconstruct_lightcone.cut.lagrangian_sky"] = digest(cut.lagrangian_sky)
        c = otswap.real_space_lightcone(cut, distances=lc_table, bias=bias_table, sigma=10.0,
                                        verbosity="silent")
        out["rsd.real_space_lightcone.cut.positions"] = digest(c.sky)

        for sigma, weighted in ((0.0, False), (8.0, False), (8.0, True)):
            c = otswap.real_space_box(box, axis=2, redshift=0.5, distances=lc_table, bias=1.5,
                                      sigma=sigma, weight_by_realizations=weighted,
                                      verbosity="silent")
            key = f"rsd.real_space_box.{sigma:g}.{'weighted' if weighted else 'plain'}"
            out[key + ".positions"] = digest(c.cartesian)
            out[key + ".n_neighbours"] = digest(c.n_neighbours)

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
        differ = sorted(k for k in expected if k in got and got[k] != expected[k])
        assert not differ, f"with {threads} threads these outputs changed: {differ}"
        assert set(got) == set(expected), (
            f"the cases differ from the recorded ones: new {sorted(set(got) - set(expected))}, "
            f"gone {sorted(set(expected) - set(got))}")


if __name__ == "__main__":
    print(json.dumps(compute(), indent=1, sort_keys=True))
