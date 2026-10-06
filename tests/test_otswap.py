"""Tests for the Python bindings.

They cover every function and property of the package, and what only exists
at the language boundary: the conversion of array-like inputs and angle
units, the read-only views on a result and their lifetime, the in-place
update made by reject_mask_crossings, the release of the GIL and the mapping
of every invalid input onto otswap.Error.

The bit-identity tests compare the bindings with the C++ library run on the
same inputs by tests/cpp_reference.cpp. They need the path of that executable
in the environment variable OTSWAP_CPP_REFERENCE, and are skipped without it.
"""

import gc
import math
import os
import subprocess
import threading
import time

import numpy as np
import pytest

import otswap


# --------------------------------------------------------------------------
# helpers

SPHERE_DEG2 = 4 * math.pi * (180 / math.pi) ** 2


def fits_card(key, value):
    """One 80-character FITS header card."""
    if isinstance(value, bool):
        text = f"{'T' if value else 'F':>20}"
    elif isinstance(value, int):
        text = f"{value:>20}"
    else:
        text = f"'{value:<8}'"
    return f"{key:<8}= {text}".ljust(80)


def fits_block(cards):
    header = "".join(cards) + "END".ljust(80)
    return header.ljust(-(-len(header) // 2880) * 2880).encode("ascii")


def write_mask(path, pixels, nside, ordering="RING"):
    """A HEALPix map as written by healpy: one BINTABLE column of 1024
    floats per row in HDU 2, with NSIDE and ORDERING in its header."""
    pixels = np.asarray(pixels, dtype=">f4")
    assert pixels.size == 12 * nside * nside and pixels.size % 1024 == 0
    nrows = pixels.size // 1024
    primary = fits_block([fits_card("SIMPLE", True), fits_card("BITPIX", 8),
                          fits_card("NAXIS", 0), fits_card("EXTEND", True)])
    table = fits_block([
        fits_card("XTENSION", "BINTABLE"), fits_card("BITPIX", 8),
        fits_card("NAXIS", 2), fits_card("NAXIS1", 4096), fits_card("NAXIS2", nrows),
        fits_card("PCOUNT", 0), fits_card("GCOUNT", 1), fits_card("TFIELDS", 1),
        fits_card("TTYPE1", "SIGNAL"), fits_card("TFORM1", "1024E"),
        fits_card("PIXTYPE", "HEALPIX"), fits_card("ORDERING", ordering),
        fits_card("NSIDE", nside), fits_card("FIRSTPIX", 0),
        fits_card("LASTPIX", pixels.size - 1)])
    data = pixels.tobytes()
    data += b"\0" * (-len(data) % 2880)
    with open(path, "wb") as f:
        f.write(primary + table + data)
    return str(path)


# NSIDE 16, RING: rings 28 to 36 are pixels [1248, 1824), the band around the
# equator where every direction with |sin(dec)| < 0.145 falls, and none with
# |sin(dec)| > 0.21.
BAND_NSIDE = 16
BAND_PIXELS = (1248, 1824)


@pytest.fixture(scope="module")
def band_mask(tmp_path_factory):
    pixels = np.ones(12 * BAND_NSIDE ** 2)
    pixels[BAND_PIXELS[0]:BAND_PIXELS[1]] = 0
    return write_mask(tmp_path_factory.mktemp("mask") / "band.fits", pixels, BAND_NSIDE)


@pytest.fixture(scope="module")
def sieve_mask(tmp_path_factory):
    """NSIDE 256, every seventh pixel unobserved: most displacements cross
    one."""
    pixels = np.ones(12 * 256 ** 2)
    pixels[::7] = 0
    return write_mask(tmp_path_factory.mktemp("mask") / "sieve.fits", pixels, 256)


@pytest.fixture(scope="module")
def full_mask32(tmp_path_factory):
    return write_mask(tmp_path_factory.mktemp("mask") / "full32.fits",
                      np.ones(12 * 32 ** 2), 32)


def box_catalogue(n_side=8, spacing=10.0, seed=1):
    """A jittered lattice of n_side^3 points, mps = spacing."""
    rng = np.random.default_rng(seed)
    g = (np.arange(n_side) + 0.5) * spacing
    lattice = np.stack(np.meshgrid(g, g, g, indexing="ij"), axis=-1).reshape(-1, 3)
    return lattice + rng.uniform(-0.3, 0.3, lattice.shape) * spacing


TABLE_ARGS = dict(omega_m=0.3, h=0.7, z_min=0.0, z_max=1.5, n_samples=4000)
RA_DEG = (10.0, 40.0)
DEC_DEG = (-15.0, 15.0)
Z_RANGE = (0.3, 0.6)


def patch_area_deg2():
    ra = np.deg2rad(RA_DEG)
    dec = np.deg2rad(DEC_DEG)
    return (ra[1] - ra[0]) * (np.sin(dec[1]) - np.sin(dec[0])) * (180 / math.pi) ** 2


def sky_catalogue(table, n, seed):
    """n directions uniform in the patch, at distances uniform in comoving
    volume between the redshifts of Z_RANGE; angles in degrees."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(*RA_DEG, n)
    sin_dec = rng.uniform(*np.sin(np.deg2rad(DEC_DEG)), n)
    d0, d1 = table.distance_at(list(Z_RANGE))
    d = np.cbrt(rng.uniform(d0 ** 3, d1 ** 3, n))
    return np.column_stack([ra, np.rad2deg(np.arcsin(sin_dec)), table.redshift_at(d)])


@pytest.fixture(scope="module")
def table():
    return otswap.DistanceTable.flat(**TABLE_ARGS)


@pytest.fixture(scope="module")
def lightcone(table):
    """Sky coordinates in degrees of 1600 tracers and 4800 randoms."""
    return sky_catalogue(table, 1600, 11), sky_catalogue(table, 4800, 12)


def lightcone_result(table, lightcone, **kwargs):
    tracers_sky, randoms_sky = lightcone
    options = dict(sky_area_deg2=patch_area_deg2(), n_bins=3, distances=table,
                   angle_unit="deg", n_realizations=2, seed=5)
    options.update(kwargs)
    return otswap.reconstruct_lightcone(tracers_sky, randoms_sky, **options)


def same_bytes(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def result_arrays(result):
    return {name: getattr(result, name) for name in
            ("displacement", "matched_random", "valid", "valid_realizations", "mean_displacement")}


# --------------------------------------------------------------------------
# the package


def test_exports():
    assert sorted(otswap.__all__) == sorted([
        "AngleUnit", "DistanceTable", "Error", "Mask", "Result", "reconstruct_box",
        "reconstruct_lightcone", "reject_mask_crossings", "to_cartesian"])
    for name in otswap.__all__:
        assert hasattr(otswap, name)
    assert otswap.AngleUnit.__args__ == ("deg", "rad")


def test_error_is_a_runtime_error():
    assert issubclass(otswap.Error, RuntimeError)
    with pytest.raises(RuntimeError):
        otswap.DistanceTable.flat(-1.0, 0.7)


def test_stub_and_marker_are_installed():
    here = os.path.dirname(otswap.__file__)
    assert os.path.isfile(os.path.join(here, "__init__.pyi"))
    assert os.path.isfile(os.path.join(here, "py.typed"))


def test_classes_have_no_public_constructor_where_the_spec_gives_none():
    with pytest.raises(TypeError):
        otswap.DistanceTable()
    with pytest.raises(TypeError):
        otswap.Result()


# --------------------------------------------------------------------------
# DistanceTable


def test_flat_defaults():
    default = otswap.DistanceTable.flat(0.3, 0.7)
    explicit = otswap.DistanceTable.flat(0.3, 0.7, w0=-1.0, wa=0.0, z_min=0.0,
                                         z_max=10.0, n_samples=50_000)
    z = np.linspace(0, 10, 101)
    assert same_bytes(default.distance_at(z), explicit.distance_at(z))
    assert default.min_redshift == 0.0 and default.max_redshift == 10.0
    assert default.has_growth_rate


def test_flat_properties_and_lookups(table):
    assert table.min_redshift == 0.0
    assert table.max_redshift == 1.5
    assert table.has_growth_rate is True
    z = np.linspace(0.01, 1.49, 50)
    d = table.distance_at(z)
    assert d.dtype == np.float64 and d.shape == z.shape
    assert np.all(np.diff(d) > 0)
    np.testing.assert_allclose(table.redshift_at(d), z, rtol=1e-10)
    f = table.growth_rate_at(z)
    assert np.all((f > 0.4) & (f < 1.0)) and np.all(np.diff(f) > 0)
    # at low z, D = c z / H0, which is 2997.92458 z in Mpc/h
    assert table.distance_at(1e-4) == pytest.approx(299792.458 / 100 * 1e-4, rel=1e-3)


def test_lookups_keep_the_shape_and_accept_array_likes(table):
    assert table.distance_at(0.5).shape == ()
    assert table.distance_at([[0.1, 0.2], [0.3, 0.4]]).shape == (2, 2)
    assert table.distance_at(np.array([0.1, 0.2], dtype=np.float32)).dtype == np.float64
    strided = np.linspace(0.1, 1.0, 20)[::3]
    assert same_bytes(table.distance_at(strided), table.distance_at(strided.copy()))
    assert table.redshift_at([]).shape == (0,)


def test_from_table():
    z = np.linspace(0.0, 2.0, 201)
    d = 3000.0 * z
    f = 0.5 + 0.1 * z
    t = otswap.DistanceTable.from_table(z, d, f)
    assert t.has_growth_rate and t.min_redshift == 0.0 and t.max_redshift == 2.0
    assert t.distance_at(1.0) == pytest.approx(3000.0)
    assert t.redshift_at(1500.0) == pytest.approx(0.5)
    assert t.growth_rate_at(1.0) == pytest.approx(0.6)
    plain = otswap.DistanceTable.from_table(list(z), tuple(d))
    assert not plain.has_growth_rate
    with pytest.raises(otswap.Error):
        plain.growth_rate_at(0.5)


@pytest.mark.parametrize("call", [
    lambda: otswap.DistanceTable.flat(-0.1, 0.7),
    lambda: otswap.DistanceTable.flat(0.3, 0.7, z_min=2.0, z_max=1.0),
    lambda: otswap.DistanceTable.flat(0.3, 0.7, n_samples=1),
    lambda: otswap.DistanceTable.flat(0.3, 0.7, n_samples=-5),
    lambda: otswap.DistanceTable.flat(0.3, 0.7, n_samples=2.5),
    lambda: otswap.DistanceTable.flat("0.3", 0.7),
    lambda: otswap.DistanceTable.from_table([0.0, 1.0, 0.5], [0.0, 1.0, 2.0]),
    lambda: otswap.DistanceTable.from_table([0.0], [0.0]),
    lambda: otswap.DistanceTable.from_table([[0.0, 1.0]], [[0.0, 1.0]]),
    lambda: otswap.DistanceTable.from_table(["a", "b"], [0.0, 1.0]),
    lambda: otswap.DistanceTable.from_table([0.0, 1.0], [0.0, 1.0], [0.5]),
])
def test_distance_table_errors(call):
    with pytest.raises(otswap.Error):
        call()


def test_lookups_raise_outside_the_table(table):
    with pytest.raises(otswap.Error):
        table.distance_at([0.5, 1.6])
    with pytest.raises(otswap.Error):
        table.redshift_at(-1.0)
    with pytest.raises(otswap.Error):
        table.growth_rate_at(2.0)
    with pytest.raises(otswap.Error):
        table.distance_at("x")


# --------------------------------------------------------------------------
# to_cartesian


def test_to_cartesian_round_trips(table, lightcone):
    sky = lightcone[0]
    xyz = otswap.to_cartesian(sky, table, angle_unit="deg")
    assert xyz.shape == sky.shape and xyz.dtype == np.float64
    r = np.linalg.norm(xyz, axis=1)
    np.testing.assert_allclose(r, table.distance_at(sky[:, 2]), rtol=1e-12)
    np.testing.assert_allclose(table.redshift_at(r), sky[:, 2], rtol=1e-9)
    ra = np.rad2deg(np.arctan2(xyz[:, 1], xyz[:, 0])) % 360
    dec = np.rad2deg(np.arcsin(xyz[:, 2] / r))
    np.testing.assert_allclose(ra, sky[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(dec, sky[:, 1], rtol=0, atol=1e-10)


def test_to_cartesian_deg_and_rad_are_identical(table, lightcone):
    sky = lightcone[0]
    rad = np.column_stack([np.deg2rad(sky[:, 0]), np.deg2rad(sky[:, 1]), sky[:, 2]])
    assert same_bytes(otswap.to_cartesian(sky, table, angle_unit="deg"),
                      otswap.to_cartesian(rad, table, angle_unit="rad"))


def test_to_cartesian_accepts_array_likes(table):
    sky = [[20.0, 5.0, 0.4], [30.0, -5.0, 0.5]]
    expected = otswap.to_cartesian(np.array(sky), table, angle_unit="deg")
    assert same_bytes(otswap.to_cartesian(sky, table, angle_unit="deg"), expected)
    fortran = np.asfortranarray(np.array(sky, dtype=np.float32))
    assert otswap.to_cartesian(fortran, table, angle_unit="deg").shape == (2, 3)
    assert otswap.to_cartesian(np.empty((0, 3)), table, angle_unit="rad").shape == (0, 3)


@pytest.mark.parametrize("sky, kwargs", [
    (np.zeros((4, 2)), dict(angle_unit="deg")),
    (np.zeros(3), dict(angle_unit="deg")),
    ([[20.0, 5.0, 0.4]], dict(angle_unit="degrees")),
    ([[20.0, 5.0, 0.4]], dict(angle_unit=None)),
    ([[20.0, 5.0, 3.0]], dict(angle_unit="deg")),
    ([[20.0, 5.0, np.nan]], dict(angle_unit="deg")),
    ([[20.0, 5.0, "a"]], dict(angle_unit="deg")),
])
def test_to_cartesian_errors(table, sky, kwargs):
    with pytest.raises(otswap.Error):
        otswap.to_cartesian(sky, table, **kwargs)


def test_to_cartesian_needs_a_table_and_a_keyword_unit(table):
    with pytest.raises(otswap.Error):
        otswap.to_cartesian([[20.0, 5.0, 0.4]], "table", angle_unit="deg")
    with pytest.raises(TypeError):
        otswap.to_cartesian([[20.0, 5.0, 0.4]], table)
    with pytest.raises(TypeError):
        otswap.to_cartesian([[20.0, 5.0, 0.4]], table, "deg")


# --------------------------------------------------------------------------
# Mask


def test_mask_properties(band_mask):
    mask = otswap.Mask(band_mask)
    assert mask.nside == BAND_NSIDE
    npix = 12 * BAND_NSIDE ** 2
    observed = npix - (BAND_PIXELS[1] - BAND_PIXELS[0])
    assert mask.sky_area_deg2 == pytest.approx(SPHERE_DEG2 * observed / npix, rel=1e-12)


def test_mask_allows(band_mask):
    mask = otswap.Mask(band_mask)
    ra = np.linspace(0, 359, 360)
    for sin_dec, expected in [(0.0, False), (0.14, False), (-0.14, False),
                              (0.22, True), (-0.22, True), (0.9, True)]:
        dec = np.full_like(ra, np.rad2deg(np.arcsin(sin_dec)))
        allowed = mask.allows(ra, dec, angle_unit="deg")
        assert allowed.dtype == np.bool_ and allowed.shape == ra.shape
        assert np.all(allowed == expected)
        assert same_bytes(mask.allows(np.deg2rad(ra), np.deg2rad(dec), angle_unit="rad"), allowed)


def test_mask_allows_broadcasts(band_mask):
    mask = otswap.Mask(band_mask)
    grid = mask.allows(np.array([0.0, 90.0, 180.0])[:, None], [0.0, 45.0], angle_unit="deg")
    assert grid.shape == (3, 2)
    assert not grid[:, 0].any() and grid[:, 1].all()
    assert mask.allows(10.0, 60.0, angle_unit="deg").shape == ()
    assert mask.allows(10.0, 60.0, angle_unit="deg")


@pytest.mark.parametrize("call", [
    lambda m: m.allows([1.0, 2.0], [1.0, 2.0, 3.0], angle_unit="deg"),
    lambda m: m.allows([1.0], [1.0], angle_unit="radians"),
    lambda m: m.allows(["a"], [1.0], angle_unit="deg"),
    lambda m: m.allows([1.0], [np.nan], angle_unit="deg"),
    lambda m: m.allows([np.inf], [1.0], angle_unit="deg"),
    lambda m: m.allows([1.0], [np.nextafter(90.0, 91.0)], angle_unit="deg"),
    lambda m: m.allows([1.0], [np.nextafter(np.pi / 2, 2.0)], angle_unit="rad"),
])
def test_mask_allows_errors(band_mask, call):
    with pytest.raises(otswap.Error):
        call(otswap.Mask(band_mask))


def test_mask_errors(tmp_path):
    with pytest.raises(otswap.Error):
        otswap.Mask(str(tmp_path / "missing.fits"))
    with pytest.raises(otswap.Error):
        otswap.Mask(42)


def test_mask_observed_above_zero(tmp_path):
    values = [0.5, 1e-30, 2.0, 1.0, 0.0, -0.0, -1.0, -1.6375e30, np.nan, -np.inf, np.inf]
    observed = [True, True, True, True, False, False, False, False, False, False, True]
    pixels = np.ones(12 * 16 ** 2)
    pixels[:len(values)] = values
    mask = otswap.Mask(write_mask(tmp_path / "values.fits", pixels, 16))
    # the centres of RING pixels 0 to 11: the first two rings of NSIDE 16
    centres = []
    for first, count in [(0, 4), (4, 8)]:
        z = 1 - (first // 4 + 1) ** 2 / (3 * 16 ** 2)
        for k in range(count):
            centres.append((360 * (k + 0.5) / count, math.degrees(math.asin(z))))
    ra, dec = np.array(centres[:len(values)]).T
    assert mask.allows(ra, dec, angle_unit="deg").tolist() == observed
    assert mask.sky_area_deg2 == pytest.approx(SPHERE_DEG2 * (pixels.size - 6) / pixels.size,
                                               rel=1e-12)


def test_mask_nested_ordering(tmp_path):
    pixels = np.ones(12 * 16 ** 2)
    pixels[:1024] = 0
    ring = otswap.Mask(write_mask(tmp_path / "ring.fits", pixels, 16, "RING"))
    nested = otswap.Mask(write_mask(tmp_path / "nested.fits", pixels, 16, "NESTED"))
    assert ring.sky_area_deg2 == nested.sky_area_deg2
    ra, dec = np.meshgrid(np.linspace(0, 359, 90), np.linspace(-89, 89, 45))
    assert (ring.allows(ra, dec, angle_unit="deg") != nested.allows(ra, dec, angle_unit="deg")).any()


# --------------------------------------------------------------------------
# reconstruct_box and Result


def test_reconstruct_box_result():
    tracers = box_catalogue()
    result = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=3, seed=7)
    assert isinstance(result, otswap.Result)
    n = len(tracers)
    assert result.n_objects == n and result.n_realizations == 3
    assert result.displacement.shape == (3, n, 3) and result.displacement.dtype == np.float64
    assert result.matched_random.shape == (3, n, 3) and result.matched_random.dtype == np.float64
    assert result.valid.shape == (3, n) and result.valid.dtype == np.bool_
    assert result.valid.all()
    assert result.valid_realizations.shape == (n,) and result.valid_realizations.dtype == np.uint32
    assert np.all(result.valid_realizations == 3)
    assert result.mean_displacement.shape == (n, 3) and result.mean_displacement.dtype == np.float64
    np.testing.assert_allclose(result.mean_displacement, result.displacement.mean(axis=0),
                               rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result.matched_random, tracers + result.displacement,
                               rtol=1e-12, atol=1e-9)
    assert result.filtered_nside == 0
    # the generated randoms fill the bounding box of the tracers
    lo, hi = tracers.min(axis=0), tracers.max(axis=0)
    assert np.all(result.matched_random >= lo) and np.all(result.matched_random <= hi)
    # each realization matches every tracer to a distinct random
    for rec in range(3):
        assert len(np.unique(result.matched_random[rec], axis=0)) == n


def test_reconstruct_box_with_randoms_uses_them():
    tracers = box_catalogue()
    rng = np.random.default_rng(3)
    randoms = rng.uniform(0, 80, (2 * len(tracers), 3))
    result = otswap.reconstruct_box(tracers, randoms, mps=10.0, n_realizations=2, seed=4)
    matched = result.matched_random.reshape(-1, 3)
    pool = {tuple(p) for p in randoms}
    assert all(tuple(p) in pool for p in matched)
    # the two realizations use disjoint subsets
    first = {tuple(p) for p in result.matched_random[0]}
    assert not first & {tuple(p) for p in result.matched_random[1]}


def test_reconstruct_box_options_and_array_likes():
    tracers = box_catalogue()
    base = otswap.reconstruct_box(tracers, mps=10.0, seed=9)
    again = otswap.reconstruct_box(tracers.tolist(), None, mps=10, seed=np.int64(9))
    assert same_bytes(base.displacement, again.displacement)
    fortran = otswap.reconstruct_box(np.asfortranarray(tracers), mps=10.0, seed=9)
    assert same_bytes(base.displacement, fortran.displacement)
    other_seed = otswap.reconstruct_box(tracers, mps=10.0, seed=10)
    assert not same_bytes(base.displacement, other_seed.displacement)
    loose = otswap.reconstruct_box(tracers, mps=10.0, seed=9, convergence=0.5)
    assert loose.displacement.shape == base.displacement.shape


@pytest.mark.parametrize("with_randoms", [False, True])
def test_box_cell_size_affects_speed_only(with_randoms):
    tracers = box_catalogue()
    randoms = None
    if with_randoms:
        randoms = np.random.default_rng(6).uniform(0, 80, (2 * len(tracers), 3))
    results = [otswap.reconstruct_box(tracers, randoms, mps=10.0, n_realizations=2, seed=9,
                                      cell_size=cell_size) for cell_size in (1.0, 2.0, 4.0, 8.0)]
    for result in results[1:]:
        for name, array in result_arrays(result).items():
            assert same_bytes(array, getattr(results[0], name)), name


def test_reconstruct_box_needs_keywords():
    tracers = box_catalogue()
    with pytest.raises(TypeError):
        otswap.reconstruct_box(tracers)
    with pytest.raises(TypeError):
        otswap.reconstruct_box(tracers, None, 10.0)


@pytest.mark.parametrize("args, kwargs", [
    ((np.zeros((10, 2)),), dict(mps=1.0)),
    ((np.zeros(30),), dict(mps=1.0)),
    ((np.empty((0, 3)),), dict(mps=1.0)),
    ((box_catalogue(), np.empty((0, 3))), dict(mps=10.0)),
    ((box_catalogue(), []), dict(mps=10.0)),
    ((box_catalogue(), np.zeros((10, 3))), dict(mps=10.0)),
    ((box_catalogue(), np.zeros((600, 4))), dict(mps=10.0)),
    ((box_catalogue(),), dict(mps=0.0)),
    ((box_catalogue(),), dict(mps="ten")),
    ((box_catalogue(),), dict(mps=10.0, n_realizations=0)),
    ((box_catalogue(),), dict(mps=10.0, n_realizations=-1)),
    ((box_catalogue(),), dict(mps=10.0, seed=-1)),
    ((box_catalogue(),), dict(mps=10.0, seed=2 ** 32)),
    ((box_catalogue(),), dict(mps=10.0, seed=1.5)),
    ((box_catalogue(),), dict(mps=10.0, convergence=-0.1)),
    ((box_catalogue(),), dict(mps=10.0, convergence=None)),
    ((box_catalogue(),), dict(mps=10.0, cell_size=0.0)),
    (([[1.0, 2.0, np.nan]] * 20,), dict(mps=1.0)),
    (([["a", "b", "c"]],), dict(mps=1.0)),
])
def test_reconstruct_box_errors(args, kwargs):
    with pytest.raises(otswap.Error):
        otswap.reconstruct_box(*args, **kwargs)


def test_result_views_are_read_only():
    result = otswap.reconstruct_box(box_catalogue(), mps=10.0, n_realizations=2, seed=1)
    for name, array in result_arrays(result).items():
        assert not array.flags.writeable, name
        with pytest.raises(ValueError):
            array[...] = 0
        with pytest.raises(ValueError):
            array.flags.writeable = True
        writable = array.copy()
        writable[...] = 0
    for name in ("n_objects", "n_realizations", "displacement", "valid", "filtered_nside"):
        with pytest.raises(AttributeError):
            setattr(result, name, 0)


def test_result_views_share_the_result_memory():
    result = otswap.reconstruct_box(box_catalogue(), mps=10.0, n_realizations=2, seed=1)
    for name, array in result_arrays(result).items():
        assert np.shares_memory(array, getattr(result, name)), name


def test_result_views_outlive_the_result():
    result = otswap.reconstruct_box(box_catalogue(), mps=10.0, n_realizations=2, seed=1)
    views = result_arrays(result)
    copies = {name: array.copy() for name, array in views.items()}
    del result
    gc.collect()
    # allocate and free memory of the same size, so that freed memory would
    # be overwritten
    for _ in range(3):
        junk = otswap.reconstruct_box(box_catalogue(seed=2), mps=10.0, n_realizations=2, seed=2)
        del junk
        gc.collect()
    for name, array in views.items():
        assert same_bytes(array, copies[name]), name


# --------------------------------------------------------------------------
# reconstruct_lightcone


def test_reconstruct_lightcone(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    result = lightcone_result(table, lightcone)
    n = len(tracers_sky)
    assert result.n_objects == n and result.n_realizations == 2
    assert result.displacement.shape == (2, n, 3)
    assert result.valid.all() and np.all(result.valid_realizations == 2)
    tracers = otswap.to_cartesian(tracers_sky, table, angle_unit="deg")
    np.testing.assert_allclose(result.matched_random, tracers + result.displacement,
                               rtol=1e-12, atol=1e-8)
    # every matched random is one of the randoms
    randoms = otswap.to_cartesian(randoms_sky, table, angle_unit="deg")
    pool = {tuple(p) for p in randoms}
    assert all(tuple(p) in pool for p in result.matched_random.reshape(-1, 3))
    # displacements are of the order of the mean particle separation
    n_per_volume = n / (patch_area_deg2() / SPHERE_DEG2 * 4 / 3 * math.pi *
                        np.diff(table.distance_at(list(Z_RANGE)) ** 3)[0])
    mps = n_per_volume ** (-1 / 3)
    assert 0.2 * mps < np.median(np.linalg.norm(result.displacement, axis=2)) < 3 * mps


def test_reconstruct_lightcone_deg_and_rad_are_identical(table, lightcone):
    tracers_sky, randoms_sky = lightcone

    def rad(sky):
        return np.column_stack([np.deg2rad(sky[:, 0]), np.deg2rad(sky[:, 1]), sky[:, 2]])

    in_deg = lightcone_result(table, lightcone)
    in_rad = lightcone_result(table, (rad(tracers_sky), rad(randoms_sky)), angle_unit="rad")
    for name, array in result_arrays(in_deg).items():
        assert same_bytes(array, getattr(in_rad, name)), name


def test_reconstruct_lightcone_with_cartesian_coordinates(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    tracers = otswap.to_cartesian(tracers_sky, table, angle_unit="deg")
    randoms = otswap.to_cartesian(randoms_sky, table, angle_unit="deg")
    converted = lightcone_result(table, lightcone)
    given = lightcone_result(table, lightcone, tracers=tracers, randoms=randoms)
    for name, array in result_arrays(converted).items():
        assert same_bytes(array, getattr(given, name)), name
    # the Cartesian coordinates are used as given
    shifted = lightcone_result(table, lightcone, tracers=tracers + 1.0, randoms=randoms + 1.0)
    assert not same_bytes(shifted.displacement, converted.displacement)


def test_lightcone_cell_size_affects_speed_only(table, lightcone):
    results = [lightcone_result(table, lightcone, cell_size=cell_size)
               for cell_size in (1.0, 2.0, 4.0, 8.0)]
    for result in results[1:]:
        for name, array in result_arrays(result).items():
            assert same_bytes(array, getattr(results[0], name)), name


def test_reconstruct_lightcone_needs_keywords(table, lightcone):
    with pytest.raises(TypeError):
        otswap.reconstruct_lightcone(*lightcone, patch_area_deg2(), 3, table, "deg")
    with pytest.raises(TypeError):
        otswap.reconstruct_lightcone(*lightcone, sky_area_deg2=patch_area_deg2(), n_bins=3,
                                     distances=table)


@pytest.mark.parametrize("kwargs", [
    dict(angle_unit="rad "),
    dict(angle_unit=0),
    dict(n_bins=0),
    dict(n_bins=-1),
    dict(n_bins=100),
    dict(sky_area_deg2=-1.0),
    dict(sky_area_deg2="big"),
    dict(distances=None),
    dict(n_realizations=5),
    dict(seed="1"),
    dict(tracers=np.zeros((1600, 3))),
    dict(randoms=np.zeros((4800, 3))),
    dict(tracers=np.zeros((1600, 3)), randoms=np.zeros((10, 3))),
    dict(tracers=np.zeros((1600, 2)), randoms=np.zeros((4800, 3))),
])
def test_reconstruct_lightcone_errors(table, lightcone, kwargs):
    with pytest.raises(otswap.Error):
        lightcone_result(table, lightcone, **kwargs)


def test_reconstruct_lightcone_input_errors(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    with pytest.raises(otswap.Error):
        lightcone_result(table, (tracers_sky[:, :2], randoms_sky))
    with pytest.raises(otswap.Error):
        lightcone_result(table, (tracers_sky, randoms_sky.ravel()))
    beyond = tracers_sky.copy()
    beyond[0, 2] = 2.0
    with pytest.raises(otswap.Error):
        lightcone_result(table, (beyond, randoms_sky))


# --------------------------------------------------------------------------
# reject_mask_crossings


def test_reject_mask_crossings_updates_in_place(table, lightcone, sieve_mask):
    mask = otswap.Mask(sieve_mask)
    result = lightcone_result(table, lightcone)
    displacement = result.displacement.copy()
    matched = result.matched_random.copy()
    valid = result.valid
    valid_realizations = result.valid_realizations
    mean = result.mean_displacement
    assert valid.all()

    assert otswap.reject_mask_crossings(result, mask) is None

    assert result.filtered_nside == 256
    assert 0 < (~valid).sum() < valid.size
    # the views taken before see the change, and are the same memory
    assert same_bytes(valid, result.valid)
    assert np.shares_memory(valid, result.valid)
    assert same_bytes(valid_realizations, result.valid_realizations)
    assert np.shares_memory(mean, result.mean_displacement)
    assert np.all(valid_realizations == valid.sum(axis=0))
    # the displacements themselves are not touched
    assert same_bytes(result.displacement, displacement)
    assert same_bytes(result.matched_random, matched)
    # the mean is over the valid realizations only, NaN with none
    none_valid = valid_realizations == 0
    assert np.all(np.isnan(mean[none_valid]))
    some = ~none_valid
    weights = valid[:, some, None]
    expected = (displacement[:, some] * weights).sum(axis=0) / valid_realizations[some, None]
    np.testing.assert_allclose(mean[some], expected, rtol=1e-12, atol=1e-12)


def test_reject_mask_crossings_threshold(table, lightcone, sieve_mask):
    mask = otswap.Mask(sieve_mask)
    strict = lightcone_result(table, lightcone)
    otswap.reject_mask_crossings(strict, mask)
    loose = lightcone_result(table, lightcone)
    otswap.reject_mask_crossings(loose, mask, max_forbidden_pixels=2)
    assert np.all(loose.valid >= strict.valid)
    assert loose.valid.sum() > strict.valid.sum()
    # it only ever marks displacements invalid
    before = strict.valid.copy()
    otswap.reject_mask_crossings(strict, mask, max_forbidden_pixels=1000)
    assert same_bytes(strict.valid, before)
    otswap.reject_mask_crossings(loose, mask)
    assert same_bytes(loose.valid, before)


def test_reject_mask_crossings_errors(table, lightcone, sieve_mask, full_mask32):
    result = lightcone_result(table, lightcone)
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask))
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings(result, otswap.Mask(full_mask32))
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask), max_forbidden_pixels=-1)
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings(result, sieve_mask)
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings("result", otswap.Mask(sieve_mask))


def test_reject_mask_crossings_with_a_mask_crossed_by_none(table, lightcone, full_mask32):
    result = lightcone_result(table, lightcone)
    otswap.reject_mask_crossings(result, otswap.Mask(full_mask32))
    assert result.valid.all() and result.filtered_nside == 32


# --------------------------------------------------------------------------
# the GIL


def largest_pause_while(call):
    """Runs call in a thread; returns how long the call took, and the longest
    interval during it in which this thread could not run."""
    elapsed = []

    def run():
        t0 = time.perf_counter()
        call()
        elapsed.append(time.perf_counter() - t0)

    worker = threading.Thread(target=run)
    last = time.perf_counter()
    largest = 0.0
    worker.start()
    while worker.is_alive():
        now = time.perf_counter()
        largest = max(largest, now - last)
        last = now
    worker.join()
    return elapsed[0], largest


def test_reconstructions_release_the_gil(table, lightcone, sieve_mask):
    tracers = box_catalogue(n_side=24, spacing=10.0)
    box = lambda: otswap.reconstruct_box(tracers, mps=10.0, n_realizations=4, seed=1)
    cone = lambda: lightcone_result(table, lightcone, n_realizations=3)
    big = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=4, seed=1)
    mask = otswap.Mask(sieve_mask)
    crossings = lambda: otswap.reject_mask_crossings(big, mask, max_forbidden_pixels=1000)
    for call in (box, cone, crossings):
        elapsed, pause = largest_pause_while(call)
        assert elapsed > 0.02, "the call is too short to tell"
        assert pause < 0.5 * elapsed


# --------------------------------------------------------------------------
# bit identity with the C++ library

CPP_REFERENCE = os.environ.get("OTSWAP_CPP_REFERENCE")

needs_cpp = pytest.mark.skipif(
    not CPP_REFERENCE or not os.access(CPP_REFERENCE, os.X_OK),
    reason="OTSWAP_CPP_REFERENCE does not name the cpp_reference executable")


def run_cpp(tmp_path, mode, arrays, **params):
    args = [CPP_REFERENCE, mode]
    for name, array in arrays.items():
        path = tmp_path / f"{name}.in"
        np.ascontiguousarray(array, dtype=np.float64).tofile(path)
        args.append(f"{name}={path}")
    args += [f"{key}={float(value)!r}" if isinstance(value, float) else f"{key}={value}"
             for key, value in params.items()]
    out = tmp_path / "out"
    args.append(f"out={out}")
    subprocess.run(args, check=True, capture_output=True, text=True)
    return out


def read_cpp_result(out, n_realizations, n_objects):
    def load(suffix, dtype, shape):
        return np.fromfile(f"{out}.{suffix}", dtype=dtype).reshape(shape)
    return {
        "displacement": load("displacement", np.float64, (n_realizations, n_objects, 3)),
        "matched_random": load("matched_random", np.float64, (n_realizations, n_objects, 3)),
        "valid": load("valid", np.uint8, (n_realizations, n_objects)).astype(np.bool_),
        "valid_realizations": load("valid_realizations", np.uint32, (n_objects,)),
        "mean_displacement": load("mean_displacement", np.float64, (n_objects, 3)),
    }


def assert_same_result(result, reference):
    for name, array in result_arrays(result).items():
        assert same_bytes(array, reference[name]), name


def to_radians(sky):
    return np.column_stack([np.deg2rad(sky[:, 0]), np.deg2rad(sky[:, 1]), sky[:, 2]])


@needs_cpp
@pytest.mark.parametrize("with_randoms", [False, True])
def test_box_matches_cpp(tmp_path, with_randoms):
    tracers = box_catalogue()
    arrays = {"tracers": tracers}
    randoms = None
    if with_randoms:
        randoms = np.random.default_rng(3).uniform(0, 80, (3 * len(tracers), 3))
        arrays["randoms"] = randoms
    params = dict(mps=10.0, n_realizations=3, convergence=1e-3, seed=21, cell_size=4.0)
    out = run_cpp(tmp_path, "box", arrays, **params)
    result = otswap.reconstruct_box(tracers, randoms, **params)
    assert_same_result(result, read_cpp_result(out, 3, len(tracers)))


@needs_cpp
@pytest.mark.parametrize("cartesian", [False, True])
def test_lightcone_matches_cpp(tmp_path, table, lightcone, sieve_mask, cartesian):
    tracers_sky, randoms_sky = lightcone
    arrays = {"tracers_sky": to_radians(tracers_sky), "randoms_sky": to_radians(randoms_sky)}
    extra = {}
    if cartesian:
        extra = dict(tracers=otswap.to_cartesian(tracers_sky, table, angle_unit="deg") * 1.001,
                     randoms=otswap.to_cartesian(randoms_sky, table, angle_unit="deg") * 1.001)
        arrays.update(extra)
    config = dict(n_realizations=2, convergence=1e-3, seed=5, cell_size=4.0)
    out = run_cpp(tmp_path, "lightcone", arrays, sky_area_deg2=patch_area_deg2(), n_bins=3,
                  **TABLE_ARGS, **config, mask=sieve_mask, max_forbidden_pixels=1)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, sky_area_deg2=patch_area_deg2(),
                                          n_bins=3, distances=table, angle_unit="deg",
                                          **extra, **config)
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask), max_forbidden_pixels=1)
    assert not result.valid.all()
    assert_same_result(result, read_cpp_result(out, 2, len(tracers_sky)))


@needs_cpp
def test_to_cartesian_matches_cpp(tmp_path, table, lightcone):
    sky = lightcone[1]
    out = run_cpp(tmp_path, "cartesian", {"sky": to_radians(sky)}, **TABLE_ARGS)
    reference = np.fromfile(f"{out}.xyz", dtype=np.float64).reshape(-1, 3)
    assert same_bytes(otswap.to_cartesian(sky, table, angle_unit="deg"), reference)
