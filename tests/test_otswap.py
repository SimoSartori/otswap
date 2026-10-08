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
        "AngleUnit", "BiasTable", "CorrectionStatus", "DistanceTable", "Error", "ExtrapolationWarning", "Mask",
        "MpsProfile", "NeighbourAverage", "RealSpaceCatalog", "Result", "SelectionCounts",
        "axis_projection", "neighbour_average", "radial_projection",
        "real_space_box", "real_space_lightcone", "recompute_means", "reconstruct_box",
        "reconstruct_lightcone", "reject_mask_crossings", "rsd_factor", "rsd_factor_box",
        "shift_along_axis", "shift_radially",
        "to_cartesian", "to_sky", "Verbosity"])
    for name in otswap.__all__:
        assert hasattr(otswap, name)
    assert otswap.AngleUnit.__args__ == ("deg", "rad")
    assert otswap.Verbosity.__args__ == ("silent", "normal", "detailed")


def test_error_is_a_runtime_error():
    assert issubclass(otswap.Error, RuntimeError)
    with pytest.raises(RuntimeError):
        otswap.DistanceTable.flat(-1.0, 0.7)


def test_stub_and_marker_are_installed():
    here = os.path.dirname(otswap.__file__)
    assert os.path.isfile(os.path.join(here, "__init__.pyi"))
    assert os.path.isfile(os.path.join(here, "io.pyi"))
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


def to_degrees(sky):
    """Right ascension and declination to degrees as the binding converts
    them, with 360 folded to 0."""
    out = np.array(sky, dtype=np.float64)
    out[:, :2] *= 180.0 / math.pi
    out[out[:, 0] >= 360.0, 0] = 0.0
    return out


def test_to_sky_inverts_to_cartesian(table, lightcone):
    sky = lightcone[0]
    back = otswap.to_sky(otswap.to_cartesian(sky, table, angle_unit="deg"), table, angle_unit="deg")
    assert back.shape == sky.shape and back.dtype == np.float64
    np.testing.assert_allclose(back[:, :2], sky[:, :2], rtol=0, atol=1e-12)
    np.testing.assert_allclose(back[:, 2], sky[:, 2], rtol=0, atol=1e-12)
    rad = otswap.to_sky(otswap.to_cartesian(sky, table, angle_unit="deg"), table, angle_unit="rad")
    assert same_bytes(back, to_degrees(rad))


def test_to_sky_rows_and_the_fold(table):
    xyz = [[1000.0, -1e-17, 0.0], [np.nan, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1000.0, 1000.0],
           [-1000.0, -1e-3, 0.0]]
    deg = otswap.to_sky(xyz, table, angle_unit="deg")
    rad = otswap.to_sky(xyz, table, angle_unit="rad")
    assert deg[0, 0] == 0.0 and rad[0, 0] == 0.0
    assert np.isnan(deg[1]).all() and np.isnan(rad[1]).all()
    assert (deg[2] == 0.0).all()
    assert deg[3, 0] == pytest.approx(90.0, abs=1e-12) and deg[3, 1] == pytest.approx(45.0, abs=1e-12)
    assert (deg[[0, 2, 3, 4], 0] < 360.0).all() and (rad[[0, 2, 3, 4], 0] < 2 * math.pi).all()
    assert deg[4, 0] == pytest.approx(180.0, abs=1e-4)
    assert otswap.to_sky(np.empty((0, 3)), table, angle_unit="deg").shape == (0, 3)


@pytest.mark.parametrize("xyz, kwargs", [
    (np.zeros((4, 2)), dict(angle_unit="deg")),
    ([[1.0, np.inf, 1.0]], dict(angle_unit="deg")),
    ([[1.0e5, 0.0, 0.0]], dict(angle_unit="deg")),
    ([[1.0, 1.0, 1.0]], dict(angle_unit="degrees")),
])
def test_to_sky_errors(table, xyz, kwargs):
    with pytest.raises(otswap.Error):
        otswap.to_sky(xyz, table, **kwargs)


def test_to_sky_needs_a_table_and_a_keyword_unit(table):
    with pytest.raises(otswap.Error):
        otswap.to_sky([[1.0, 1.0, 1.0]], "table", angle_unit="deg")
    with pytest.raises(TypeError):
        otswap.to_sky([[1.0, 1.0, 1.0]], table, "deg")
    with pytest.raises(otswap.Error, match="wider redshift range"):
        otswap.to_sky([[1.0, 1.0, 1.0], [1.0e5, 0.0, 0.0]], table, angle_unit="rad")


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


# Mask.from_array

# A RING or NESTED map of NSIDE 16 holding integers from -4 to 6.
TWIN_VALUES = (np.arange(12 * 16 ** 2) * 7919) % 11 - 4


def same_mask(a, b):
    """NSIDE, sky area, and allows on a 181 x 91 grid, the poles included."""
    ra, dec = np.meshgrid(np.linspace(0, 360, 181), np.linspace(-90, 90, 91))
    return (a.nside == b.nside and a.sky_area_deg2 == b.sky_area_deg2
            and same_bytes(a.allows(ra, dec, angle_unit="deg"), b.allows(ra, dec, angle_unit="deg")))


@pytest.mark.parametrize("nest", [False, True])
@pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int32, np.uint8, np.bool_])
def test_mask_from_array_matches_its_fits_twin(tmp_path, dtype, nest):
    raw = TWIN_VALUES if np.dtype(dtype).kind in "if" else np.clip(TWIN_VALUES, 0, None)
    values = raw.astype(dtype)
    fits = otswap.Mask(write_mask(tmp_path / "twin.fits", values.astype(np.float64), 16,
                                  "NESTED" if nest else "RING"))
    assert same_mask(otswap.Mask.from_array(values, nest=nest), fits)


def test_mask_from_array_follows_the_ordering():
    ring = otswap.Mask.from_array(TWIN_VALUES)
    nested = otswap.Mask.from_array(TWIN_VALUES, nest=True)
    assert ring.nside == nested.nside == 16
    assert not same_mask(ring, nested)
    assert same_mask(ring, otswap.Mask.from_array(TWIN_VALUES, False))


def test_mask_from_array_observed_above_zero():
    values = [0.5, 1e-30, 1e-300, 1e300, np.inf, 0.0, -0.0, -1.0, -1.6375e30, np.nan, -np.inf]
    observed = [True, True, True, True, True, False, False, False, False, False, False]
    pixels = np.ones(12 * 16 ** 2)
    pixels[:len(values)] = values
    mask = otswap.Mask.from_array(pixels)
    centres = []
    for first, count in [(0, 4), (4, 8)]:
        z = 1 - (first // 4 + 1) ** 2 / (3 * 16 ** 2)
        for k in range(count):
            centres.append((360 * (k + 0.5) / count, math.degrees(math.asin(z))))
    ra, dec = np.array(centres[:len(values)]).T
    assert mask.allows(ra, dec, angle_unit="deg").tolist() == observed
    assert mask.sky_area_deg2 == pytest.approx(SPHERE_DEG2 * (pixels.size - 6) / pixels.size,
                                               rel=1e-12)


def test_mask_from_array_input_forms():
    values = TWIN_VALUES.astype(np.float64)
    reference = otswap.Mask.from_array(values.copy())

    assert same_mask(otswap.Mask.from_array(values.tolist()), reference)

    read_only = values.copy()
    read_only.flags.writeable = False
    assert same_mask(otswap.Mask.from_array(read_only), reference)

    wide = np.zeros(2 * values.size)
    wide[::2] = values
    strided = wide[::2]
    assert not strided.flags.c_contiguous
    assert same_mask(otswap.Mask.from_array(strided), reference)

    # read once: changing the source afterwards changes nothing
    source = values.copy()
    mask = otswap.Mask.from_array(source)
    source[:] = 0
    assert same_mask(mask, reference)
    mask = otswap.Mask.from_array(strided)
    wide[:] = 0
    assert same_mask(mask, reference)


def test_mask_from_array_reads_float64_in_place():
    import tracemalloc

    def traced_peak(values):
        tracemalloc.start()
        try:
            before = tracemalloc.get_traced_memory()[0]
            otswap.Mask.from_array(values)
            return tracemalloc.get_traced_memory()[1] - before
        finally:
            tracemalloc.stop()

    values = np.ones(12 * 256 ** 2)
    assert traced_peak(values) < values.nbytes // 8
    # any other dtype is converted to float64 first
    assert traced_peak(values.astype(np.float32)) >= values.nbytes


@pytest.mark.parametrize("call", [
    lambda: otswap.Mask.from_array(np.ones((12, 16))),
    lambda: otswap.Mask.from_array(np.ones(12, dtype=np.complex128)),
    lambda: otswap.Mask.from_array(["a"] * 12),
    lambda: otswap.Mask.from_array([1.0, None] * 6),
    lambda: otswap.Mask.from_array(np.float64(1.0)),
    lambda: otswap.Mask.from_array([]),
    lambda: otswap.Mask.from_array(np.ones(13)),
    lambda: otswap.Mask.from_array(np.ones(24)),
    lambda: otswap.Mask.from_array(np.ones(12 * 3 ** 2), nest=True),
    lambda: otswap.Mask.from_array(np.ones(12), nest="yes"),
    lambda: otswap.Mask.from_array(np.ones(12), nest=1),
    lambda: otswap.Mask.from_array([[1.0, 2.0], [3.0]]),
])
def test_mask_from_array_errors(call):
    with pytest.raises(otswap.Error):
        call()


def test_mask_from_array_nside_three_is_ring_only():
    assert otswap.Mask.from_array(np.ones(12 * 3 ** 2)).nside == 3


def test_mask_fits_constructor_takes_no_ordering(band_mask):
    with pytest.raises(TypeError):
        otswap.Mask(band_mask, nest=True)


def test_stub_lists_from_array():
    stub = os.path.join(os.path.dirname(otswap.__file__), "__init__.pyi")
    with open(stub) as f:
        text = f.read()
    assert "def from_array(values: ArrayLike, nest: bool = False) -> \"Mask\":" in text


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
        otswap.reconstruct_box(tracers, None, 10.0)


def test_reconstruct_box_computes_the_mps_when_not_given():
    tracers = box_catalogue()
    computed = otswap.reconstruct_box(tracers, n_realizations=2, seed=3)
    volume = np.prod(tracers.max(axis=0) - tracers.min(axis=0))
    assert computed.mps == pytest.approx((volume / len(tracers)) ** (1 / 3), rel=1e-13)
    given = otswap.reconstruct_box(tracers, n_realizations=2, seed=3, mps=computed.mps)
    assert same_bytes(computed.displacement, given.displacement) and given.mps == computed.mps
    flat = tracers.copy()
    flat[:, 2] = 1.0
    with pytest.raises(otswap.Error, match="span no volume"):
        otswap.reconstruct_box(flat)


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
# the mask in reconstruct_lightcone


def nine_mask():
    """NSIDE 64, RING, every ninth pixel unobserved."""
    pixels = np.ones(12 * 64 ** 2)
    pixels[::9] = 0
    return otswap.Mask.from_array(pixels)


MASK_CUT = (0.35, 0.55)


def test_lightcone_mask_flags_rows_and_area(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    mask = nine_mask()
    options = dict(n_bins=1, distances=table, angle_unit="deg", n_realizations=2, seed=5,
                   redshift_cut=MASK_CUT, verbosity="silent")
    plain = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask,
                                         reject_crossings=False, **options)

    def outside(sky):
        cut = ~((sky[:, 2] >= MASK_CUT[0]) & (sky[:, 2] <= MASK_CUT[1]))
        return cut, ~mask.allows(sky[:, 0], sky[:, 1], angle_unit="deg")

    cut, masked = outside(tracers_sky)
    flags = plain.outside_mask
    assert flags.dtype == np.bool_ and flags.shape == (1600,) and not flags.flags.writeable
    assert np.array_equal(flags, masked) and np.array_equal(plain.outside_redshift_cut, cut)
    assert (cut & masked).any() and (masked & ~cut).any() and (cut & ~masked).any()

    keep = ~cut & ~masked
    random_cut, random_masked = outside(randoms_sky)
    kept = otswap.reconstruct_lightcone(tracers_sky[keep], randoms_sky[~random_cut & ~random_masked],
                                        sky_area_deg2=mask.sky_area_deg2, n_bins=1, distances=table,
                                        angle_unit="deg", n_realizations=2, seed=5)
    for name, array in result_arrays(kept).items():
        full = getattr(plain, name)
        rows = full[:, keep] if full.ndim > 1 and full.shape[1] == 1600 else full[keep]
        assert same_bytes(rows, array), name
    assert np.isnan(plain.displacement[:, ~keep]).all() and not plain.valid[:, ~keep].any()


def test_lightcone_mask_filters_the_result(table, lightcone):
    mask = nine_mask()
    options = dict(mask=mask, n_bins=1, distances=table, angle_unit="deg", n_realizations=2,
                   seed=5, verbosity="silent")
    filtered = otswap.reconstruct_lightcone(*lightcone, **options)
    by_hand = otswap.reconstruct_lightcone(*lightcone, reject_crossings=False, **options)
    assert by_hand.filtered_nside == 0 and filtered.filtered_nside == 64
    otswap.reject_mask_crossings(by_hand, mask)
    assert np.array_equal(filtered.valid, by_hand.valid)
    assert 0 < (~filtered.valid[:, ~filtered.outside_mask]).sum()

    two = otswap.reconstruct_lightcone(*lightcone, max_unobserved_pixels_crossed=2, **options)
    by_hand_two = otswap.reconstruct_lightcone(*lightcone, reject_crossings=False, **options)
    otswap.reject_mask_crossings(by_hand_two, mask, max_unobserved_pixels_crossed=2)
    assert np.array_equal(two.valid, by_hand_two.valid)
    assert two.selection.max_unobserved_pixels_crossed == 2


def test_lagrangian_sky(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    mask = nine_mask()
    options = dict(mask=mask, n_bins=1, distances=table, n_realizations=2, seed=5,
                   redshift_cut=MASK_CUT, verbosity="silent")
    deg = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, angle_unit="deg", **options)
    rad = otswap.reconstruct_lightcone(to_radians(tracers_sky), to_radians(randoms_sky),
                                       angle_unit="rad", **options)
    sky = deg.lagrangian_sky
    assert sky.shape == (1600, 3) and sky.dtype == np.float64 and not sky.flags.writeable
    assert same_bytes(sky, to_degrees(rad.lagrangian_sky))

    position = otswap.to_cartesian(tracers_sky, table, angle_unit="deg") + deg.mean_displacement
    assert same_bytes(sky, otswap.to_sky(position, table, angle_unit="deg"))
    none = deg.outside_mask | deg.outside_redshift_cut | (deg.valid_realizations == 0)
    assert np.isnan(sky[none]).all() and np.isfinite(sky[~none]).all()
    assert (deg.valid_realizations[~deg.outside_mask & ~deg.outside_redshift_cut] == 0).any()
    assert ((sky[~none, 0] >= 0) & (sky[~none, 0] < 360)).all()

    otswap.reject_mask_crossings(deg, mask)
    assert same_bytes(deg.lagrangian_sky, sky), "a filter that rejects nothing keeps it"
    unfiltered = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, angle_unit="deg",
                                              reject_crossings=False, **options)
    before = unfiltered.lagrangian_sky
    assert not np.array_equal(before, sky, equal_nan=True)
    otswap.reject_mask_crossings(unfiltered, mask)
    after = unfiltered.lagrangian_sky
    assert same_bytes(after, sky), "one that rejects recomputes it, as the filter of the call"
    assert same_bytes(after, otswap.to_sky(unfiltered.lagrangian, table, angle_unit="deg"))
    assert not np.array_equal(before, after, equal_nan=True), "the earlier copy is not updated"
    box = otswap.reconstruct_box(box_catalogue(), mps=10.0)
    assert box.lagrangian_sky is None


def test_lightcone_selection_counts(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    mask = nine_mask()
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=1,
                                          distances=table, angle_unit="deg", n_realizations=2,
                                          seed=5, redshift_cut=MASK_CUT, verbosity="silent")
    s = result.selection
    assert isinstance(s, otswap.SelectionCounts)
    cut = result.outside_redshift_cut
    masked = result.outside_mask
    assert s.redshift_cut == MASK_CUT and s.mask_applied
    assert (s.tracers, s.tracers_outside_redshift_cut, s.tracers_outside_mask,
            s.tracers_outside_both) == (1600, cut.sum(), masked.sum(), (cut & masked).sum())
    random_masked = ~mask.allows(randoms_sky[:, 0], randoms_sky[:, 1], angle_unit="deg")
    assert s.randoms == 4800 and s.randoms_outside_mask == random_masked.sum()
    kept = 1600 - cut.sum() - masked.sum() + (cut & masked).sum()
    assert s.displacements == 2 * kept
    assert s.displacements_crossing_mask == 2 * kept - result.valid.sum()
    assert s.max_unobserved_pixels_crossed == 0
    lines = repr(s).split("\n")
    assert lines == [
        f"otswap: kept {kept} of 1600 tracers: {cut.sum()} outside the redshift cut [0.35, 0.55], "
        f"{masked.sum()} outside the mask ({(cut & masked).sum()} outside both)",
        f"otswap: kept {4800 - s.randoms_outside_redshift_cut - s.randoms_outside_mask + s.randoms_outside_both}"
        f" of 4800 randoms: {s.randoms_outside_redshift_cut} outside the redshift cut [0.35, 0.55], "
        f"{s.randoms_outside_mask} outside the mask ({s.randoms_outside_both} outside both)",
        f"otswap: rejected {s.displacements_crossing_mask} of {s.displacements} displacements "
        "crossing more than 0 unobserved pixels"]
    # the counts keep the result alive
    del result
    gc.collect()
    assert s.tracers == 1600

    plain = lightcone_result(table, lightcone)
    p = plain.selection
    assert p.redshift_cut is None and not p.mask_applied
    assert p.max_unobserved_pixels_crossed is None and p.displacements == 0
    assert p.tracers == 1600 and p.tracers_outside_mask == 0
    assert not plain.outside_mask.any() and repr(p) == "SelectionCounts(no selection)"
    box = otswap.reconstruct_box(box_catalogue(), mps=10.0, seed=3)
    assert box.selection.tracers == 0 and not box.outside_mask.any()
    with pytest.raises(TypeError):
        otswap.SelectionCounts()


def call_line(name, result):
    return (f"otswap: {name}: {result.n_realizations} realizations of {result.n_objects} tracers in "
            f"{result.elapsed_seconds:.2f} s; {int((result.valid_realizations == 0).sum())} without a "
            "valid displacement")


def test_lightcone_report_is_printed_from_python(table, lightcone, capfd):
    mask = nine_mask()
    options = dict(n_bins=1, distances=table, angle_unit="deg", n_realizations=2, seed=5)
    result = otswap.reconstruct_lightcone(*lightcone, mask=mask, redshift_cut=MASK_CUT, **options)
    out, err = capfd.readouterr()
    assert out == repr(result.selection) + "\n" + call_line("reconstructLightcone", result) + "\n"
    assert len(out.splitlines()) == 4
    assert err == "", "nothing is written to the process's stderr"

    otswap.reconstruct_lightcone(*lightcone, mask=mask, verbosity="silent", **options)
    assert capfd.readouterr() == ("", "")
    area = otswap.reconstruct_lightcone(*lightcone, sky_area_deg2=patch_area_deg2(), **options)
    assert capfd.readouterr() == (call_line("reconstructLightcone", area) + "\n", "")

    cut = otswap.reconstruct_lightcone(*lightcone, sky_area_deg2=patch_area_deg2(),
                                       redshift_cut=MASK_CUT, **options)
    out, err = capfd.readouterr()
    assert out.splitlines() == repr(cut.selection).split("\n") + [call_line("reconstructLightcone", cut)]
    assert "outside the redshift cut [0.35, 0.55]" in out and "mask" not in out and err == ""

    detailed = otswap.reconstruct_lightcone(*lightcone, sky_area_deg2=patch_area_deg2(),
                                            verbosity="detailed", **options)
    p = detailed.mps_profile
    assert capfd.readouterr() == (
        f"otswap: mps(z) from {p.mps[0]:g} to {p.mps[0]:g} Mpc/h in 1 bin over "
        f"[{p.redshift_min:g}, {p.redshift_max:g}]\n" + call_line("reconstructLightcone", detailed) + "\n", "")


def test_box_and_correction_reports(table, lightcone, capfd):
    tracers = box_catalogue()
    box = otswap.reconstruct_box(tracers, n_realizations=2, seed=1, verbosity="detailed")
    out, err = capfd.readouterr()
    assert out == (f"otswap: mean particle separation {box.mps:g} Mpc/h (from the tracers' bounding box)\n"
                   + call_line("reconstructBox", box) + "\n") and err == ""
    otswap.reconstruct_box(tracers, mps=10.0, verbosity="detailed")
    assert capfd.readouterr()[0].startswith("otswap: mean particle separation 10 Mpc/h (given)\n")

    r = lightcone_result(table, lightcone, n_bins=2, redshift_cut=(0.35, 0.55), verbosity="silent")
    capfd.readouterr()
    with pytest.warns(otswap.ExtrapolationWarning):
        c = otswap.real_space_lightcone(r, distances=table, bias=otswap.BiasTable([0.4, 0.5], [1.2, 1.6]),
                                        sigma=10.0, verbosity="detailed")
    out, err = capfd.readouterr()
    lines = out.splitlines()
    uncorrected = int((c.status >= 2).sum())
    left_out = int(r.outside_redshift_cut.sum())
    moved = int((c.status == 1).sum())
    assert moved == int((~np.isnan(c.sky[:, 2]) & (r.valid_realizations == 0)).sum())
    assert lines == [
        f"otswap: corrected {1600 - uncorrected} tracers, {moved} of them moved with the average of their "
        "neighbours (no valid realization)",
        f"otswap: left {uncorrected} uncorrected: {left_out} left out of the reconstruction, "
        f"{uncorrected - left_out} with no valid tracer within 3 sigma",
        f"otswap: realSpaceLightcone: corrected {1600 - uncorrected} of 1600 tracers in "
        f"{c.elapsed_seconds:.2f} s; {uncorrected} left uncorrected"], "no b(z) line: it is the warning"
    assert err == ""
    with pytest.warns(otswap.ExtrapolationWarning):
        otswap.real_space_lightcone(r, distances=table, bias=otswap.BiasTable([0.4, 0.5], [1.2, 1.6]),
                                    sigma=10.0, verbosity="silent")
    assert capfd.readouterr() == ("", ""), "silent prints nothing, and still warns"
    b = otswap.real_space_box(box, axis=1, redshift=0.5, distances=table, bias=1.5, sigma=8.0)
    n, u = len(tracers), int((b.status >= 2).sum())
    assert capfd.readouterr() == (f"otswap: realSpaceBox: corrected {n - u} of {n} tracers in "
                                  f"{b.elapsed_seconds:.2f} s; {u} left uncorrected\n", "")
    for bad in (True, "loud", None):
        with pytest.raises(otswap.Error, match="verbosity must be"):
            otswap.reconstruct_box(tracers, mps=10.0, verbosity=bad)


def test_lightcone_cartesian_with_mask(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    mask = nine_mask()
    options = dict(mask=mask, n_bins=1, distances=table, angle_unit="deg", n_realizations=2,
                   seed=5, redshift_cut=MASK_CUT, verbosity="silent")
    sky = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, **options)
    cartesian = otswap.reconstruct_lightcone(
        tracers_sky, randoms_sky, tracers=otswap.to_cartesian(tracers_sky, table, angle_unit="deg"),
        randoms=otswap.to_cartesian(randoms_sky, table, angle_unit="deg"), **options)
    for name, array in result_arrays(sky).items():
        assert same_bytes(array, getattr(cartesian, name)), name
    assert same_bytes(sky.outside_mask, cartesian.outside_mask)
    assert repr(sky.selection) == repr(cartesian.selection)


@pytest.mark.parametrize("kwargs, text", [
    (dict(sky_area_deg2=900.0, mask="mask"), "not both"),
    (dict(), "one is required"),
    (dict(mask="mask.fits"), "mask must be an otswap.Mask"),
    (dict(mask=None), "one is required"),
    (dict(mask=True, sky_area_deg2=None), "mask must be an otswap.Mask"),
    (dict(sky_area_deg2=900.0, reject_crossings=1), "reject_crossings must be True or False"),
    (dict(sky_area_deg2=900.0, verbosity=True), "verbosity must be"),
    (dict(sky_area_deg2=900.0, verbosity="loud"), "verbosity must be"),
    (dict(sky_area_deg2=900.0, max_unobserved_pixels_crossed=-1),
     "max_unobserved_pixels_crossed is -1"),
])
def test_lightcone_mask_argument_errors(table, lightcone, kwargs, text):
    kwargs = {k: (nine_mask() if v == "mask" else v) for k, v in kwargs.items()}
    with pytest.raises(otswap.Error, match=text):
        otswap.reconstruct_lightcone(*lightcone, n_bins=1, distances=table, angle_unit="deg",
                                     **kwargs)


def test_declinations_are_checked_in_the_unit_given(table, lightcone):
    tracers_sky, randoms_sky = lightcone
    above = np.nextafter(90.0, 91.0)
    bad = tracers_sky.copy()
    bad[7, 1] = above
    options = dict(sky_area_deg2=patch_area_deg2(), n_bins=1, distances=table)
    with pytest.raises(otswap.Error, match=r"tracers_sky holds a declination of 90\.000000 degrees "
                                           r"at object 7, outside \[-90, 90\]"):
        otswap.reconstruct_lightcone(bad, randoms_sky, angle_unit="deg", **options)
    bad_randoms = randoms_sky.copy()
    bad_randoms[3, 1] = -above
    with pytest.raises(otswap.Error, match=r"randoms_sky holds a declination of -90\.000000 degrees "
                                           r"at object 3"):
        otswap.reconstruct_lightcone(tracers_sky, bad_randoms, angle_unit="deg", **options)
    with pytest.raises(otswap.Error, match=r"randoms_sky holds a declination"):
        otswap.reconstruct_lightcone(tracers_sky, bad_randoms, angle_unit="deg", mask=nine_mask(),
                                     n_bins=1, distances=table)
    with pytest.raises(otswap.Error, match=r"^sky holds a declination of 90\.000000 degrees"):
        otswap.to_cartesian(bad[:10], table, angle_unit="deg")
    with pytest.raises(otswap.Error, match=r"the declination is 90\.000000 degrees; it must lie in "
                                           r"\[-90, 90\]"):
        nine_mask().allows([1.0], [above], angle_unit="deg")

    radians = np.column_stack([np.deg2rad(bad[:, 0]), np.deg2rad(bad[:, 1]), bad[:, 2]])
    radians[7, 1] = np.nextafter(np.pi / 2, 2.0)
    with pytest.raises(otswap.Error, match=r"the tracer sky array holds a declination of 1\.570796 "
                                           r"at object 7, outside \[-pi/2, pi/2\]"):
        otswap.reconstruct_lightcone(radians, to_radians(randoms_sky), angle_unit="rad", **options)
    with pytest.raises(otswap.Error, match=r"the sky array holds a declination of 1\.570796"):
        otswap.to_cartesian(radians[:10], table, angle_unit="rad")

    poles = [[10.0, 90.0, 0.4], [20.0, -90.0, 0.4]]
    assert np.isfinite(otswap.to_cartesian(poles, table, angle_unit="deg")).all()
    assert otswap.to_cartesian(np.deg2rad(poles) * [1, 1, 0] + [0, 0, 0.4], table,
                               angle_unit="rad").shape == (2, 3)


def test_stub_lists_to_sky_and_lagrangian_sky():
    stub = os.path.join(os.path.dirname(otswap.__file__), "__init__.pyi")
    with open(stub) as f:
        text = f.read()
    assert "def to_sky(" in text
    assert "    def lagrangian_sky(self) -> Optional[NDArray[np.float64]]:" in text


def test_stub_lists_the_mask_keywords():
    stub = os.path.join(os.path.dirname(otswap.__file__), "__init__.pyi")
    with open(stub) as f:
        text = f.read()
    for line in ("    sky_area_deg2: Optional[float] = None,", "    mask: Optional[Mask] = None,",
                 "    reject_crossings: bool = True,", "    max_unobserved_pixels_crossed: int = 0,",
                 "    verbosity: Verbosity = \"normal\",", "class SelectionCounts:",
                 "    def outside_mask(self) -> NDArray[np.bool_]:"):
        assert line in text, line
    assert "max_forbidden_pixels" not in text


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
    otswap.reject_mask_crossings(loose, mask, max_unobserved_pixels_crossed=2)
    assert np.all(loose.valid >= strict.valid)
    assert loose.valid.sum() > strict.valid.sum()
    # it only ever marks displacements invalid
    before = strict.valid.copy()
    otswap.reject_mask_crossings(strict, mask, max_unobserved_pixels_crossed=1000)
    assert same_bytes(strict.valid, before)
    otswap.reject_mask_crossings(loose, mask)
    assert same_bytes(loose.valid, before)


def test_reject_mask_crossings_errors(table, lightcone, sieve_mask, full_mask32):
    result = lightcone_result(table, lightcone)
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask))
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings(result, otswap.Mask(full_mask32))
    with pytest.raises(otswap.Error):
        otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask), max_unobserved_pixels_crossed=-1)
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
    crossings = lambda: otswap.reject_mask_crossings(big, mask, max_unobserved_pixels_crossed=1000)
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
                  **TABLE_ARGS, **config, mask=sieve_mask, max_unobserved_pixels_crossed=1)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, sky_area_deg2=patch_area_deg2(),
                                          n_bins=3, distances=table, angle_unit="deg",
                                          **extra, **config)
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask), max_unobserved_pixels_crossed=1)
    assert not result.valid.all()
    assert_same_result(result, read_cpp_result(out, 2, len(tracers_sky)))


@needs_cpp
@pytest.mark.parametrize("cartesian", [False, True])
def test_lightcone_mask_overload_matches_cpp(tmp_path, table, lightcone, sieve_mask, cartesian):
    tracers_sky, randoms_sky = lightcone
    arrays = {"tracers_sky": to_radians(tracers_sky), "randoms_sky": to_radians(randoms_sky)}
    extra = {}
    if cartesian:
        extra = dict(tracers=otswap.to_cartesian(tracers_sky, table, angle_unit="deg") * 1.001,
                     randoms=otswap.to_cartesian(randoms_sky, table, angle_unit="deg") * 1.001)
        arrays.update(extra)
    config = dict(n_realizations=2, convergence=1e-3, seed=5, cell_size=4.0)
    out = run_cpp(tmp_path, "lightcone", arrays, n_bins=3, **TABLE_ARGS, **config,
                  mask=sieve_mask, max_unobserved_pixels_crossed=1)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=otswap.Mask(sieve_mask),
                                          n_bins=3, distances=table, angle_unit="deg",
                                          max_unobserved_pixels_crossed=1, verbosity="silent",
                                          **extra, **config)
    assert result.outside_mask.any() and not result.valid.all()
    assert_same_result(result, read_cpp_result(out, 2, len(tracers_sky)))
    assert same_bytes(result.outside_mask.view(np.uint8),
                      np.fromfile(f"{out}.outside_mask", dtype=np.uint8))
    reference = np.fromfile(f"{out}.lagrangian_sky", dtype=np.float64).reshape(-1, 3)
    assert same_bytes(result.lagrangian_sky, to_degrees(reference))


@needs_cpp
def test_to_cartesian_matches_cpp(tmp_path, table, lightcone):
    sky = lightcone[1]
    out = run_cpp(tmp_path, "cartesian", {"sky": to_radians(sky)}, **TABLE_ARGS)
    reference = np.fromfile(f"{out}.xyz", dtype=np.float64).reshape(-1, 3)
    assert same_bytes(otswap.to_cartesian(sky, table, angle_unit="deg"), reference)


@needs_cpp
def test_to_sky_matches_cpp(tmp_path, table, lightcone):
    xyz = otswap.to_cartesian(lightcone[1], table, angle_unit="deg") * 0.999
    xyz[7] = np.nan
    out = run_cpp(tmp_path, "sky", {"cartesian": xyz}, **TABLE_ARGS)
    reference = np.fromfile(f"{out}.sky", dtype=np.float64).reshape(-1, 3)
    assert same_bytes(otswap.to_sky(xyz, table, angle_unit="rad"), reference)
    assert same_bytes(otswap.to_sky(xyz, table, angle_unit="deg"), to_degrees(reference))


# --------------------------------------------------------------------------
# io

INT64_MIN, INT64_MAX = -2 ** 63, 2 ** 63 - 1


def sample_columns():
    ids = np.array([INT64_MIN, INT64_MIN + 1, -1, 0, 2 ** 53 + 1, INT64_MAX - 1, INT64_MAX],
                   dtype=np.int64)
    x = np.array([0.5, -1.25, np.nan, 3.0, 1e-300, -0.0, 7.123456789012345])
    n = np.arange(7, dtype=np.float64) - 3
    return x, n, ids, [otswap.io.Column("x", x, description="a value", unit="Mpc/h"),
                       otswap.io.Column("n", n, type="J"),
                       otswap.io.Column("id", ids, type="K")]


@pytest.mark.parametrize("name", ["table.dat", "table.fits"])
def test_io_round_trip(tmp_path, name):
    x, n, ids, columns = sample_columns()
    path = tmp_path / name
    otswap.io.write(path, columns, precision=17)
    back = otswap.io.read(path, ["x", "n"], integer_columns=["id"])
    assert back.n_rows == 7 and back.n_columns == 2 and back.n_integer_columns == 1
    assert back.values.shape == (7, 2) and back.integers.shape == (7, 1)
    assert back.integers.dtype == np.int64 and same_bytes(back.integers[:, 0], ids)
    assert same_bytes(back.values[:, 0], x) or (
        np.array_equal(back.values[:, 0], x, equal_nan=True))
    assert np.signbit(back.values[5, 0])
    assert same_bytes(back.values[:, 1], n)
    with pytest.raises(ValueError):
        back.values[0, 0] = 1.0


def test_io_ascii_layout(tmp_path):
    path = tmp_path / "t.dat"
    otswap.io.write(str(path), [otswap.io.Column("tracX", [1.0, 2.0], description="position",
                                                 unit="Mpc/h"),
                                otswap.io.Column("n", [3, 4], type="J"),
                                otswap.io.Column("id", [5, 6], type="K")],
                    keywords=[("PRODUCT", "test", "what it is"), ("SEED", "12345", "")])
    assert path.read_text() == ("## PRODUCT = test / what it is\n## SEED = 12345\n"
                                "# tracX   position [Mpc/h]\n#\n###   tracX   n   id\n"
                                "1 3 5\n2 4 6\n")
    np.testing.assert_array_equal(np.loadtxt(path), [[1, 3, 5], [2, 4, 6]])
    assert same_bytes(otswap.io.read(path, [0, "2"]).values, np.array([[1.0, 5.0], [2.0, 6.0]]))


def test_io_precision(tmp_path):
    rng = np.random.default_rng(4)
    x = rng.standard_normal(500) * 10.0 ** rng.integers(-30, 30, 500)
    path = tmp_path / "p.dat"
    otswap.io.write(path, [otswap.io.Column("x", x)], precision=17)
    assert same_bytes(otswap.io.read(path, ["x"]).values[:, 0], x)
    otswap.io.write(path, [otswap.io.Column("x", [1 / 3])])
    assert path.read_text().endswith("0.333333333\n")
    for bad in (0, 18, -1, 1.5, "9"):
        with pytest.raises(otswap.Error):
            otswap.io.write(path, [otswap.io.Column("x", [1.0])], precision=bad)


def test_io_column(tmp_path):
    c = otswap.io.Column("id", np.array([1.0, -2.0, 2.0 ** 62]), type="K")
    assert c.data.dtype == np.int64 and c.data.tolist() == [1, -2, 2 ** 62]
    assert (c.name, c.type, c.description, c.unit) == ("id", "K", "", "")
    assert otswap.io.Column("b", [True, False], type="K").data.tolist() == [1, 0]
    assert otswap.io.Column("u", np.array([2 ** 63 - 1], dtype=np.uint64), type="K").data[0] == INT64_MAX
    d = otswap.io.Column("x", [1, 2], type="D", description="d", unit="u")
    assert d.data.dtype == np.float64 and (d.description, d.unit) == ("d", "u")
    for data, kind in (([1.5], "K"), ([2.0 ** 63], "K"), (np.array([2 ** 63], dtype=np.uint64), "K"),
                       (["a"], "K"), ([np.nan], "K"), ([[1, 2]], "D"), ([[1, 2]], "K"), ("abc", "D")):
        with pytest.raises(otswap.Error):
            otswap.io.Column("bad", data, type=kind)
    for kwargs in (dict(type="X"), dict(type="DD"), dict(type=None), dict(description=1), dict(unit=None)):
        with pytest.raises(otswap.Error):
            otswap.io.Column("bad", [1.0], **kwargs)
    with pytest.raises(otswap.Error):
        otswap.io.Column(3, [1.0])
    with pytest.raises(TypeError):
        otswap.io.Column("x", [1.0], "K")


def test_io_errors(tmp_path):
    path = tmp_path / "e.dat"
    for value in (np.nan, 1.5, 2.0 ** 31):
        with pytest.raises(otswap.Error, match=r"'J' column 'n' holds .* at row 1"):
            otswap.io.write(path, [otswap.io.Column("n", [1.0, value], type="J")])
    with pytest.raises(otswap.Error, match="holds 1 entries"):
        otswap.io.write(path, [otswap.io.Column("a", [1.0, 2.0]), otswap.io.Column("b", [1.0])])
    for columns in ([], "abc", [1.0], None):
        with pytest.raises(otswap.Error):
            otswap.io.write(path, columns)
    for keywords in ([("A", "1")], [("a", "1", "")], ["A"], [("A", 1, "")], 3):
        with pytest.raises(otswap.Error):
            otswap.io.write(path, [otswap.io.Column("a", [1.0])], keywords=keywords)
    otswap.io.write(path, [otswap.io.Column("a", [1.0])])
    for kwargs in (dict(columns="a"), dict(columns=["b"]), dict(columns=[-1]), dict(columns=[True]),
                   dict(columns=["a"], delimiter=",,"), dict(columns=["a"], comment=""),
                   dict(columns=["a"], integer_columns=["a"]), dict(columns=[])):
        with pytest.raises(otswap.Error):
            otswap.io.read(path, **kwargs)
    with pytest.raises(otswap.Error):
        otswap.io.read(3, ["a"])
    with pytest.raises(otswap.Error, match="cannot open"):
        otswap.io.read(tmp_path / "missing.fits", ["a"])
    path.write_text("1 2.5\n")
    with pytest.raises(otswap.Error, match="holds 2.5 at line 1, which is not an integer"):
        otswap.io.read(path, [0], integer_columns=[1])
    path.write_text("% c\n1;2\n")
    assert otswap.io.read(path, [1], delimiter=";", comment="%").values.tolist() == [[2.0]]


def test_io_table_keeps_its_arrays(tmp_path):
    path = tmp_path / "k.fits"
    otswap.io.write(path, [otswap.io.Column("a", [1.0, 2.0]), otswap.io.Column("b", [3, 4], type="K")])
    values = otswap.io.read(path, ["a"], integer_columns=["b"]).integers
    gc.collect()
    assert values.tolist() == [[3], [4]]


def test_bias_table(tmp_path):
    b = otswap.BiasTable([0.5, 1.0, 1.5], [1.2, 1.6, 2.1])
    assert b.redshift.tolist() == [0.5, 1.0, 1.5] and b.bias.tolist() == [1.2, 1.6, 2.1]
    with pytest.raises(ValueError):
        b.bias[0] = 2.0
    for redshift, bias in (([0.5], [1.0]), ([0.5, 0.4], [1.0, 1.0]), ([0.5, 1.0], [1.0, -1.0]),
                           ([0.5, 1.0], [1.0]), ([0.5, np.nan], [1.0, 1.0]), ([[0.5, 1.0]], [1.0, 1.0])):
        with pytest.raises(otswap.Error):
            otswap.BiasTable(redshift, bias)
    path = tmp_path / "bias.dat"
    path.write_text("# z b\n0.5 1.2\n1.0 1.6\n")
    read = otswap.io.read_bias_table(path)
    assert isinstance(read, otswap.BiasTable)
    assert read.redshift.tolist() == [0.5, 1.0] and read.bias.tolist() == [1.2, 1.6]
    fits = tmp_path / "bias.fits"
    otswap.io.write(fits, [otswap.io.Column("REDSHIFT", [0.5, 1.0]), otswap.io.Column("BIAS", [1.2, 1.6])])
    assert otswap.io.read_bias_table(str(fits)).bias.tolist() == [1.2, 1.6]
    path.write_text("0.5 1.2\n0.4 1.6\n")
    with pytest.raises(otswap.Error, match="data row 2"):
        otswap.io.read_bias_table(path)


def test_io_stub_lists_the_module():
    stub = os.path.join(os.path.dirname(otswap.__file__), "io.pyi")
    with open(stub) as f:
        text = f.read()
    for line in ("class Column:", "class Table:", "def read(", "def write(", "def read_bias_table(",
                 "def write_displacements(", "def write_displacement_field(",
                 "def write_real_space_catalog(", "def write_mps_profile("):
        assert line in text, line
    assert sorted(otswap.io.__all__) == ["Column", "Table", "read", "read_bias_table", "write",
                                         "write_displacement_field", "write_displacements",
                                         "write_mps_profile", "write_real_space_catalog"]


@needs_cpp
@pytest.mark.parametrize("name", ["table.dat", "table.fits"])
@pytest.mark.parametrize("precision", [9, 17])
def test_io_write_matches_cpp(tmp_path, name, precision):
    x, n, ids, _ = sample_columns()
    np.ascontiguousarray(x).tofile(tmp_path / "x.in")
    np.ascontiguousarray(n).tofile(tmp_path / "n.in")
    np.ascontiguousarray(ids).tofile(tmp_path / "id.in")
    cpp = tmp_path / ("cpp_" + name)
    subprocess.run([CPP_REFERENCE, "table", f"x={tmp_path / 'x.in'}", f"n={tmp_path / 'n.in'}",
                    f"id={tmp_path / 'id.in'}", f"precision={precision}", f"out={cpp}"],
                   check=True, capture_output=True, text=True)
    py = tmp_path / ("py_" + name)
    otswap.io.write(py, [otswap.io.Column("x", x, description="a value", unit="Mpc/h"),
                         otswap.io.Column("n", n, type="J", description="a count"),
                         otswap.io.Column("id", ids, type="K")],
                    precision=precision,
                    keywords=[("PRODUCT", "test table", "written by cpp_reference"), ("SEED", "12345", "")])
    assert py.read_bytes() == cpp.read_bytes()


# --------------------------------------------------------------------------
# what a result and a catalogue store


def test_result_fields_box():
    tracers = box_catalogue()
    r = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2, seed=7)
    assert r.geometry == "box" and r.angle_unit is None and r.seed == 7
    assert same_bytes(r.tracers, tracers) and not r.tracers.flags.writeable
    assert same_bytes(r.lagrangian, tracers + r.mean_displacement)
    assert r.tracers_sky is None and r.lagrangian_sky is None and r.mps_profile is None
    assert r.mps == 10.0 and r.elapsed_seconds > 0
    drawn = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2)
    assert drawn.seed != 0
    again = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2, seed=drawn.seed)
    assert same_bytes(again.displacement, drawn.displacement)


@pytest.mark.parametrize("unit", ["deg", "rad"])
def test_result_fields_lightcone(table, lightcone, unit):
    tracers_sky, randoms_sky = lightcone
    given = tracers_sky if unit == "deg" else to_radians(tracers_sky)
    randoms = randoms_sky if unit == "deg" else to_radians(randoms_sky)
    r = otswap.reconstruct_lightcone(given, randoms, sky_area_deg2=patch_area_deg2(), n_bins=2,
                                     distances=table, angle_unit=unit, n_realizations=2, seed=5,
                                     redshift_cut=(0.35, 0.55), verbosity="silent")
    assert r.geometry == "lightcone" and r.angle_unit == unit and r.mps is None
    assert same_bytes(r.tracers, otswap.to_cartesian(given, table, angle_unit=unit))
    sky = r.tracers_sky
    assert not sky.flags.writeable and r.tracers_sky is sky or unit == "rad"
    if unit == "rad":
        assert same_bytes(sky, given)
    else:
        np.testing.assert_allclose(sky, given, rtol=1e-15, atol=0)
        assert same_bytes(sky, to_degrees(to_radians(given)))
    assert same_bytes(r.lagrangian, r.tracers + r.mean_displacement)
    p = r.mps_profile
    kept = given[~r.outside_redshift_cut, 2]
    assert p.redshift.shape == (2,) and p.count.dtype == np.uint64 and p.count.sum() == kept.size
    assert p.redshift_min == kept.min() and p.redshift_max == kept.max()
    assert p.representative in p.mps
    assert same_bytes(p.at(p.redshift), p.mps)
    assert p.at([[0.4, 0.5]]).shape == (1, 2)
    with pytest.raises(ValueError):
        p.mps[0] = 1.0
    assert r.elapsed_seconds > 0


def test_selection_follows_a_later_filter(table, lightcone, sieve_mask):
    result = lightcone_result(table, lightcone)
    assert result.selection.max_unobserved_pixels_crossed is None
    before = int(result.valid.sum())
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask), max_unobserved_pixels_crossed=3)
    after3 = int(result.valid.sum())
    otswap.reject_mask_crossings(result, otswap.Mask(sieve_mask))
    after0 = int(result.valid.sum())
    s = result.selection
    assert s.max_unobserved_pixels_crossed == 0 and s.displacements == before
    assert s.displacements_crossing_mask == before - after0 and after0 < after3 < before
    assert repr(s) == (f"otswap: rejected {before - after0} of {before} displacements crossing "
                       "more than 0 unobserved pixels")


def test_result_from_arrays(table, lightcone, sieve_mask):
    r = lightcone_result(table, lightcone, n_bins=2, redshift_cut=(0.35, 0.55), verbosity="silent")
    otswap.reject_mask_crossings(r, otswap.Mask(sieve_mask))
    rebuilt = otswap.Result.from_arrays(r.displacement, r.matched_random, r.valid, r.tracers,
                                        tracers_sky=r.tracers_sky, angle_unit="deg", distances=table,
                                        outside_redshift_cut=r.outside_redshift_cut, seed=r.seed)
    for name in ("mean_displacement", "valid_realizations", "lagrangian", "lagrangian_sky",
                 "displacement", "matched_random", "valid", "tracers", "outside_redshift_cut",
                 "outside_mask"):
        assert same_bytes(getattr(rebuilt, name), getattr(r, name)), name
    assert rebuilt.geometry == "lightcone" and rebuilt.seed == r.seed and rebuilt.angle_unit == "deg"
    assert rebuilt.n_realizations == 2 and rebuilt.n_objects == 1600 and rebuilt.mps_profile is None

    box = otswap.reconstruct_box(box_catalogue(), mps=10.0, n_realizations=3, seed=2)
    valid = box.valid.copy()
    valid[1, :20] = False
    valid[:, 5] = False
    partial = otswap.Result.from_arrays(box.displacement, box.matched_random, valid, box.tracers)
    assert partial.geometry == "box" and partial.lagrangian_sky is None
    expected = (box.displacement * valid[..., None]).sum(axis=0)
    some = partial.valid_realizations > 0
    np.testing.assert_allclose(partial.mean_displacement[some],
                               expected[some] / partial.valid_realizations[some, None], rtol=1e-12)
    assert np.isnan(partial.mean_displacement[5]).all() and np.isnan(partial.lagrangian[5]).all()
    assert otswap.recompute_means(partial) is None
    assert same_bytes(partial.mean_displacement, partial.mean_displacement.copy())

    d, m, v, t = box.displacement, box.matched_random, box.valid, box.tracers
    for args, kwargs in (((d[0], m, v, t), {}), ((d, m[:, :-1], v, t), {}), ((d, m, v[:, :-1], t), {}),
                         ((d, m, v.astype(float) * 2, t), {}), ((d, m, v, t[:-1]), {}),
                         ((d, m, v, t), dict(angle_unit="deg")), ((d, m, v, t), dict(distances=table)),
                         ((d, m, v, t), dict(tracers_sky=t, angle_unit="deg")),
                         ((d, m, v, t), dict(tracers_sky=t, distances=table)),
                         ((d, m, v, t), dict(outside_mask=np.zeros(3))),
                         ((d[:, :0], m[:, :0], v[:, :0], t[:0]), {})):
        with pytest.raises(otswap.Error):
            otswap.Result.from_arrays(*args, **kwargs)
    with pytest.raises(otswap.Error):
        otswap.recompute_means("result")


def test_catalogue_stores_shift_factor_and_valid_realizations(table, lightcone):
    tracers_sky = lightcone[0]
    r = lightcone_result(table, lightcone, n_bins=2, redshift_cut=(0.35, 0.55), verbosity="silent")
    with __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", otswap.ExtrapolationWarning)
        bias = otswap.BiasTable([0.3, 0.6], [1.2, 1.6])
        c = otswap.real_space_lightcone(r, distances=table, bias=bias, sigma=10.0, verbosity="silent")
        factor = otswap.rsd_factor(tracers_sky[:, 2], table, bias=bias)
    out = r.outside_redshift_cut
    assert c.geometry == "lightcone" and c.elapsed_seconds >= 0
    assert same_bytes(c.valid_realizations, r.valid_realizations)
    assert np.isnan(c.factor[out]).all() and same_bytes(c.factor[~out], factor[~out])
    moved = ~np.isnan(c.sky[:, 2])
    assert np.isnan(c.shift[~moved]).all() and np.isfinite(c.shift[moved]).all()
    np.testing.assert_allclose(table.distance_at(c.sky[moved, 2]) - table.distance_at(tracers_sky[moved, 2]),
                               c.shift[moved], atol=1e-6)
    np.testing.assert_allclose(np.linalg.norm(c.cartesian[moved], axis=1),
                               np.linalg.norm(r.tracers[moved], axis=1) + c.shift[moved], rtol=1e-13)
    tracers = box_catalogue()
    box = otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2, seed=1)
    b = otswap.real_space_box(box, axis=2, redshift=0.5, distances=table, bias=1.5, sigma=8.0,
                              verbosity="silent")
    assert b.geometry == "box" and (b.factor == otswap.rsd_factor_box(0.5, table, bias=1.5)).all()
    np.testing.assert_allclose(b.cartesian[:, 2] - tracers[:, 2], b.shift, atol=1e-12)


# --------------------------------------------------------------------------
# the writers


def test_write_displacements(tmp_path, table, lightcone, sieve_mask):
    r = lightcone_result(table, lightcone, n_bins=2, redshift_cut=(0.35, 0.55), verbosity="silent")
    path = tmp_path / "d.dat"
    otswap.io.write_displacements(path, r)
    text = path.read_text()
    assert text.startswith("## PRODUCT = displacements / what the file holds\n## OTSWAPV = 0.1.0")
    assert "## SEED = 5 /" in text and "## ZCUTMIN = 0.35 /" in text
    names = ["tracRA", "tracDec", "tracRed", "lagrRA", "lagrDec", "lagrRed", "tracX", "tracY", "tracZ",
             "lagrX", "lagrY", "lagrZ", "displX", "displY", "displZ"]
    t = otswap.io.read(path, names, integer_columns=["nValidRec", "outsideRedshiftCut", "outsideMask"])
    assert t.n_rows == r.n_objects
    np.testing.assert_allclose(t.values[:, :3], lightcone[0], rtol=1e-8)
    np.testing.assert_allclose(t.values[:, 3:6], r.lagrangian_sky, rtol=1e-8)
    np.testing.assert_allclose(t.values[:, 9:12], r.lagrangian, rtol=1e-8)
    assert np.array_equal(t.integers[:, 0], r.valid_realizations)
    assert np.array_equal(t.integers[:, 1], r.outside_redshift_cut)

    otswap.io.write_displacements(path, r, groups=["index", "displacement"])
    assert otswap.io.read(path, ["displZ"], integer_columns=["index"]).integers[:, 0].tolist() == \
        list(range(r.n_objects))
    box = otswap.reconstruct_box(box_catalogue(), n_realizations=2, seed=1, verbosity="silent")
    for groups, match in ((["tracer_sky"], "describes a lightcone"), (["tracer", "tracer"], "twice"),
                          (["nothing"], "unknown group"), ([], "empty"), ("tracer", "sequence")):
        with pytest.raises(otswap.Error, match=match):
            otswap.io.write_displacements(path, box, groups=groups)
    with pytest.raises(otswap.Error):
        otswap.io.write_displacements(path, "result")


@pytest.mark.parametrize("name", ["field.dat", "field.fits"])
def test_write_displacement_field_round_trip(tmp_path, table, lightcone, sieve_mask, name):
    r = lightcone_result(table, lightcone, n_bins=2, redshift_cut=(0.35, 0.55), verbosity="silent")
    otswap.reject_mask_crossings(r, otswap.Mask(sieve_mask))
    path = tmp_path / name
    otswap.io.write_displacement_field(path, r)
    cols = ["tracX", "tracY", "tracZ", "lagrX", "lagrY", "lagrZ", "displX", "displY", "displZ"]
    t = otswap.io.read(path, cols, integer_columns=["realization", "index", "valid",
                                                    "outsideRedshiftCut", "outsideMask"])
    nrec, n = r.n_realizations, r.n_objects
    v = t.values.reshape(nrec, n, 9)
    k = t.integers.reshape(nrec, n, 5)
    assert (k[:, :, 0] == np.arange(nrec)[:, None]).all() and (k[:, :, 1] == np.arange(n)).all()
    back = otswap.Result.from_arrays(v[:, :, 6:9], v[:, :, 3:6], k[:, :, 2], v[0, :, 0:3],
                                     tracers_sky=r.tracers_sky, angle_unit="deg", distances=table,
                                     outside_redshift_cut=k[0, :, 3], outside_mask=k[0, :, 4], seed=r.seed)
    for field in ("displacement", "matched_random", "valid", "tracers", "mean_displacement",
                  "valid_realizations", "lagrangian", "lagrangian_sky"):
        assert same_bytes(getattr(back, field), getattr(r, field)), field
    finite = ~np.isnan(v[:, :, 6])
    assert np.array_equal(v[:, :, 6][finite], (v[:, :, 3] - v[:, :, 0])[finite])


def test_write_catalogue_and_profile(tmp_path, table, lightcone):
    r = lightcone_result(table, lightcone, n_bins=2, verbosity="silent")
    c = otswap.real_space_lightcone(r, distances=table, bias=otswap.BiasTable([0.3, 0.6], [1.2, 1.6]),
                                    sigma=10.0, verbosity="silent")
    path = tmp_path / "c.fits"
    otswap.io.write_real_space_catalog(path, c, groups=["index", "sky", "shift", "status"])
    t = otswap.io.read(path, ["tracRA", "tracDec", "tracRed", "shift", "rsdFactor"],
                       integer_columns=["index", "status"])
    np.testing.assert_allclose(t.values[:, :3], c.sky, rtol=1e-15)
    assert same_bytes(t.values[:, 3], c.shift) and same_bytes(t.values[:, 4], c.factor)
    assert np.array_equal(t.integers[:, 1], c.status)
    assert c.sigma == 10.0 and c.weight_by_realizations is False and c.axis is None
    otswap.io.write_real_space_catalog(path.with_suffix(".dat"), c)
    assert "## SIGMA = 10 /" in path.with_suffix(".dat").read_text()
    with pytest.raises(otswap.Error, match="unknown group"):
        otswap.io.write_real_space_catalog(path, c, groups=["positions"])

    profile = tmp_path / "mps.dat"
    otswap.io.write_mps_profile(profile, r.mps_profile)
    t = otswap.io.read(profile, ["redshift", "MPS"], integer_columns=["nTracers"])
    np.testing.assert_allclose(t.values[:, 1], r.mps_profile.mps, rtol=1e-8)
    assert np.array_equal(t.integers[:, 0], r.mps_profile.count)
    with pytest.raises(otswap.Error):
        otswap.io.write_mps_profile(profile, r)

    box = otswap.reconstruct_box(box_catalogue(), n_realizations=2, seed=1, verbosity="silent")
    b = otswap.real_space_box(box, axis=1, redshift=0.5, distances=table, bias=1.5, sigma=8.0,
                              verbosity="silent")
    assert (b.axis, b.box_redshift, b.box_bias) == (1, 0.5, 1.5)
    otswap.io.write_real_space_catalog(tmp_path / "b.dat", b)
    text = (tmp_path / "b.dat").read_text()
    assert "## AXIS = 1 /" in text and "###   tracX   tracY   tracZ   nValidRec   nNeighbours" in text


@needs_cpp
def test_writers_match_cpp(tmp_path, table, lightcone, sieve_mask):
    tracers_sky, randoms_sky = lightcone
    arrays = {"tracers_sky": to_radians(tracers_sky), "randoms_sky": to_radians(randoms_sky)}
    config = dict(n_realizations=2, convergence=1e-3, seed=5, cell_size=4.0)
    out = run_cpp(tmp_path, "lightcone", arrays, n_bins=3, **TABLE_ARGS, **config,
                  mask=sieve_mask, max_unobserved_pixels_crossed=1, write=1)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=otswap.Mask(sieve_mask),
                                          n_bins=3, distances=table, angle_unit="deg",
                                          max_unobserved_pixels_crossed=1, verbosity="silent", **config)
    otswap.io.write_displacements(tmp_path / "py.displacements.dat", result)
    otswap.io.write_displacement_field(tmp_path / "py.field.fits", result)
    otswap.io.write_mps_profile(tmp_path / "py.mps.dat", result.mps_profile)
    for suffix in ("displacements.dat", "field.fits", "mps.dat"):
        assert (tmp_path / f"py.{suffix}").read_bytes() == open(f"{out}.{suffix}", "rb").read(), suffix
