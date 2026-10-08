# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests for the Python bindings of the redshift-space correction.

They cover what exists at the language boundary: the shapes and dtypes of
every return value, the keyword arguments, the angle unit of the lightcone
output (that of the result corrected), the read-only views of
RealSpaceCatalog, the ExtrapolationWarning issued once per call, and the
mapping of invalid input onto otswap.Error.
The numerics are tested on the C++ side, in tests/test_rsd.cpp.
"""

import math
import warnings

import numpy as np
import pytest

import otswap

TABLE_ARGS = dict(omega_m=0.3, h=0.7, z_min=0.0, z_max=1.5, n_samples=4000)
BIAS = otswap.BiasTable([0.2, 0.8], [1.2, 1.8])
NARROW = otswap.BiasTable([0.4, 0.5], [1.2, 1.3])


@pytest.fixture(scope="module")
def table():
    return otswap.DistanceTable.flat(**TABLE_ARGS)


@pytest.fixture(scope="module")
def lightcone(table):
    """Tracers and randoms in degrees, and the reconstruction of the tracers."""
    rng = np.random.default_rng(3)

    def sky(n):
        return np.column_stack([rng.uniform(10, 40, n), rng.uniform(-15, 15, n),
                                rng.uniform(0.3, 0.6, n)])

    tracers, randoms = sky(900), sky(2700)
    result = otswap.reconstruct_lightcone(tracers, randoms, sky_area_deg2=870.0, n_bins=2,
                                          distances=table, angle_unit="deg",
                                          n_realizations=2, seed=9, verbosity="silent")
    return tracers, result


@pytest.fixture(scope="module")
def box():
    rng = np.random.default_rng(4)
    tracers = rng.uniform(0, 80, (512, 3))
    return tracers, otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2, seed=3,
                                           verbosity="silent")


def test_exports_and_warning_class():
    for name in ("ExtrapolationWarning", "RealSpaceCatalog", "BiasTable", "NeighbourAverage",
                 "radial_projection", "axis_projection", "neighbour_average", "rsd_factor",
                 "rsd_factor_box", "shift_radially", "shift_along_axis", "real_space_lightcone",
                 "real_space_box"):
        assert name in otswap.__all__ and hasattr(otswap, name)
    for gone in ("line_of_sight_projection", "shift_along_line_of_sight"):
        assert not hasattr(otswap, gone)
    assert issubclass(otswap.ExtrapolationWarning, UserWarning)
    with pytest.raises(TypeError):
        otswap.RealSpaceCatalog()


def test_projections():
    d = np.array([[1.0, 2.0, 3.0], [np.nan, 0.0, 0.0]])
    p = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 1.0]])
    radial = otswap.radial_projection(d, p)
    assert radial.shape == (2,) and radial.dtype == np.float64
    assert radial[0] == pytest.approx(11 / 5) and np.isnan(radial[1])
    assert otswap.axis_projection(d, 1).tolist()[0] == 2.0
    with pytest.raises(otswap.Error):
        otswap.axis_projection(d, 3)
    with pytest.raises(otswap.Error):
        otswap.radial_projection(d, np.zeros((2, 3)))
    with pytest.raises(otswap.Error):
        otswap.radial_projection(d, p[:1])
    with pytest.raises(TypeError):
        otswap.radial_projection(d, positions=p, axis=0)


def test_neighbour_average():
    rng = np.random.default_rng(5)
    p = rng.uniform(0, 50, (300, 3))
    v = rng.uniform(-1, 1, 300)
    r = rng.integers(0, 3, 300)
    a = otswap.neighbour_average(p, v, r, sigma=5.0)
    assert isinstance(a, otswap.NeighbourAverage) and isinstance(a, tuple)
    assert a.values.shape == (300,) and a.values.dtype == np.float64
    values, nn, nr = a
    assert values is a.values and nn is a.n_neighbours and nr is a.n_realizations_averaged
    assert nn.dtype == np.uint32 and nr.dtype == np.uint32 and nn.shape == nr.shape == (300,)
    assert np.array_equal(otswap.neighbour_average(p, v, r, sigma=0.0).values, v)
    w = otswap.neighbour_average(p, v, r, sigma=5.0, weight_by_realizations=True)
    assert not np.array_equal(w.values, a.values, equal_nan=True)
    with pytest.raises(TypeError):
        otswap.neighbour_average(p, v, r, 5.0)
    with pytest.raises(TypeError):
        otswap.neighbour_average(p, v, r, sigma=5.0, diagnostics=True)


@pytest.mark.parametrize("kwargs", [
    dict(valid_realizations=[1.5] * 4),
    dict(valid_realizations=[-1] * 4),
    dict(valid_realizations=[1] * 3),
    dict(sigma=-1.0),
    dict(sigma=math.nan),
    dict(weight_by_realizations=1),
    dict(values=[1.0, math.nan, 1.0, 1.0]),
])
def test_neighbour_average_errors(kwargs):
    args = dict(positions=np.arange(12.0).reshape(4, 3), values=[1.0] * 4,
                valid_realizations=[1] * 4, sigma=1.0)
    args.update(kwargs)
    positions, values, valid = args.pop("positions"), args.pop("values"), args.pop("valid_realizations")
    with pytest.raises(otswap.Error):
        otswap.neighbour_average(positions, values, valid, **args)


def test_rsd_factor_and_the_warning(table):
    inside = otswap.rsd_factor([0.3, 0.5], table, bias=BIAS)
    f = table.growth_rate_at([0.3, 0.5])
    b = np.interp([0.3, 0.5], BIAS.redshift, BIAS.bias)
    assert inside == pytest.approx(f / (b + 3 * f / 5), rel=1e-14)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        otswap.rsd_factor([0.3, 0.5], table, bias=BIAS)
    assert caught == []

    with pytest.warns(otswap.ExtrapolationWarning, match="2 of 3") as caught:
        otswap.rsd_factor([0.1, 0.5, 1.0], table, bias=BIAS)
    assert len(caught) == 1

    with warnings.catch_warnings():
        warnings.simplefilter("error", otswap.ExtrapolationWarning)
        with pytest.raises(otswap.ExtrapolationWarning):
            otswap.rsd_factor([1.0], table, bias=BIAS)

    with pytest.raises(otswap.Error):
        otswap.rsd_factor([0.0], table, bias=otswap.BiasTable([0.5, 1.0], [0.2, 1.2]))
    with pytest.raises(otswap.Error):
        otswap.rsd_factor([2.0], table, bias=BIAS)
    for bias in ([0.2, 0.8], None, (BIAS.redshift, BIAS.bias)):
        with pytest.raises(otswap.Error, match="BiasTable"):
            otswap.rsd_factor([0.5], table, bias=bias)
    with pytest.raises(TypeError):
        otswap.rsd_factor([0.5], table, bias_redshift=[0.2, 0.8], bias=BIAS)
    no_growth = otswap.DistanceTable.from_table([0.0, 1.0], [0.0, 3000.0])
    with pytest.raises(otswap.Error):
        otswap.rsd_factor([0.5], no_growth, bias=BIAS)


def test_rsd_factor_box(table):
    f = float(table.growth_rate_at(0.5))
    assert otswap.rsd_factor_box(0.5, table, bias=2.0) == pytest.approx(f / (2.0 + 0.6 * f), rel=1e-14)
    for bias in (0.0, -1.0, math.nan):
        with pytest.raises(otswap.Error):
            otswap.rsd_factor_box(0.5, table, bias=bias)


def test_shifts():
    p = np.array([[3.0, 4.0, 12.0], [1.0, 1.0, 1.0]])
    radial = otswap.shift_radially(p, [2.0, np.nan])
    assert radial.shape == (2, 3)
    assert np.linalg.norm(radial[0]) == pytest.approx(15.0, rel=1e-15)
    assert np.isnan(radial[1]).all()
    axis = otswap.shift_along_axis(p, [2.0, -1.0], 2)
    assert axis.tolist() == [[3.0, 4.0, 14.0], [1.0, 1.0, 0.0]]
    with pytest.raises(otswap.Error):
        otswap.shift_radially(p, [np.inf, 0.0])
    with pytest.raises(otswap.Error):
        otswap.shift_along_axis(p, [1.0, 1.0], 3)


def test_real_space_lightcone(table, lightcone):
    tracers, result = lightcone
    c = otswap.real_space_lightcone(result, distances=table, bias=BIAS, sigma=10.0,
                                    verbosity="silent")
    assert c.n_objects == len(tracers) and c.geometry == "lightcone" and c.angle_unit == "deg"
    assert c.sky.shape == (len(tracers), 3) and not c.sky.flags.writeable
    assert c.sky is c.sky, "converted once and kept"
    ok = ~np.isnan(c.sky[:, 2])
    np.testing.assert_allclose(c.sky[ok, :2], tracers[ok, :2], rtol=1e-15, atol=0,
                               err_msg="angles as given, within the rounding of the unit conversion")
    assert ((c.sky[ok, 0] >= 0) & (c.sky[ok, 0] < 360)).all()
    assert c.status.dtype == np.uint8 and set(np.unique(c.status)) <= {0, 1, 2}
    assert np.array_equal(np.flatnonzero(~ok), np.flatnonzero(c.status >= 2))
    assert np.isnan(c.cartesian[~ok]).all() and np.isfinite(c.cartesian[ok]).all()
    np.testing.assert_allclose(c.cartesian[ok], otswap.to_cartesian(c.sky[ok], table, angle_unit="deg"),
                               rtol=1e-9)
    assert c.n_neighbours.dtype == np.uint32 and c.n_realizations_averaged.dtype == np.uint32
    assert c.n_extrapolated == 0

    rad = tracers.copy()
    rad[:, :2] = np.deg2rad(rad[:, :2])
    rad_result = otswap.Result.from_arrays(result.displacement, result.matched_random, result.valid,
                                           result.tracers, tracers_sky=rad, angle_unit="rad",
                                           distances=table)
    cr = otswap.real_space_lightcone(rad_result, distances=table, bias=BIAS, sigma=10.0,
                                     verbosity="silent")
    assert cr.angle_unit == "rad"
    assert np.array_equal(cr.sky[:, 2], c.sky[:, 2], equal_nan=True)
    assert np.array_equal(cr.sky[ok, :2], rad[ok, :2]), "in radians, the input bits"
    assert np.array_equal(cr.cartesian, c.cartesian, equal_nan=True)


def test_real_space_lightcone_warns_once(table, lightcone):
    _, result = lightcone
    with pytest.warns(otswap.ExtrapolationWarning) as caught:
        c = otswap.real_space_lightcone(result, distances=table, bias=NARROW, sigma=0.0,
                                        verbosity="silent")
    assert len(caught) == 1 and c.n_extrapolated > 0
    assert f"{c.n_extrapolated} of 900 redshifts" in str(caught[0].message)


def test_real_space_box(table, box):
    tracers, result = box
    c = otswap.real_space_box(result, axis=2, redshift=0.5, distances=table, bias=1.5, sigma=8.0,
                              verbosity="silent")
    assert c.cartesian.shape == tracers.shape and c.sky is None and c.angle_unit is None
    assert c.geometry == "box"
    ok = ~np.isnan(c.cartesian[:, 0])
    assert np.array_equal(c.cartesian[ok, :2], tracers[ok, :2]), "only the axis coordinate moves"
    assert (c.status == 0).all() and c.n_extrapolated == 0


@pytest.mark.parametrize("kwargs", [
    dict(axis=3), dict(bias=0.0), dict(sigma=-1.0), dict(redshift=2.0),
    dict(weight_by_realizations="no"), dict(verbosity="loud"),
])
def test_real_space_box_errors(table, box, kwargs):
    _, result = box
    args = dict(axis=2, redshift=0.5, distances=table, bias=1.5, sigma=8.0)
    args.update(kwargs)
    with pytest.raises(otswap.Error):
        otswap.real_space_box(result, **args)


def test_real_space_lightcone_errors(table, lightcone, box):
    tracers, result = lightcone
    base = dict(distances=table, bias=BIAS, sigma=10.0, verbosity="silent")
    with pytest.raises(otswap.Error):
        otswap.real_space_lightcone("result", **base)
    with pytest.raises(otswap.Error, match="corrects a lightcone result"):
        otswap.real_space_lightcone(box[1], **base)
    with pytest.raises(otswap.Error, match="corrects a box result"):
        otswap.real_space_box(result, axis=2, redshift=0.5, distances=table, bias=1.5, sigma=8.0)
    with pytest.raises(otswap.Error, match="BiasTable"):
        otswap.real_space_lightcone(result, **dict(base, bias=[1.2, 1.8]))
    with pytest.raises(TypeError):
        otswap.real_space_lightcone(result, tracers, **base)
    # Every tracer placed at the table's last redshift: about half the shifts
    # point outwards, past the table.
    short = otswap.DistanceTable.flat(0.3, 0.7, z_max=0.6, n_samples=2000)
    at_edge = tracers.copy()
    at_edge[:, 2] = 0.6
    edge = otswap.Result.from_arrays(result.displacement, result.matched_random, result.valid,
                                     otswap.to_cartesian(at_edge, short, angle_unit="deg"),
                                     tracers_sky=at_edge, angle_unit="deg", distances=short)
    with pytest.raises(otswap.Error, match="wider redshift range"):
        otswap.real_space_lightcone(edge, **dict(base, distances=short, sigma=0.0))


def test_result_agrees_with_its_tracers(table, lightcone, box):
    """matched_random - displacement gives back the tracers the result
    carries, within rounding."""
    tracers, result = lightcone
    cart = otswap.to_cartesian(tracers, table, angle_unit="deg")
    assert np.array_equal(result.tracers, cart)
    back = result.matched_random - result.displacement
    assert np.allclose(back, cart[None], rtol=0, atol=1e-9)
    t, r = box
    assert np.array_equal(r.tracers, t)
    assert np.allclose(r.matched_random - r.displacement, t[None], rtol=0, atol=1e-9)


def test_redshift_cut(table):
    rng = np.random.default_rng(8)

    def sky(n):
        return np.column_stack([rng.uniform(10, 40, n), rng.uniform(-15, 15, n),
                                rng.uniform(0.2, 0.7, n)])

    tracers, randoms = sky(900), sky(3600)
    options = dict(sky_area_deg2=870.0, n_bins=1, distances=table, angle_unit="deg",
                   n_realizations=2, seed=4, verbosity="silent")
    cut = otswap.reconstruct_lightcone(tracers, randoms, redshift_cut=(0.3, 0.6), **options)
    inside = (tracers[:, 2] >= 0.3) & (tracers[:, 2] <= 0.6)
    flags = cut.outside_redshift_cut
    assert flags.dtype == np.bool_ and flags.shape == (900,) and not flags.flags.writeable
    assert np.array_equal(flags, ~inside)

    kept_randoms = randoms[(randoms[:, 2] >= 0.3) & (randoms[:, 2] <= 0.6)]
    kept = otswap.reconstruct_lightcone(tracers[inside], kept_randoms, **options)
    assert np.array_equal(cut.displacement[:, inside], kept.displacement)
    assert np.isnan(cut.displacement[:, ~inside]).all()
    assert not cut.valid[:, ~inside].any() and (cut.valid_realizations[~inside] == 0).all()

    plain = otswap.reconstruct_lightcone(tracers, randoms, **options)
    assert not plain.outside_redshift_cut.any()

    c = otswap.real_space_lightcone(cut, distances=table, bias=BIAS, sigma=10.0, verbosity="silent")
    assert (c.status[~inside] == otswap.CorrectionStatus.LEFT_OUT).all()
    assert (c.status[inside] != otswap.CorrectionStatus.LEFT_OUT).all()
    assert [int(v) for v in otswap.CorrectionStatus] == [0, 1, 2, 3]
    assert (c.n_neighbours[~inside] == 0).all() and np.isnan(c.factor[~inside]).all()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        otswap.real_space_lightcone(cut, distances=table, bias=NARROW, sigma=10.0, verbosity="silent")
    assert len(caught) == 1 and f" of {inside.sum()} redshifts" in str(caught[0].message)

    for bad in [(0.3,), (0.3, 0.5, 0.6), (math.nan, 0.6), (0.6, 0.3), "low-high"]:
        with pytest.raises(otswap.Error):
            otswap.reconstruct_lightcone(tracers, randoms, redshift_cut=bad, **options)


def test_mask(table):
    """Masked tracers are left out of the correction as cut ones are, and the
    extrapolation warning counts the tracers kept only."""
    rng = np.random.default_rng(8)

    def sky(n):
        return np.column_stack([rng.uniform(10, 40, n), rng.uniform(-15, 15, n),
                                rng.uniform(0.3, 0.6, n)])

    tracers, randoms = sky(900), sky(3600)
    pixels = np.ones(12 * 64 ** 2)
    pixels[::9] = 0
    mask = otswap.Mask.from_array(pixels)
    result = otswap.reconstruct_lightcone(tracers, randoms, mask=mask, n_bins=1, distances=table,
                                          angle_unit="deg", n_realizations=2, seed=4,
                                          verbosity="silent")
    masked = ~mask.allows(tracers[:, 0], tracers[:, 1], angle_unit="deg")
    assert np.array_equal(result.outside_mask, masked) and masked.any()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        c = otswap.real_space_lightcone(result, distances=table, bias=NARROW, sigma=10.0,
                                        verbosity="silent")
    assert (c.status[masked] == 3).all()
    assert (c.n_neighbours[masked] == 0).all() and np.isnan(c.sky[masked]).all()
    assert np.isnan(c.cartesian[masked]).all()
    assert len(caught) == 1 and f" of {(~masked).sum()} redshifts" in str(caught[0].message)
