# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests for the Python bindings of the redshift-space correction.

They cover what exists at the language boundary: the shapes and dtypes of
every return value, the keyword arguments, the angle unit of the lightcone
output, the read-only views of RealSpaceCatalog, the ExtrapolationWarning
issued once per call, and the mapping of invalid input onto otswap.Error.
The numerics are tested on the C++ side, in tests/test_rsd.cpp.
"""

import math
import warnings

import numpy as np
import pytest

import otswap

TABLE_ARGS = dict(omega_m=0.3, h=0.7, z_min=0.0, z_max=1.5, n_samples=4000)
BIAS_Z = [0.2, 0.8]
BIAS = [1.2, 1.8]


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
                                          n_realizations=2, seed=9)
    return tracers, result


@pytest.fixture(scope="module")
def box():
    rng = np.random.default_rng(4)
    tracers = rng.uniform(0, 80, (512, 3))
    return tracers, otswap.reconstruct_box(tracers, mps=10.0, n_realizations=2, seed=3)


def test_exports_and_warning_class():
    for name in ("ExtrapolationWarning", "RealSpaceCatalog", "line_of_sight_projection",
                 "neighbour_average", "rsd_factor", "rsd_factor_box", "shift_along_line_of_sight",
                 "real_space_lightcone", "real_space_box"):
        assert name in otswap.__all__ and hasattr(otswap, name)
    assert issubclass(otswap.ExtrapolationWarning, UserWarning)
    with pytest.raises(TypeError):
        otswap.RealSpaceCatalog()


def test_line_of_sight_projection():
    d = np.array([[1.0, 2.0, 3.0], [np.nan, 0.0, 0.0]])
    p = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 1.0]])
    radial = otswap.line_of_sight_projection(d, positions=p)
    assert radial.shape == (2,) and radial.dtype == np.float64
    assert radial[0] == pytest.approx(11 / 5) and np.isnan(radial[1])
    assert otswap.line_of_sight_projection(d, axis=1).tolist()[0] == 2.0
    for kwargs in ({}, dict(positions=p, axis=0)):
        with pytest.raises(otswap.Error):
            otswap.line_of_sight_projection(d, **kwargs)
    with pytest.raises(otswap.Error):
        otswap.line_of_sight_projection(d, axis=3)
    with pytest.raises(otswap.Error):
        otswap.line_of_sight_projection(d, positions=np.zeros((2, 3)))


def test_neighbour_average():
    rng = np.random.default_rng(5)
    p = rng.uniform(0, 50, (300, 3))
    v = rng.uniform(-1, 1, 300)
    r = rng.integers(0, 3, 300)
    a = otswap.neighbour_average(p, v, r, sigma=5.0)
    assert a.shape == (300,) and a.dtype == np.float64
    b, nn, nr = otswap.neighbour_average(p, v, r, sigma=5.0, diagnostics=True)
    assert np.array_equal(a, b, equal_nan=True)
    assert nn.dtype == np.uint32 and nr.dtype == np.uint32 and nn.shape == nr.shape == (300,)
    assert np.array_equal(otswap.neighbour_average(p, v, r, sigma=0.0), v)
    w = otswap.neighbour_average(p, v, r, sigma=5.0, weight_by_realizations=True)
    assert not np.array_equal(w, a, equal_nan=True)
    with pytest.raises(TypeError):
        otswap.neighbour_average(p, v, r, 5.0)


@pytest.mark.parametrize("kwargs", [
    dict(valid_realizations=[1.5] * 4),
    dict(valid_realizations=[-1] * 4),
    dict(valid_realizations=[1] * 3),
    dict(sigma=-1.0),
    dict(sigma=math.nan),
    dict(weight_by_realizations=1),
    dict(diagnostics="yes"),
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
    inside = otswap.rsd_factor([0.3, 0.5], table, bias_redshift=BIAS_Z, bias=BIAS)
    f = table.growth_rate_at([0.3, 0.5])
    b = np.interp([0.3, 0.5], BIAS_Z, BIAS)
    assert inside == pytest.approx(f / (b + 3 * f / 5), rel=1e-14)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        otswap.rsd_factor([0.3, 0.5], table, bias_redshift=BIAS_Z, bias=BIAS)
    assert caught == []

    with pytest.warns(otswap.ExtrapolationWarning, match="2 of 3") as caught:
        otswap.rsd_factor([0.1, 0.5, 1.0], table, bias_redshift=BIAS_Z, bias=BIAS)
    assert len(caught) == 1

    with warnings.catch_warnings():
        warnings.simplefilter("error", otswap.ExtrapolationWarning)
        with pytest.raises(otswap.ExtrapolationWarning):
            otswap.rsd_factor([1.0], table, bias_redshift=BIAS_Z, bias=BIAS)

    with pytest.raises(otswap.Error):
        otswap.rsd_factor([0.0], table, bias_redshift=[0.5, 1.0], bias=[0.2, 1.2])
    with pytest.raises(otswap.Error):
        otswap.rsd_factor([0.5], table, bias_redshift=[0.5], bias=[1.0])
    with pytest.raises(otswap.Error):
        otswap.rsd_factor([2.0], table, bias_redshift=BIAS_Z, bias=BIAS)
    no_growth = otswap.DistanceTable.from_table([0.0, 1.0], [0.0, 3000.0])
    with pytest.raises(otswap.Error):
        otswap.rsd_factor([0.5], no_growth, bias_redshift=BIAS_Z, bias=BIAS)


def test_rsd_factor_box(table):
    f = float(table.growth_rate_at(0.5))
    assert otswap.rsd_factor_box(0.5, table, bias=2.0) == pytest.approx(f / (2.0 + 0.6 * f), rel=1e-14)
    for bias in (0.0, -1.0, math.nan):
        with pytest.raises(otswap.Error):
            otswap.rsd_factor_box(0.5, table, bias=bias)


def test_shift_along_line_of_sight():
    p = np.array([[3.0, 4.0, 12.0], [1.0, 1.0, 1.0]])
    radial = otswap.shift_along_line_of_sight(p, [2.0, np.nan])
    assert radial.shape == (2, 3)
    assert np.linalg.norm(radial[0]) == pytest.approx(15.0, rel=1e-15)
    assert np.isnan(radial[1]).all()
    axis = otswap.shift_along_line_of_sight(p, [2.0, -1.0], axis=2)
    assert axis.tolist() == [[3.0, 4.0, 14.0], [1.0, 1.0, 0.0]]
    with pytest.raises(otswap.Error):
        otswap.shift_along_line_of_sight(p, [np.inf, 0.0])


def test_real_space_lightcone(table, lightcone):
    tracers, result = lightcone
    c = otswap.real_space_lightcone(result, tracers, distances=table, bias_redshift=BIAS_Z,
                                    bias=BIAS, sigma=10.0, angle_unit="deg")
    assert c.n_objects == len(tracers)
    assert c.positions.shape == (len(tracers), 3)
    assert not c.positions.flags.writeable
    ok = ~np.isnan(c.positions[:, 2])
    assert np.array_equal(c.positions[ok, :2], tracers[ok, :2]), "angles come back as given, in degrees"
    assert np.array_equal(np.flatnonzero(~ok), c.uncorrected)
    assert c.uncorrected.dtype == np.int64
    assert c.n_neighbours.dtype == np.uint32 and c.n_realizations_averaged.dtype == np.uint32

    rad = tracers.copy()
    rad[:, :2] = np.deg2rad(rad[:, :2])
    cr = otswap.real_space_lightcone(result, rad, distances=table, bias_redshift=BIAS_Z,
                                     bias=BIAS, sigma=10.0, angle_unit="rad")
    assert np.array_equal(cr.positions[:, 2], c.positions[:, 2], equal_nan=True)
    assert np.array_equal(cr.positions[ok, :2], rad[ok, :2])


def test_real_space_lightcone_warns_once(table, lightcone):
    tracers, result = lightcone
    with pytest.warns(otswap.ExtrapolationWarning) as caught:
        otswap.real_space_lightcone(result, tracers, distances=table, bias_redshift=[0.4, 0.5],
                                    bias=[1.4, 1.5], sigma=0.0, angle_unit="deg")
    assert len(caught) == 1


def test_real_space_box(table, box):
    tracers, result = box
    c = otswap.real_space_box(result, tracers, axis=2, redshift=0.5, distances=table, bias=1.5,
                              sigma=8.0)
    assert c.positions.shape == tracers.shape
    ok = ~np.isnan(c.positions[:, 0])
    assert np.array_equal(c.positions[ok, :2], tracers[ok, :2]), "only the axis coordinate moves"
    assert len(c.uncorrected) == 0
    assert c.uncorrected.dtype == np.int64


@pytest.mark.parametrize("kwargs", [
    dict(axis=3), dict(bias=0.0), dict(sigma=-1.0), dict(redshift=2.0),
    dict(weight_by_realizations="no"), dict(tracers=np.zeros((5, 3))),
])
def test_real_space_box_errors(table, box, kwargs):
    tracers, result = box
    args = dict(tracers=tracers, axis=2, redshift=0.5, distances=table, bias=1.5, sigma=8.0)
    args.update(kwargs)
    t = args.pop("tracers")
    with pytest.raises(otswap.Error):
        otswap.real_space_box(result, t, **args)


def test_real_space_lightcone_errors(table, lightcone):
    tracers, result = lightcone
    base = dict(distances=table, bias_redshift=BIAS_Z, bias=BIAS, sigma=10.0, angle_unit="deg")
    with pytest.raises(otswap.Error):
        otswap.real_space_lightcone(result, tracers[:-1], **base)
    with pytest.raises(otswap.Error):
        otswap.real_space_lightcone(result, tracers, **dict(base, angle_unit="radians"))
    with pytest.raises(otswap.Error):
        otswap.real_space_lightcone("result", tracers, **base)
    # Every tracer placed at the table's last redshift: about half the shifts
    # point outwards, past the table. Only the count is checked against the
    # result, so the moved tracers are accepted.
    short = otswap.DistanceTable.flat(0.3, 0.7, z_max=0.6, n_samples=2000)
    at_edge = tracers.copy()
    at_edge[:, 2] = 0.6
    with pytest.raises(otswap.Error, match="wider redshift range"):
        otswap.real_space_lightcone(result, at_edge, **dict(base, distances=short, sigma=0.0))


def test_result_agrees_with_its_tracers(table, lightcone, box):
    """matched_random - displacement gives back the tracers the result was
    reconstructed from: a check the library leaves to its callers."""
    tracers, result = lightcone
    cart = otswap.to_cartesian(tracers, table, angle_unit="deg")
    back = result.matched_random - result.displacement
    assert np.allclose(back, cart[None], rtol=0, atol=1e-9)
    t, r = box
    assert np.allclose(r.matched_random - r.displacement, t[None], rtol=0, atol=1e-9)


def test_redshift_cut(table):
    rng = np.random.default_rng(8)

    def sky(n):
        return np.column_stack([rng.uniform(10, 40, n), rng.uniform(-15, 15, n),
                                rng.uniform(0.2, 0.7, n)])

    tracers, randoms = sky(900), sky(3600)
    options = dict(sky_area_deg2=870.0, n_bins=1, distances=table, angle_unit="deg",
                   n_realizations=2, seed=4)
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

    c = otswap.real_space_lightcone(cut, tracers, distances=table, bias_redshift=[0.2, 0.8],
                                    bias=[1.2, 1.8], sigma=10.0, angle_unit="deg")
    assert set(np.flatnonzero(~inside)) <= set(c.uncorrected.tolist())
    assert (c.n_neighbours[~inside] == 0).all()

    for bad in [(0.3,), (0.3, 0.5, 0.6), (math.nan, 0.6), (0.6, 0.3), "low-high"]:
        with pytest.raises(otswap.Error):
            otswap.reconstruct_lightcone(tracers, randoms, redshift_cut=bad, **options)
