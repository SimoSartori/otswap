# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests for the use of otswap's types from another nanobind module.

_otswap_crosstest (tests/crossmodule/crosstest.cpp) binds no type: it takes
and returns Config, CorrectionConfig, RedshiftCut, BoxLimits, Mask,
DistanceTable, MpsProfile, BiasTable, Result and RealSpaceCatalog as their
C++ types, which nanobind finds in the registry otswap fills. These tests pass
objects made in Python to it and read in Python the objects it makes: their
types, their fields on both sides, the angle unit of a Result or catalogue
made in C++ (radians), the weak-keyed tables that hold the unit, the size
checks of the array getters on a malformed Result, and the errors.

The module is built with the tests and the Python bindings, into the
crossmodule/ directory of the build tree, and never installed. The tests
need that directory in the environment variable OTSWAP_CROSSTEST_DIR, and
are skipped without it.
"""

import gc
import math
import os
import re
import sys

import numpy as np
import pytest

import otswap

CROSSTEST_DIR = os.environ.get("OTSWAP_CROSSTEST_DIR")

if CROSSTEST_DIR:
    sys.path.insert(0, CROSSTEST_DIR)
    import _otswap_crosstest as cross
else:
    cross = None

pytestmark = pytest.mark.skipif(
    cross is None,
    reason="OTSWAP_CROSSTEST_DIR does not name the directory of the _otswap_crosstest module")

QUIET = otswap.Config(n_realizations=2, seed=3, verbosity="silent")
QUIET_CORRECTION = otswap.CorrectionConfig(verbosity="silent")
SIGMA = 10.0


def same_bytes(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def full_match(text):
    return "^" + re.escape(text) + "$"


@pytest.fixture(scope="module")
def table():
    """The table cross.make_table() builds."""
    return otswap.DistanceTable.flat(0.3, 0.7, z_min=0.0, z_max=1.5, n_samples=4000)


@pytest.fixture(scope="module")
def bias():
    return otswap.BiasTable([0.2, 0.8], [1.2, 1.8])


@pytest.fixture(scope="module")
def lattice():
    """The lattice of cross.small_box_result()."""
    axis = 10.0 * np.arange(6)
    return np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)


@pytest.fixture(scope="module")
def sky():
    rng = np.random.default_rng(5)

    def sample(n):
        return np.column_stack([rng.uniform(10, 40, n), rng.uniform(-15, 15, n),
                                rng.uniform(0.3, 0.6, n)])

    return sample(1200), sample(3600)


def lightcone_result(sky, table):
    """A lightcone result made by otswap, in degrees."""
    tracers, randoms = sky
    return otswap.reconstruct_lightcone(tracers, randoms, sky_area_deg2=870.0, n_bins=2,
                                        distances=table, angle_unit="deg", config=QUIET)


@pytest.fixture(scope="module")
def box_result(lattice):
    return otswap.reconstruct_box(lattice, mps=10.0, config=QUIET)


@pytest.fixture(scope="module")
def lc_result(sky, table):
    return lightcone_result(sky, table)


# --------------------------------------------------------------------------
# types and configuration


@pytest.mark.parametrize("make, args, cls", [
    ("default_config", (), "Config"),
    ("default_correction", (), "CorrectionConfig"),
    ("make_cut", (0.2, 0.8), "RedshiftCut"),
    ("make_box", (), "BoxLimits"),
    ("make_mask", (), "Mask"),
    ("make_table", (), "DistanceTable"),
    ("make_bias_table", (), "BiasTable"),
    ("small_box_result", (), "Result"),
    ("small_lightcone_result", (), "Result"),
    ("small_box_catalog", (), "RealSpaceCatalog"),
    ("small_lightcone_catalog", (), "RealSpaceCatalog"),
])
def test_returned_types_are_otswap_types(make, args, cls):
    assert type(getattr(cross, make)(*args)) is getattr(otswap, cls)


def test_python_defaults_are_the_cpp_ones():
    assert otswap.Config() == cross.default_config()
    assert otswap.CorrectionConfig() == cross.default_correction()
    assert cross.config_fields(otswap.Config())["log_is_clog"]


def test_config_both_ways():
    config = otswap.Config(n_realizations=3, convergence=0.01, seed=5, cell_size=2.5,
                           reject_crossings=False, max_unobserved_pixels_crossed=2,
                           verbosity="detailed")
    assert cross.config_fields(config) == dict(
        n_realizations=3, convergence=0.01, seed=5, cell_size=2.5, reject_crossings=False,
        max_unobserved_pixels_crossed=2, verbosity="detailed", log_is_clog=True)
    changed = cross.changed_config(config)
    assert changed == otswap.Config(n_realizations=4, convergence=0.02, seed=6, cell_size=3.5,
                                    reject_crossings=True, max_unobserved_pixels_crossed=3,
                                    verbosity="silent")
    assert (changed.n_realizations, changed.verbosity) == (4, "silent")
    assert config.n_realizations == 3


def test_correction_both_ways():
    correction = otswap.CorrectionConfig(weight_by_realizations=True, verbosity="silent")
    assert cross.correction_fields(correction) == dict(weight_by_realizations=True,
                                                       verbosity="silent", log_is_clog=True)
    assert cross.changed_correction(correction) == otswap.CorrectionConfig(
        weight_by_realizations=False, verbosity="normal")


def test_redshift_cut_both_ways_and_implicit_conversion():
    assert cross.cut_bounds(otswap.RedshiftCut(0.2, 0.8)) == (0.2, 0.8)
    assert cross.cut_bounds(otswap.RedshiftCut()) == (-math.inf, math.inf)
    assert cross.cut_bounds((0.2, 0.8)) == (0.2, 0.8)
    cut = cross.make_cut(0.1, 0.9)
    assert (cut.min, cut.max) == (0.1, 0.9)
    for wrong in ((0.1, 0.2, 0.3), "x", None):
        with pytest.raises(TypeError):
            cross.cut_bounds(wrong)


def test_box_limits_both_ways_and_implicit_conversion():
    corners = ((0.0, 0.0, 0.0), (1.0, 2.0, 3.0))
    assert cross.box_corners(otswap.BoxLimits(*corners)) == corners
    assert cross.box_corners(((0, 0, 0), (1, 2, 3))) == corners
    box = cross.make_box()
    assert (box.lower, box.upper) == ((0.0, 0.0, 0.0), (50.0, 60.0, 70.0))
    for wrong in (((0, 0), (1, 1)), "x", None):
        with pytest.raises(TypeError):
            cross.box_corners(wrong)


def test_tables_and_mask_both_ways(table, bias, sky):
    assert cross.mask_nside(otswap.Mask.from_array(np.ones(12 * 4 ** 2))) == 4
    mask = cross.make_mask()
    assert mask.nside == 2
    assert mask.sky_area_deg2 == pytest.approx(47 / 48 * 4 * math.pi * (180 / math.pi) ** 2)
    assert cross.distance_at(table, 0.5) == table.distance_at(0.5)
    assert same_bytes(cross.make_table().distance_at([0.1, 0.5, 1.2]),
                      table.distance_at([0.1, 0.5, 1.2]))
    assert cross.bias_nodes(bias) == ([0.2, 0.8], [1.2, 1.8])
    made = cross.make_bias_table()
    assert same_bytes(made.redshift, bias.redshift) and same_bytes(made.bias, bias.bias)
    profile = otswap.measure_mps_profile(sky[0], sky_area_deg2=870.0, n_bins=3, distances=table,
                                         angle_unit="deg", verbosity="silent")
    assert cross.profile_bins(profile) == 3


# --------------------------------------------------------------------------
# Result


def test_otswap_results_read_by_cpp(box_result, lc_result, sky):
    fields = cross.result_fields(box_result)
    assert fields == dict(n_objects=216, n_realizations=2, geometry="box", seed=3,
                          has_distances=False, first_tracers_sky=None, log_is_clog=True)

    before = lc_result.tracers_sky.copy()
    fields = cross.result_fields(lc_result)
    assert (fields["geometry"], fields["n_objects"], fields["has_distances"]) == ("lightcone", 1200, True)
    # C++ sees radians, converted from the degrees of the call
    expected = (*np.deg2rad(sky[0][0, :2]), sky[0][0, 2])
    assert fields["first_tracers_sky"] == expected
    assert lc_result.angle_unit == "deg"
    assert same_bytes(lc_result.tracers_sky, before)


def test_cpp_box_result_used_by_otswap(tmp_path, table):
    result = cross.small_box_result()
    assert result.angle_unit is None and result.geometry == "box"
    assert (result.n_objects, result.n_realizations) == (216, 2)
    assert (result.config.seed, result.config.n_realizations) == (7, 2)
    assert result.tracers_sky is None and result.lagrangian_sky is None and result.distances is None
    assert result.mean_displacement.shape == (216, 3)
    # the mps of the bounding box of the tracers (drawn randoms), [0, 50]^3
    assert result.mps == pytest.approx((50.0 ** 3 / 216) ** (1 / 3))
    catalog = otswap.real_space_box(result, axis=2, redshift=0.5, distances=table, bias=1.5,
                                    sigma=SIGMA, correction=QUIET_CORRECTION)
    assert catalog.n_objects == 216 and catalog.angle_unit is None
    otswap.io.write_displacement_mean(tmp_path / "mean.dat", result)
    assert (tmp_path / "mean.dat").stat().st_size > 0


def test_cpp_lightcone_result_used_by_otswap():
    result = cross.small_lightcone_result()
    fields = cross.result_fields(result)
    assert result.angle_unit == "rad" and result.geometry == "lightcone"
    assert result not in otswap._otswap._result_state
    sky = result.tracers_sky
    assert sky.shape == (1200, 3) and not sky.flags.writeable
    assert tuple(sky[0]) == fields["first_tracers_sky"]
    assert result.distances is not None
    otswap.reject_mask_crossings(result, cross.make_mask())
    assert result.filtered_nside == 2
    otswap.recompute_means(result)
    assert result.lagrangian_sky.shape == (1200, 3)
    assert result.angle_unit == "rad"


# --------------------------------------------------------------------------
# RealSpaceCatalog


def test_otswap_catalogues_read_by_cpp(box_result, lc_result, table, bias):
    catalog = otswap.real_space_box(box_result, axis=2, redshift=0.5, distances=table, bias=1.5,
                                    sigma=SIGMA, correction=QUIET_CORRECTION)
    fields = cross.catalog_fields(catalog)
    assert (fields["geometry"], fields["axis"], fields["sigma"], fields["first_sky"]) == ("box", 2, SIGMA, None)
    assert sum(fields["status_counts"]) == 216

    catalog = otswap.real_space_lightcone(lc_result, distances=table, bias=bias, sigma=SIGMA,
                                          correction=QUIET_CORRECTION)
    assert catalog.angle_unit == "deg"
    fields = cross.catalog_fields(catalog)
    assert fields["geometry"] == "lightcone" and sum(fields["status_counts"]) == 1200
    # C++ sees radians, the catalogue returns degrees
    first = np.array(fields["first_sky"])
    np.testing.assert_allclose(np.concatenate([np.rad2deg(first[:2]), first[2:]]), catalog.sky[0],
                               rtol=1e-14, equal_nan=True)


def test_cpp_catalogues_read_by_otswap():
    catalog = cross.small_box_catalog()
    assert catalog.geometry == "box" and catalog.angle_unit is None and catalog.sky is None
    assert (catalog.axis, catalog.sigma, catalog.n_objects) == (2, SIGMA, 216)
    assert catalog.cartesian.shape == (216, 3)

    catalog = cross.small_lightcone_catalog()
    fields = cross.catalog_fields(catalog)
    assert catalog.geometry == "lightcone" and catalog.angle_unit == "rad"
    assert catalog not in otswap._otswap._catalog_state
    assert same_bytes(catalog.sky[0], np.array(fields["first_sky"]))
    assert catalog.sky.shape == (1200, 3) and not catalog.sky.flags.writeable


def test_catalogue_of_a_cpp_result_is_in_radians(bias):
    catalog = otswap.real_space_lightcone(cross.small_lightcone_result(), distances=cross.make_table(),
                                          bias=bias, sigma=SIGMA, correction=QUIET_CORRECTION)
    assert catalog.angle_unit == "rad"
    assert same_bytes(catalog.sky, cross.small_lightcone_catalog().sky)


@pytest.mark.parametrize("extension", [".dat", ".fits"])
def test_catalogue_writer_same_bytes_both_ways(tmp_path, extension):
    made_in_cpp = tmp_path / ("cpp" + extension)
    made_in_python = tmp_path / ("python" + extension)
    otswap.io.write_real_space_catalog(made_in_cpp, cross.small_box_catalog())
    catalog = otswap.real_space_box(cross.small_box_result(), axis=2, redshift=0.5,
                                    distances=cross.make_table(), bias=1.5, sigma=SIGMA,
                                    correction=QUIET_CORRECTION)
    otswap.io.write_real_space_catalog(made_in_python, catalog)
    assert made_in_cpp.read_bytes() == made_in_python.read_bytes()


# --------------------------------------------------------------------------
# errors


def test_error_thrown_by_the_other_module(lattice):
    with pytest.raises(otswap.Error) as raised:
        cross.throw_config_error()
    with pytest.raises(otswap.Error) as same:
        otswap.reconstruct_box(lattice, config=otswap.Config(n_realizations=0))
    assert str(raised.value) == str(same.value) == \
        "the number of realizations is 0; it must be at least 1"


def test_broken_results_give_the_library_messages(tmp_path, table, bias):
    result = cross.broken_result("n_realizations")
    with pytest.raises(otswap.Error, match=full_match(
            "the result holds nRealizations = 2 and config.nRealizations = 3; they must be equal")):
        otswap.io.write_displacement_mean(tmp_path / "mean.dat", result)

    result = cross.broken_result("mean_displacement")
    with pytest.raises(otswap.Error, match=full_match(
            "cannot write meanDisplacement: it holds 3597 entries where 3600 are expected")):
        otswap.io.write_displacement_mean(tmp_path / "mean.dat", result)
    with pytest.raises(otswap.Error, match=full_match(
            "the result is malformed: meanDisplacement and validRealizations hold 3597 and 1200 "
            "entries for 1200 objects")):
        otswap.real_space_lightcone(result, distances=table, bias=bias, sigma=SIGMA,
                                    correction=QUIET_CORRECTION)

    result = cross.broken_result("displacement")
    with pytest.raises(otswap.Error, match=full_match(
            "the result is malformed: displacement holds 7197 entries for 2 realizations of "
            "1200 objects")):
        otswap.reject_mask_crossings(result, cross.make_mask())

    result = cross.broken_result("tracers")
    with pytest.raises(otswap.Error, match=full_match(
            "the result is malformed: tracers holds 3597 entries for 1200 objects; the arcs start "
            "at the tracers, which a result returned by reconstructLightcone carries")):
        otswap.reject_mask_crossings(result, cross.make_mask())


def test_getters_of_a_malformed_result_raise():
    result = cross.broken_result("mean_displacement")
    with pytest.raises(otswap.Error, match=full_match(
            "the result is malformed: mean_displacement holds 3597 entries for 1200 objects")):
        result.mean_displacement
    assert result.tracers.shape == (1200, 3)

    result = cross.broken_result("displacement")
    with pytest.raises(otswap.Error, match=full_match(
            "the result is malformed: displacement holds 7197 entries for 2 realizations of "
            "1200 objects")):
        result.displacement
    assert result.valid.shape == (2, 1200)


# --------------------------------------------------------------------------
# the state tables


def test_state_tables_hold_no_entry_after_the_object_goes(sky, table, bias):
    results, catalogs = otswap._otswap._result_state, otswap._otswap._catalog_state
    gc.collect()
    before = len(results), len(catalogs)

    result = lightcone_result(sky, table)
    assert result in results and len(results) == before[0] + 1
    catalog = otswap.real_space_lightcone(result, distances=table, bias=bias, sigma=SIGMA,
                                          correction=QUIET_CORRECTION)
    assert catalog in catalogs and len(catalogs) == before[1] + 1
    sky_deg, catalog_deg = result.tracers_sky, catalog.sky
    assert result.angle_unit == catalog.angle_unit == "deg"

    del result, catalog
    gc.collect()
    assert (len(results), len(catalogs)) == before
    assert sky_deg.shape == catalog_deg.shape == (1200, 3)


def test_box_results_have_no_entry(box_result):
    assert box_result not in otswap._otswap._result_state


@pytest.mark.parametrize("make", [
    lambda box, lc: lc,
    lambda box, lc: cross.small_lightcone_result(),
    lambda box, lc: cross.small_box_catalog(),
    lambda box, lc: otswap.real_space_box(box, axis=0, redshift=0.5, distances=cross.make_table(),
                                          bias=1.5, sigma=SIGMA, correction=QUIET_CORRECTION),
])
def test_no_new_attribute(make, box_result, lc_result):
    target = make(box_result, lc_result)
    with pytest.raises(AttributeError):
        target.note = "x"
    assert not hasattr(target, "__dict__")
