// SPDX-License-Identifier: GPL-2.0-or-later
/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * This program is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as   *
 * published by the Free Software Foundation; either version 2 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU General Public License for more details.                     *
 *                                                                  *
 * You should have received a copy of the GNU General Public        *
 * License along with this program; if not, write to the Free       *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,      *
 * Boston, MA 02110-1301 USA.                                       *
 ********************************************************************/

/**
 *  @file tests/crossmodule/crosstest.cpp
 *
 *  @brief _otswap_crosstest, a nanobind module for the tests only, built
 *  as another package would build one on top of otswap: it binds no type,
 *  and takes and returns the shared types as their C++ types, which
 *  nanobind finds in the registry otswap filled at its import. Its
 *  functions read what C++ sees of an object made in Python, and make
 *  objects in C++ for Python to read. tests/test_crossmodule.py uses it.
 *  Built with the tests and the Python bindings; never installed.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "otswap/OT.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

  constexpr double kPi = 3.14159265358979323846;

  const char* verbosity_name (const otswap::Verbosity v)
  {
    if (v == otswap::Verbosity::Silent) return "silent";
    if (v == otswap::Verbosity::Detailed) return "detailed";
    return "normal";
  }

  /// Silent, writing nothing.
  otswap::Config quiet_config (const unsigned seed)
  {
    otswap::Config config;
    config.nRealizations = 2;
    config.seed = seed;
    config.verbosity = otswap::Verbosity::Silent;
    config.log = nullptr;
    return config;
  }

  otswap::CorrectionConfig quiet_correction ()
  {
    otswap::CorrectionConfig config;
    config.verbosity = otswap::Verbosity::Silent;
    config.log = nullptr;
    return config;
  }

  /// A 6 x 6 x 6 lattice, 10 Mpc/h apart.
  std::vector<double> lattice ()
  {
    std::vector<double> t;
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        for (int k = 0; k < 6; ++k)
          t.insert(t.end(), {10. * i, 10. * j, 10. * k});
    return t;
  }

  /// n points of a fixed low-discrepancy sequence, from index first, on the
  /// patch 10 < RA < 40 deg, -15 < Dec < 15 deg, 0.3 < z < 0.6; radians.
  std::vector<double> sky_sample (const std::size_t n, const std::size_t first)
  {
    const double a1 = 0.8191725133961645, a2 = 0.6710436067037893, a3 = 0.5497004779019703;
    const double s = std::sin(15. * kPi / 180.);
    std::vector<double> sky;
    for (std::size_t i = first; i < first + n; ++i) {
      const double u = std::fmod(0.5 + a1 * i, 1.), v = std::fmod(0.5 + a2 * i, 1.),
                   w = std::fmod(0.5 + a3 * i, 1.);
      sky.insert(sky.end(), {(10. + 30. * u) * kPi / 180., std::asin(-s + 2. * s * v), 0.3 + 0.3 * w});
    }
    return sky;
  }

  /// The area of that patch, in square degrees.
  double patch_area ()
  {
    const double deg = 180. / kPi;
    return 30. / deg * 2. * std::sin(15. / deg) * deg * deg;
  }

  otswap::DistanceTable make_table ()
  {
    return otswap::DistanceTable(0.3, 0.7, -1., 0., 0., 1.5, 4000);
  }

  otswap::BiasTable make_bias_table ()
  {
    return otswap::BiasTable{{0.2, 0.8}, {1.2, 1.8}};
  }

  otswap::Result small_box_result ()
  {
    return otswap::reconstructBox(lattice(), {}, quiet_config(7));
  }

  otswap::Result small_lightcone_result ()
  {
    return otswap::reconstructLightcone(sky_sample(1200, 0), sky_sample(3600, 5000), patch_area(), 2,
                                        make_table(), quiet_config(9));
  }

  nb::object row_or_none (const std::vector<double>& flat)
  {
    if (flat.size() < 3) return nb::none();
    return nb::make_tuple(flat[0], flat[1], flat[2]);
  }

}

NB_MODULE(_otswap_crosstest, m)
{
  m.doc() = "Test-only module that uses otswap's types from another nanobind module.";

  // registers the shared types, which this module does not bind
  const nb::module_ otswap_module = nb::module_::import_("otswap");

  // otswap's own translator may not recognise an otswap::Error thrown in
  // this module (with libc++ the two modules do not share its type
  // information), so this module translates its own; otswap.Error is held
  // for the life of the process
  nb::register_exception_translator(
    [] (const std::exception_ptr& p, void* error) {
      try {
        std::rethrow_exception(p);
      }
      catch (const otswap::Error& e) {
        PyErr_SetString(static_cast<PyObject*>(error), e.what());
      }
    }, nb::object(otswap_module.attr("Error")).release().ptr());

  // ------------------------------------------------------------ configuration

  m.def("default_config", [] { return otswap::Config(); });

  m.def("changed_config", [] (const otswap::Config& c) {
    otswap::Config out = c;
    out.nRealizations += 1;
    out.convergence *= 2.;
    out.seed += 1;
    out.cellSize += 1.;
    out.rejectCrossings = !out.rejectCrossings;
    out.maxUnobservedPixelsCrossed += 1;
    out.verbosity = c.verbosity == otswap::Verbosity::Silent ? otswap::Verbosity::Normal
                  : c.verbosity == otswap::Verbosity::Normal ? otswap::Verbosity::Detailed
                  : otswap::Verbosity::Silent;
    return out;
  });

  m.def("config_fields", [] (const otswap::Config& c) {
    nb::dict d;
    d["n_realizations"] = c.nRealizations;
    d["convergence"] = c.convergence;
    d["seed"] = c.seed;
    d["cell_size"] = c.cellSize;
    d["reject_crossings"] = c.rejectCrossings;
    d["max_unobserved_pixels_crossed"] = c.maxUnobservedPixelsCrossed;
    d["verbosity"] = verbosity_name(c.verbosity);
    d["log_is_clog"] = c.log == &std::clog;
    return d;
  });

  m.def("default_correction", [] { return otswap::CorrectionConfig(); });

  m.def("changed_correction", [] (const otswap::CorrectionConfig& c) {
    otswap::CorrectionConfig out = c;
    out.weightByRealizations = !out.weightByRealizations;
    out.verbosity = c.verbosity == otswap::Verbosity::Silent ? otswap::Verbosity::Normal
                  : c.verbosity == otswap::Verbosity::Normal ? otswap::Verbosity::Detailed
                  : otswap::Verbosity::Silent;
    return out;
  });

  m.def("correction_fields", [] (const otswap::CorrectionConfig& c) {
    nb::dict d;
    d["weight_by_realizations"] = c.weightByRealizations;
    d["verbosity"] = verbosity_name(c.verbosity);
    d["log_is_clog"] = c.log == &std::clog;
    return d;
  });

  m.def("cut_bounds", [] (const otswap::RedshiftCut& c) { return std::make_pair(c.min, c.max); });

  m.def("make_cut", [] (const double min, const double max) { return otswap::RedshiftCut{min, max}; });

  m.def("box_corners", [] (const otswap::BoxLimits& b) {
    return nb::make_tuple(nb::make_tuple(b.lower[0], b.lower[1], b.lower[2]),
                          nb::make_tuple(b.upper[0], b.upper[1], b.upper[2]));
  });

  m.def("make_box", [] {
    otswap::BoxLimits b;
    b.lower = {0., 0., 0.};
    b.upper = {50., 60., 70.};
    return b;
  });

  // ---------------------------------------------------------- tables, mask

  m.def("mask_nside", [] (const otswap::Mask& mask) { return mask.nside(); });

  // NSIDE 2, every pixel observed but the first
  m.def("make_mask", [] {
    std::vector<double> values(48, 1.);
    values[0] = 0.;
    return otswap::Mask(values, otswap::PixelOrdering::Ring);
  });

  m.def("distance_at", [] (const otswap::DistanceTable& t, const double z) { return t.distanceAt(z); });

  m.def("make_table", &make_table);

  m.def("profile_bins", [] (const otswap::MpsProfile& p) { return p.redshift.size(); });

  m.def("bias_nodes", [] (const otswap::BiasTable& b) { return nb::make_tuple(b.redshift, b.bias); });

  m.def("make_bias_table", &make_bias_table);

  // ------------------------------------------------------ Result, catalogue

  m.def("result_fields", [] (const otswap::Result& r) {
    nb::dict d;
    d["n_objects"] = r.nObjects;
    d["n_realizations"] = r.nRealizations;
    d["geometry"] = r.geometry == otswap::Geometry::Box ? "box" : "lightcone";
    d["seed"] = r.config.seed;
    d["has_distances"] = (bool)r.distances;
    d["first_tracers_sky"] = row_or_none(r.tracersSky);
    d["log_is_clog"] = r.config.log == &std::clog;
    return d;
  });

  m.def("small_box_result", &small_box_result);

  m.def("small_lightcone_result", &small_lightcone_result);

  m.def("catalog_fields", [] (const otswap::RealSpaceCatalog& c) {
    nb::dict d;
    d["n_objects"] = c.nObjects;
    d["geometry"] = c.geometry == otswap::Geometry::Box ? "box" : "lightcone";
    d["axis"] = c.axis;
    d["sigma"] = c.sigma;
    d["first_sky"] = row_or_none(c.sky);
    std::vector<std::size_t> counts(4, 0);
    for (const otswap::CorrectionStatus s : c.status) ++counts[(std::size_t)s];
    d["status_counts"] = nb::make_tuple(counts[0], counts[1], counts[2], counts[3]);
    return d;
  });

  m.def("small_box_catalog", [] {
    return otswap::realSpaceBox(small_box_result(), 2, 0.5, make_table(), 1.5, 10., quiet_correction());
  });

  m.def("small_lightcone_catalog", [] {
    return otswap::realSpaceLightcone(small_lightcone_result(), make_table(), make_bias_table(), 10.,
                                      quiet_correction());
  });

  // ----------------------------------------------------------------- errors

  m.def("throw_config_error", [] {
    otswap::Config config = quiet_config(7);
    config.nRealizations = 0;
    (void)otswap::reconstructBox(lattice(), {}, config);
  });

  // small_lightcone_result() with one member broken, as otswap's own
  // constructors cannot make it
  m.def("broken_result", [] (const std::string& what) {
    otswap::Result r = small_lightcone_result();
    if (what == "n_realizations") r.config.nRealizations += 1;
    else if (what == "mean_displacement") r.meanDisplacement.resize(r.meanDisplacement.size() - 3);
    else if (what == "displacement") r.displacement.resize(r.displacement.size() - 3);
    else if (what == "tracers") r.tracers.resize(r.tracers.size() - 3);
    else throw otswap::Error("broken_result: unknown member " + what);
    return r;
  }, "what"_a);
}
