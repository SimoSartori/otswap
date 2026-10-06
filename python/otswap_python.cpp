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
 *  @file python/otswap_python.cpp
 *
 *  @brief nanobind bindings of include/otswap/OT.h, as specified by
 *  python/otswap/__init__.pyi. The io layer is not bound.
 *
 *  Every argument is taken as a Python object and converted here, so that
 *  any invalid argument raises otswap.Error. Arrays are converted with
 *  numpy.asarray to C-contiguous float64 and their shapes checked; the C++
 *  library then works on its own copy. Right ascension and declination are
 *  converted from degrees here, so the library only ever sees radians.
 *  Result arrays are read-only views on the C++ result, kept alive by the
 *  Python object that holds it.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cstdint>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "otswap/OT.h"

namespace nb = nanobind;
using namespace nb::literals;

// valid_realizations is exposed as uint32 and valid as bool over the uint8
// bytes of the C++ result.
static_assert(sizeof(unsigned) == 4, "valid_realizations is exposed as uint32");
static_assert(sizeof(bool) == 1 && sizeof(std::uint8_t) == 1, "valid is exposed as bool");

namespace {

  // pi/180 as a double: the factor numpy.deg2rad multiplies by, so that an
  // angle converted here and one converted with numpy.deg2rad are the same
  // double.
  constexpr double kDegToRad = 3.14159265358979323846 / 180.;

  using CArray = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;

  std::string repr_of (nb::handle value)
  {
    return nb::repr(value).c_str();
  }

  std::string shape_of (const CArray& a)
  {
    std::string s = "(";
    for (std::size_t i = 0; i < a.ndim(); ++i)
      s += (i ? ", " : "") + std::to_string(a.shape(i));
    return s + (a.ndim() == 1 ? ",)" : ")");
  }

  // ----------------------------------------------------- argument conversion

  // The argument as a C-contiguous float64 numpy array of the same shape.
  nb::object float64_array (nb::handle value, const std::string& name)
  {
    const nb::module_ np = nb::module_::import_("numpy");
    try {
      return np.attr("asarray")(value, "dtype"_a = np.attr("float64"), "order"_a = "C");
    }
    catch (nb::python_error& e) {
      throw otswap::Error(name + " cannot be converted to a float64 array: " + nb::str(e.value()).c_str());
    }
  }

  CArray as_float64 (nb::handle value, const std::string& name)
  {
    return nb::cast<CArray>(float64_array(value, name));
  }

  // An (N, 3) array, copied into the flat layout the library takes.
  std::vector<double> rows3 (nb::handle value, const std::string& name)
  {
    const CArray a = as_float64(value, name);
    if (a.ndim() != 2 || a.shape(1) != 3)
      throw otswap::Error(name + " must have shape (N, 3); it has shape " + shape_of(a));
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  // A one-dimensional array.
  std::vector<double> column (nb::handle value, const std::string& name)
  {
    const CArray a = as_float64(value, name);
    if (a.ndim() != 1)
      throw otswap::Error(name + " must be one-dimensional; it has shape " + shape_of(a));
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  double to_double (nb::handle value, const std::string& name)
  {
    try {
      return nb::cast<double>(value);
    }
    catch (nb::cast_error&) {
      throw otswap::Error(name + " must be a real number; got " + repr_of(value));
    }
  }

  unsigned to_unsigned (nb::handle value, const std::string& name)
  {
    long long v;
    try {
      v = nb::cast<long long>(value);
    }
    catch (nb::cast_error&) {
      throw otswap::Error(name + " must be an integer; got " + repr_of(value));
    }
    if (v < 0 || v > 0xffffffffLL)
      throw otswap::Error(name + " is " + std::to_string(v) + "; it must lie in [0, 4294967295]");
    return (unsigned)v;
  }

  std::string to_string (nb::handle value, const std::string& name)
  {
    if (!nb::isinstance<nb::str>(value))
      throw otswap::Error(name + " must be a str; got " + repr_of(value));
    return nb::str(value).c_str();
  }

  // An instance of a bound class; the reference stays valid while the
  // caller holds the argument.
  template <typename T>
  T& instance (nb::handle value, const std::string& name, const std::string& type)
  {
    if (!nb::isinstance<T>(value))
      throw otswap::Error(name + " must be an otswap." + type + "; got " + repr_of(value));
    return nb::cast<T&>(value);
  }

  bool in_degrees (nb::handle unit)
  {
    if (nb::isinstance<nb::str>(unit)) {
      const std::string u = nb::str(unit).c_str();
      if (u == "deg") return true;
      if (u == "rad") return false;
    }
    throw otswap::Error("angle_unit must be \"deg\" or \"rad\"; got " + repr_of(unit));
  }

  // Right ascension and declination of a flat sky array, from degrees to
  // radians; redshifts are left alone.
  void sky_to_radians (std::vector<double>& sky)
  {
    for (std::size_t i = 0; i + 2 < sky.size(); i += 3) {
      sky[i]   *= kDegToRad;
      sky[i+1] *= kDegToRad;
    }
  }

  otswap::Config make_config (nb::handle nRealizations, nb::handle convergence,
                              nb::handle seed, nb::handle cellSize)
  {
    otswap::Config config;
    config.nRealizations = to_unsigned(nRealizations, "n_realizations");
    config.convergence = to_double(convergence, "convergence");
    config.seed = to_unsigned(seed, "seed");
    config.cellSize = to_double(cellSize, "cell_size");
    return config;
  }

  // --------------------------------------------------------- returned arrays

  std::vector<std::size_t> shape_vector (const CArray& a)
  {
    std::vector<std::size_t> s(a.ndim());
    for (std::size_t i = 0; i < a.ndim(); ++i) s[i] = a.shape(i);
    return s;
  }

  // A new numpy array that owns values.
  template <typename T>
  nb::ndarray<nb::numpy, T> owned (std::vector<T>&& values, const std::vector<std::size_t>& shape)
  {
    auto* held = new std::vector<T>(std::move(values));
    nb::capsule owner(held, [] (void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T>(held->data(), shape.size(), shape.data(), owner);
  }

  // A DistanceTable lookup applied to every element, keeping the shape.
  template <typename F>
  nb::ndarray<nb::numpy, double> elementwise (nb::handle value, const std::string& name, F lookup)
  {
    const CArray a = as_float64(value, name);
    std::vector<double> out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) out[i] = lookup(a.data()[i]);
    return owned(std::move(out), shape_vector(a));
  }

  // A read-only view on an array of a Result, which keeps the Python object
  // holding the Result alive.
  template <typename T>
  nb::ndarray<nb::numpy, const T> view (nb::handle owner, const T* data,
                                        std::initializer_list<std::size_t> shape)
  {
    return nb::ndarray<nb::numpy, const T>(data, shape, owner);
  }

  using ResultHandle = nb::handle_t<otswap::Result>;

  const otswap::Result& result_of (ResultHandle self)
  {
    return nb::cast<const otswap::Result&>(self);
  }

}

// Every argument is declared .none(), so that None reaches the conversions
// above and raises otswap.Error like any other invalid value.
NB_MODULE(_otswap, m)
{
  m.doc() = "Compiled part of the otswap package; import otswap instead.";

  nb::exception<otswap::Error> error(m, "Error", PyExc_RuntimeError);
  error.attr("__doc__") = "Raised by every otswap function on invalid input or failure.";

  // ------------------------------------------------------------------ Result

  nb::class_<otswap::Result>(m, "Result", "Displacement field produced by a reconstruction.")
    .def_prop_ro("n_objects", [] (ResultHandle self) { return result_of(self).nObjects; })
    .def_prop_ro("n_realizations", [] (ResultHandle self) { return result_of(self).nRealizations; })
    .def_prop_ro("displacement", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<double>(self, r.displacement.data(), {r.nRealizations, r.nObjects, 3});
      }, "Displacements, shape (n_realizations, n_objects, 3).")
    .def_prop_ro("matched_random", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<double>(self, r.matchedRandom.data(), {r.nRealizations, r.nObjects, 3});
      }, "Cartesian position of the random matched to each tracer, shape "
         "(n_realizations, n_objects, 3).")
    .def_prop_ro("valid", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.valid.data()),
                          {r.nRealizations, r.nObjects});
      }, "Whether each displacement is valid, shape (n_realizations, n_objects).")
    .def_prop_ro("valid_realizations", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(r.validRealizations.data()),
                                   {r.nObjects});
      }, "Number of valid realizations per tracer, shape (n_objects,).")
    .def_prop_ro("mean_displacement", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<double>(self, r.meanDisplacement.data(), {r.nObjects, 3});
      }, "Displacement averaged over the valid realizations, shape (n_objects, 3).")
    .def_prop_ro("filtered_nside", [] (ResultHandle self) { return result_of(self).filteredNside; },
                 "NSIDE of the mask the result was filtered against, 0 if none.");

  // ----------------------------------------------------------- DistanceTable

  nb::class_<otswap::DistanceTable>(m, "DistanceTable",
      "Sampled relation between redshift, comoving distance and, optionally, the linear "
      "growth rate, interpolated in both directions. Immutable.")
    .def_static("flat",
        [] (nb::handle omegaM, nb::handle h, nb::handle w0, nb::handle wa,
            nb::handle zMin, nb::handle zMax, nb::handle nSamples) {
          return otswap::DistanceTable(to_double(omegaM, "omega_m"), to_double(h, "h"),
                                       to_double(w0, "w0"), to_double(wa, "wa"),
                                       to_double(zMin, "z_min"), to_double(zMax, "z_max"),
                                       to_unsigned(nSamples, "n_samples"));
        },
        "omega_m"_a.none(), "h"_a.none(), "w0"_a.none() = -1., "wa"_a.none() = 0.,
        "z_min"_a.none() = 0., "z_max"_a.none() = 10., "n_samples"_a.none() = 50000,
        "Tabulate a flat cosmology with dark energy w(a) = w0 + wa (1 - a).")
    .def_static("from_table",
        [] (nb::handle redshift, nb::handle distance, nb::handle growthRate) {
          std::vector<double> z = column(redshift, "redshift");
          std::vector<double> d = column(distance, "distance");
          std::vector<double> f;
          if (!growthRate.is_none()) f = column(growthRate, "growth_rate");
          return otswap::DistanceTable(std::move(z), std::move(d), std::move(f));
        },
        "redshift"_a.none(), "distance"_a.none(), "growth_rate"_a = nb::none(),
        "Use a table computed by the caller.")
    .def_prop_ro("min_redshift", &otswap::DistanceTable::minRedshift)
    .def_prop_ro("max_redshift", &otswap::DistanceTable::maxRedshift)
    .def_prop_ro("has_growth_rate", &otswap::DistanceTable::hasGrowthRate)
    .def("distance_at", [] (const otswap::DistanceTable& t, nb::handle z) {
          return elementwise(z, "z", [&t] (const double v) { return t.distanceAt(v); });
        }, "z"_a.none(), "Comoving distance at each redshift. Raises outside the table.")
    .def("redshift_at", [] (const otswap::DistanceTable& t, nb::handle distance) {
          return elementwise(distance, "distance", [&t] (const double v) { return t.redshiftAt(v); });
        }, "distance"_a.none(), "Redshift at each comoving distance. Raises outside the table.")
    .def("growth_rate_at", [] (const otswap::DistanceTable& t, nb::handle z) {
          return elementwise(z, "z", [&t] (const double v) { return t.growthRateAt(v); });
        }, "z"_a.none(), "Linear growth rate at each redshift.");

  m.def("to_cartesian",
      [] (nb::handle sky, nb::handle distances, nb::handle angleUnit) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        std::vector<double> s = rows3(sky, "sky");
        if (degrees) sky_to_radians(s);
        const std::size_t n = s.size() / 3;
        return owned(otswap::toCartesian(s, table), {n, 3});
      },
      "sky"_a.none(), "distances"_a.none(), nb::kw_only(), "angle_unit"_a.none(),
      "Convert sky coordinates, shape (N, 3), to Cartesian ones in Mpc/h.");

  // -------------------------------------------------------------------- Mask

  nb::class_<otswap::Mask>(m, "Mask",
                           "HEALPix mask read from a FITS file: a pixel is observed when its "
                           "value is greater than 0. Immutable.")
    .def("__init__", [] (otswap::Mask* self, nb::handle fitsFile) {
          new (self) otswap::Mask(to_string(fitsFile, "fits_file"));
        }, "fits_file"_a.none())
    .def_prop_ro("nside", &otswap::Mask::nside)
    .def_prop_ro("sky_area_deg2", &otswap::Mask::skyAreaDeg2,
                 "Area covered by the observed pixels, in square degrees.")
    .def("allows",
        [] (const otswap::Mask& mask, nb::handle ra, nb::handle dec, nb::handle angleUnit) {
          const bool degrees = in_degrees(angleUnit);
          const nb::object r0 = float64_array(ra, "ra"), d0 = float64_array(dec, "dec");
          const nb::module_ np = nb::module_::import_("numpy");
          nb::object pair;
          try {
            pair = np.attr("broadcast_arrays")(r0, d0);
          }
          catch (nb::python_error& e) {
            throw otswap::Error(std::string("ra and dec cannot be broadcast together: ") +
                                nb::str(e.value()).c_str());
          }
          const CArray r = as_float64(pair[0], "ra"), d = as_float64(pair[1], "dec");
          auto* held = new bool[r.size() ? r.size() : 1];
          nb::capsule owner(held, [] (void* p) noexcept { delete[] static_cast<bool*>(p); });
          for (std::size_t i = 0; i < r.size(); ++i)
            held[i] = degrees ? mask.allows(r.data()[i] * kDegToRad, d.data()[i] * kDegToRad)
                              : mask.allows(r.data()[i], d.data()[i]);
          const std::vector<std::size_t> shape = shape_vector(r);
          return nb::ndarray<nb::numpy, bool>(held, shape.size(), shape.data(), owner);
        },
        "ra"_a.none(), "dec"_a.none(), nb::kw_only(), "angle_unit"_a.none(),
        "Whether each direction falls in an observed pixel.");

  // ---------------------------------------------------------- Reconstruction

  m.def("reconstruct_box",
      [] (nb::handle tracers, nb::handle randoms, nb::handle mps, nb::handle nRealizations,
          nb::handle convergence, nb::handle seed, nb::handle cellSize) {
        const std::vector<double> t = rows3(tracers, "tracers");
        std::vector<double> r;
        if (!randoms.is_none()) {
          r = rows3(randoms, "randoms");
          if (r.empty())
            throw otswap::Error("randoms is empty; omit it to have the randoms drawn "
                                "in the bounding box of the tracers");
        }
        const double separation = to_double(mps, "mps");
        const otswap::Config config = make_config(nRealizations, convergence, seed, cellSize);
        nb::gil_scoped_release released;
        return otswap::reconstructBox(t, r, separation, config);
      },
      "tracers"_a.none(), "randoms"_a = nb::none(), nb::kw_only(), "mps"_a.none(),
      "n_realizations"_a.none() = 1, "convergence"_a.none() = 1.e-3, "seed"_a.none() = 0,
      "cell_size"_a.none() = 4.,
      "Reconstruct in box geometry, with a constant mean particle separation.");

  m.def("reconstruct_lightcone",
      [] (nb::handle tracersSky, nb::handle randomsSky, nb::handle skyAreaDeg2, nb::handle nBins,
          nb::handle distances, nb::handle angleUnit, nb::handle tracers, nb::handle randoms,
          nb::handle nRealizations, nb::handle convergence, nb::handle seed, nb::handle cellSize) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        std::vector<double> ts = rows3(tracersSky, "tracers_sky");
        std::vector<double> rs = rows3(randomsSky, "randoms_sky");
        if (degrees) { sky_to_radians(ts); sky_to_radians(rs); }
        if (tracers.is_none() != randoms.is_none())
          throw otswap::Error("give both tracers and randoms, or neither");
        const double area = to_double(skyAreaDeg2, "sky_area_deg2");
        const unsigned bins = to_unsigned(nBins, "n_bins");
        const otswap::Config config = make_config(nRealizations, convergence, seed, cellSize);
        if (tracers.is_none()) {
          nb::gil_scoped_release released;
          return otswap::reconstructLightcone(ts, rs, area, bins, table, config);
        }
        const std::vector<double> t = rows3(tracers, "tracers");
        const std::vector<double> r = rows3(randoms, "randoms");
        nb::gil_scoped_release released;
        return otswap::reconstructLightcone(t, r, ts, rs, area, bins, table, config);
      },
      "tracers_sky"_a.none(), "randoms_sky"_a.none(), nb::kw_only(), "sky_area_deg2"_a.none(),
      "n_bins"_a.none(), "distances"_a.none(), "angle_unit"_a.none(), "tracers"_a = nb::none(),
      "randoms"_a = nb::none(), "n_realizations"_a.none() = 1, "convergence"_a.none() = 1.e-3,
      "seed"_a.none() = 0, "cell_size"_a.none() = 4.,
      "Reconstruct in lightcone geometry.");

  m.def("reject_mask_crossings",
      [] (nb::handle result, nb::handle mask, nb::handle maxForbiddenPixels) {
        otswap::Result& r = instance<otswap::Result>(result, "result", "Result");
        const otswap::Mask& k = instance<otswap::Mask>(mask, "mask", "Mask");
        const unsigned limit = to_unsigned(maxForbiddenPixels, "max_forbidden_pixels");
        nb::gil_scoped_release released;
        otswap::rejectMaskCrossings(r, k, limit);
      },
      "result"_a.none(), "mask"_a.none(), "max_forbidden_pixels"_a.none() = 0,
      "Mark as invalid the displacements whose path crosses more than max_forbidden_pixels "
      "distinct unobserved pixels of the mask. Updates result in place.");
}
