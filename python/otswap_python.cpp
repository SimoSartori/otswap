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
 *  Python object that holds it. The selection report of a lightcone
 *  reconstruction is printed from Python, to sys.stdout, which a notebook
 *  shows, rather than by the library to the process's stderr.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <sstream>
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
// uncorrected is exposed as int64 over the size_t indices of the C++ catalogue.
static_assert(sizeof(std::size_t) == 8, "uncorrected is exposed as int64");

namespace {

  // pi/180 as a double: the factor numpy.deg2rad multiplies by, so that an
  // angle converted here and one converted with numpy.deg2rad are the same
  // double.
  constexpr double kDegToRad = 3.14159265358979323846 / 180.;

  // 180/pi as a double: the factor numpy.rad2deg multiplies by.
  constexpr double kRadToDeg = 180. / 3.14159265358979323846;

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

  bool to_bool (nb::handle value, const std::string& name)
  {
    if (!PyBool_Check(value.ptr()))
      throw otswap::Error(name + " must be True or False; got " + repr_of(value));
    return value.ptr() == Py_True;
  }

  // A one-dimensional array of counts: each entry an integer in
  // [0, 4294967295].
  std::vector<unsigned> counts (nb::handle value, const std::string& name)
  {
    const std::vector<double> c = column(value, name);
    std::vector<unsigned> out(c.size());
    for (std::size_t i = 0; i < c.size(); ++i) {
      if (!(c[i] >= 0. && c[i] <= 4294967295.) || c[i] != std::floor(c[i]))
        throw otswap::Error(name + " holds " + std::to_string(c[i]) + " at entry " +
                            std::to_string(i) + "; every entry must be an integer in [0, 4294967295]");
      out[i] = (unsigned)c[i];
    }
    return out;
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

  /// Right ascension and declination of a flat sky array, from degrees to
  /// radians; redshifts are left alone. A finite declination outside
  /// [-90, 90] is refused first, in degrees; 90 degrees converts to the
  /// double pi/2, so the library's own check in radians agrees.
  void sky_to_radians (std::vector<double>& sky, const std::string& name)
  {
    for (std::size_t i = 0; i + 2 < sky.size(); i += 3)
      if (std::isfinite(sky[i+1]) && std::fabs(sky[i+1]) > 90.)
        throw otswap::Error(name + " holds a declination of " + std::to_string(sky[i+1]) +
                            " degrees at object " + std::to_string(i / 3) + ", outside [-90, 90]");
    for (std::size_t i = 0; i + 2 < sky.size(); i += 3) {
      sky[i]   *= kDegToRad;
      sky[i+1] *= kDegToRad;
    }
  }

  /// Right ascension and declination of a flat sky array, from radians to
  /// degrees; redshifts are left alone. A right ascension just below 2 pi
  /// can round to 360 in the product, and is folded to 0, so that it stays
  /// in [0, 360).
  void sky_to_degrees (std::vector<double>& sky)
  {
    for (std::size_t i = 0; i + 2 < sky.size(); i += 3) {
      sky[i]   *= kRadToDeg;
      sky[i+1] *= kRadToDeg;
      if (sky[i] >= 360.) sky[i] = 0.;
    }
  }

  // None, or a (min, max) pair.
  otswap::RedshiftCut to_cut (nb::handle value)
  {
    otswap::RedshiftCut cut;
    if (value.is_none()) return cut;
    const std::vector<double> c = column(value, "redshift_cut");
    if (c.size() != 2)
      throw otswap::Error("redshift_cut must be a (min, max) pair; it holds " +
                          std::to_string(c.size()) + " values");
    cut.min = c[0];
    cut.max = c[1];
    return cut;
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

  // ------------------------------------------------------------- warnings

  // otswap.ExtrapolationWarning, created with the module.
  PyObject* extrapolation_warning = nullptr;

  // Raise one ExtrapolationWarning for a call that extrapolated b(z) at
  // count of total redshifts. Needs the GIL.
  void warn_extrapolation (const std::size_t count, const std::size_t total,
                           const std::vector<double>& biasRedshift)
  {
    if (count == 0) return;
    std::ostringstream message;
    message << "b(z) extrapolated at " << count << " of " << total
            << " redshifts, outside the bias table's range [" << biasRedshift.front() << ", "
            << biasRedshift.back() << "]";
    if (PyErr_WarnEx(extrapolation_warning, message.str().c_str(), 1) < 0)
      throw nb::python_error();
  }

  using ResultHandle = nb::handle_t<otswap::Result>;
  using CatalogHandle = nb::handle_t<otswap::RealSpaceCatalog>;

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

  extrapolation_warning = PyErr_NewExceptionWithDoc(
    "otswap._otswap.ExtrapolationWarning",
    "Issued when b(z) is extrapolated beyond its table.", PyExc_UserWarning, nullptr);
  if (extrapolation_warning == nullptr) throw nb::python_error();
  m.attr("ExtrapolationWarning") = nb::handle(extrapolation_warning);

  // --------------------------------------------------------- SelectionCounts

  nb::class_<otswap::SelectionCounts>(m, "SelectionCounts",
      "What a lightcone reconstruction left out, and why. Read-only.")
    .def_prop_ro("redshift_cut", [] (const otswap::SelectionCounts& c) -> nb::object {
        if (!c.redshiftCutApplied) return nb::none();
        return nb::make_tuple(c.redshiftCut.min, c.redshiftCut.max);
      }, "The (min, max) of the redshift cut, or None when no cut was applied.")
    .def_ro("mask_applied", &otswap::SelectionCounts::maskApplied)
    .def_ro("tracers", &otswap::SelectionCounts::tracers)
    .def_ro("tracers_outside_redshift_cut", &otswap::SelectionCounts::tracersOutsideRedshiftCut)
    .def_ro("tracers_outside_mask", &otswap::SelectionCounts::tracersOutsideMask)
    .def_ro("tracers_outside_both", &otswap::SelectionCounts::tracersOutsideBoth)
    .def_ro("randoms", &otswap::SelectionCounts::randoms)
    .def_ro("randoms_outside_redshift_cut", &otswap::SelectionCounts::randomsOutsideRedshiftCut)
    .def_ro("randoms_outside_mask", &otswap::SelectionCounts::randomsOutsideMask)
    .def_ro("randoms_outside_both", &otswap::SelectionCounts::randomsOutsideBoth)
    .def_prop_ro("max_unobserved_pixels_crossed", [] (const otswap::SelectionCounts& c) -> nb::object {
        if (!c.crossingsRejected) return nb::none();
        return nb::int_(c.maxUnobservedPixelsCrossed);
      }, "The threshold of the mask filter the call applied, or None when it applied none.")
    .def_ro("displacements", &otswap::SelectionCounts::displacements)
    .def_ro("displacements_crossing_mask", &otswap::SelectionCounts::displacementsCrossingMask)
    .def("__repr__", [] (const otswap::SelectionCounts& c) {
        std::string text = c.message();
        if (text.empty()) return nb::str("SelectionCounts(no selection)");
        text.pop_back();
        return nb::str(text.c_str());
      });

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
                 "NSIDE of the mask the result was filtered against, 0 if none.")
    .def_prop_ro("outside_redshift_cut", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.outsideRedshiftCut.data()),
                          {r.nObjects});
      }, "Whether each tracer lay outside the redshift cut, shape (n_objects,).")
    .def_prop_ro("outside_mask", [] (ResultHandle self) {
        const otswap::Result& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.outsideMask.data()), {r.nObjects});
      }, "Whether each tracer fell on an unobserved pixel of the mask, shape (n_objects,).")
    .def_prop_ro("selection", [] (const otswap::Result& r) -> const otswap::SelectionCounts& {
        return r.selection;
      }, nb::rv_policy::reference_internal, "What the reconstruction left out, and why.")
    .def_prop_ro("lagrangian_sky", [] (ResultHandle self) -> nb::object {
        const otswap::Result& r = result_of(self);
        if (r.lagrangianSky.empty()) return nb::none();
        return nb::cast(view<double>(self, r.lagrangianSky.data(), {r.nObjects, 3}));
      }, "Sky coordinates of each tracer's mean Lagrangian position, shape (n_objects, 3), in "
         "the angle_unit of the reconstruct_lightcone call; None for a box result and after a "
         "reject_mask_crossings call that rejected any displacement.");

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
        if (degrees) sky_to_radians(s, "sky");
        const std::size_t n = s.size() / 3;
        return owned(otswap::toCartesian(s, table), {n, 3});
      },
      "sky"_a.none(), "distances"_a.none(), nb::kw_only(), "angle_unit"_a.none(),
      "Convert sky coordinates, shape (N, 3), to Cartesian ones in Mpc/h.");

  m.def("to_sky",
      [] (nb::handle cartesian, nb::handle distances, nb::handle angleUnit) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        const std::vector<double> c = rows3(cartesian, "cartesian");
        std::vector<double> sky = otswap::toSky(c, table);
        if (degrees) sky_to_degrees(sky);
        const std::size_t n = sky.size() / 3;
        return owned(std::move(sky), {n, 3});
      },
      "cartesian"_a.none(), "distances"_a.none(), nb::kw_only(), "angle_unit"_a.none(),
      "Convert Cartesian coordinates, shape (N, 3), in Mpc/h, to sky coordinates.");

  // -------------------------------------------------------------------- Mask

  nb::class_<otswap::Mask>(m, "Mask",
                           "HEALPix mask read from a FITS file or built from a full-sky map in "
                           "memory: a pixel is observed when its value is greater than 0. "
                           "Immutable.")
    .def("__init__", [] (otswap::Mask* self, nb::handle fitsFile) {
          new (self) otswap::Mask(to_string(fitsFile, "fits_file"));
        }, "fits_file"_a.none())
    .def_static("from_array",
        [] (nb::handle values, nb::handle nest) {
          const bool nested = to_bool(nest, "nest");
          // The dtype is checked before any conversion, so that complex
          // values, strings or objects are refused rather than cast. A
          // C-contiguous float64 array is then read in place; any other is
          // converted to one.
          const nb::module_ np = nb::module_::import_("numpy");
          nb::object raw;
          try {
            raw = np.attr("asarray")(values);
          }
          catch (nb::python_error& e) {
            throw otswap::Error(std::string("values cannot be converted to an array: ") +
                                nb::str(e.value()).c_str());
          }
          const std::string kind = nb::str(raw.attr("dtype").attr("kind")).c_str();
          if (kind != "b" && kind != "i" && kind != "u" && kind != "f")
            throw otswap::Error("values must hold real numbers; its dtype is " +
                                repr_of(raw.attr("dtype")));
          const CArray a = as_float64(raw, "values");
          if (a.ndim() != 1)
            throw otswap::Error("values must be one-dimensional; it has shape " + shape_of(a));
          return otswap::Mask(a.data(), a.size(),
                              nested ? otswap::PixelOrdering::Nested : otswap::PixelOrdering::Ring);
        },
        "values"_a.none(), "nest"_a.none() = false,
        "A full-sky HEALPix map given in memory, one value per pixel; nest=True for NESTED.")
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
          for (std::size_t i = 0; i < r.size(); ++i) {
            const double dec = d.data()[i];
            if (degrees && std::isfinite(dec) && std::fabs(dec) > 90.)
              throw otswap::Error("the declination is " + std::to_string(dec) +
                                  " degrees; it must lie in [-90, 90]");
            held[i] = degrees ? mask.allows(r.data()[i] * kDegToRad, dec * kDegToRad)
                              : mask.allows(r.data()[i], dec);
          }
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
      [] (nb::handle tracersSky, nb::handle randomsSky, nb::handle skyAreaDeg2, nb::handle mask,
          nb::handle nBins, nb::handle distances, nb::handle angleUnit, nb::handle tracers,
          nb::handle randoms, nb::handle nRealizations, nb::handle convergence, nb::handle seed,
          nb::handle cellSize, nb::handle redshiftCut, nb::handle rejectCrossings,
          nb::handle maxUnobservedPixelsCrossed, nb::handle verbose) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        std::vector<double> ts = rows3(tracersSky, "tracers_sky");
        std::vector<double> rs = rows3(randomsSky, "randoms_sky");
        if (degrees) { sky_to_radians(ts, "tracers_sky"); sky_to_radians(rs, "randoms_sky"); }
        if (tracers.is_none() != randoms.is_none())
          throw otswap::Error("give both tracers and randoms, or neither");
        if (!skyAreaDeg2.is_none() && !mask.is_none())
          throw otswap::Error("give sky_area_deg2 or mask, not both; with a mask the sky area is "
                              "the mask's");
        if (skyAreaDeg2.is_none() && mask.is_none())
          throw otswap::Error("give sky_area_deg2 or mask; one is required");
        const otswap::Mask* footprint = mask.is_none() ? nullptr
                                       : &instance<otswap::Mask>(mask, "mask", "Mask");
        const double area = footprint ? 0. : to_double(skyAreaDeg2, "sky_area_deg2");
        const unsigned bins = to_unsigned(nBins, "n_bins");
        otswap::Config config = make_config(nRealizations, convergence, seed, cellSize);
        config.rejectCrossings = to_bool(rejectCrossings, "reject_crossings");
        config.maxUnobservedPixelsCrossed =
          to_unsigned(maxUnobservedPixelsCrossed, "max_unobserved_pixels_crossed");
        const bool report = to_bool(verbose, "verbose");
        const otswap::RedshiftCut cut = to_cut(redshiftCut);
        std::vector<double> t, r;
        if (!tracers.is_none()) {
          t = rows3(tracers, "tracers");
          r = rows3(randoms, "randoms");
        }

        config.verbose = false;
        otswap::Result result;
        {
          nb::gil_scoped_release released;
          if (tracers.is_none())
            result = footprint ? otswap::reconstructLightcone(ts, rs, *footprint, bins, table, config, cut)
                               : otswap::reconstructLightcone(ts, rs, area, bins, table, config, cut);
          else
            result = footprint ? otswap::reconstructLightcone(t, r, ts, rs, *footprint, bins, table, config, cut)
                               : otswap::reconstructLightcone(t, r, ts, rs, area, bins, table, config, cut);
        }
        if (report) {
          const std::string message = result.selection.message();
          const nb::object print = nb::module_::import_("builtins").attr("print");
          std::size_t start = 0;
          for (std::size_t end; (end = message.find('\n', start)) != std::string::npos; start = end + 1)
            print(nb::str(message.data() + start, end - start));
        }
        if (degrees) sky_to_degrees(result.lagrangianSky);
        return result;
      },
      "tracers_sky"_a.none(), "randoms_sky"_a.none(), nb::kw_only(), "sky_area_deg2"_a = nb::none(),
      "mask"_a = nb::none(), "n_bins"_a.none(), "distances"_a.none(), "angle_unit"_a.none(),
      "tracers"_a = nb::none(), "randoms"_a = nb::none(), "n_realizations"_a.none() = 1,
      "convergence"_a.none() = 1.e-3, "seed"_a.none() = 0, "cell_size"_a.none() = 4.,
      "redshift_cut"_a = nb::none(), "reject_crossings"_a.none() = true,
      "max_unobserved_pixels_crossed"_a.none() = 0, "verbose"_a.none() = true,
      "Reconstruct in lightcone geometry.");

  m.def("reject_mask_crossings",
      [] (nb::handle result, nb::handle mask, nb::handle maxUnobservedPixelsCrossed) {
        otswap::Result& r = instance<otswap::Result>(result, "result", "Result");
        const otswap::Mask& k = instance<otswap::Mask>(mask, "mask", "Mask");
        const unsigned limit = to_unsigned(maxUnobservedPixelsCrossed, "max_unobserved_pixels_crossed");
        nb::gil_scoped_release released;
        otswap::rejectMaskCrossings(r, k, limit);
      },
      "result"_a.none(), "mask"_a.none(), "max_unobserved_pixels_crossed"_a.none() = 0,
      "Mark as invalid the displacements whose path crosses more than "
      "max_unobserved_pixels_crossed distinct unobserved pixels of the mask. Updates result in "
      "place.");

  // ----------------------------------------------- Redshift-space correction

  nb::class_<otswap::RealSpaceCatalog>(m, "RealSpaceCatalog",
      "A catalogue moved to real space; row i describes input object i.")
    .def_prop_ro("n_objects", [] (CatalogHandle self) {
        return nb::cast<const otswap::RealSpaceCatalog&>(self).nObjects;
      })
    .def_prop_ro("positions", [] (CatalogHandle self) {
        const auto& c = nb::cast<const otswap::RealSpaceCatalog&>(self);
        return view<double>(self, c.positions.data(), {c.nObjects, 3});
      }, "Corrected positions, shape (n_objects, 3); NaN rows for the uncorrected.")
    .def_prop_ro("n_neighbours", [] (CatalogHandle self) {
        const auto& c = nb::cast<const otswap::RealSpaceCatalog&>(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(c.nNeighbours.data()),
                                   {c.nObjects});
      }, "Valid tracers averaged for each tracer, shape (n_objects,).")
    .def_prop_ro("n_realizations_averaged", [] (CatalogHandle self) {
        const auto& c = nb::cast<const otswap::RealSpaceCatalog&>(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(c.nRealizationsAveraged.data()),
                                   {c.nObjects});
      }, "Sum of the valid realizations of those tracers, shape (n_objects,).")
    .def_prop_ro("uncorrected", [] (CatalogHandle self) -> nb::object {
        const auto& c = nb::cast<const otswap::RealSpaceCatalog&>(self);
        if (c.uncorrected.empty())
          return nb::module_::import_("numpy").attr("empty")(0, "dtype"_a = "int64");
        return nb::cast(view<std::int64_t>(self, reinterpret_cast<const std::int64_t*>(c.uncorrected.data()),
                                           {c.uncorrected.size()}));
      }, "Indices of the tracers left without a correction, increasing.");

  m.def("line_of_sight_projection",
      [] (nb::handle displacement, nb::handle positions, nb::handle axis) {
        const std::vector<double> d = rows3(displacement, "displacement");
        if (positions.is_none() == axis.is_none())
          throw otswap::Error("give exactly one of positions, for a radial line of sight, "
                              "and axis, for a box");
        std::vector<double> out;
        if (!axis.is_none()) {
          out = otswap::lineOfSightProjection(d, to_unsigned(axis, "axis"));
        }
        else {
          const std::vector<double> p = rows3(positions, "positions");
          out = otswap::lineOfSightProjection(p, d);
        }
        const std::size_t n = out.size();
        return owned(std::move(out), {n});
      },
      "displacement"_a.none(), nb::kw_only(), "positions"_a = nb::none(), "axis"_a = nb::none(),
      "Component of each displacement along its line of sight, shape (N,).");

  m.def("neighbour_average",
      [] (nb::handle positions, nb::handle values, nb::handle validRealizations, nb::handle sigma,
          nb::handle weightByRealizations, nb::handle diagnostics) -> nb::object {
        const std::vector<double> p = rows3(positions, "positions");
        const std::vector<double> v = column(values, "values");
        const std::vector<unsigned> r = counts(validRealizations, "valid_realizations");
        const double s = to_double(sigma, "sigma");
        const bool weight = to_bool(weightByRealizations, "weight_by_realizations");
        const bool withDiagnostics = to_bool(diagnostics, "diagnostics");
        std::vector<double> average;
        std::vector<unsigned> nNeighbours, nRealizations;
        {
          nb::gil_scoped_release released;
          average = otswap::neighbourAverage(p, v, r, s, weight, nNeighbours, nRealizations);
        }
        const std::size_t n = average.size();
        nb::object a = nb::cast(owned(std::move(average), {n}));
        if (!withDiagnostics) return a;
        return nb::make_tuple(a, owned(std::move(nNeighbours), {n}), owned(std::move(nRealizations), {n}));
      },
      "positions"_a.none(), "values"_a.none(), "valid_realizations"_a.none(), nb::kw_only(),
      "sigma"_a.none(), "weight_by_realizations"_a.none() = false, "diagnostics"_a.none() = false,
      "Gaussian average of values over the valid neighbours of each object.");

  m.def("rsd_factor",
      [] (nb::handle redshift, nb::handle distances, nb::handle biasRedshift, nb::handle bias) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const std::vector<double> z = column(redshift, "redshift");
        const std::vector<double> zb = column(biasRedshift, "bias_redshift");
        const std::vector<double> b = column(bias, "bias");
        std::size_t nExtrapolated = 0;
        std::vector<double> factor;
        {
          nb::gil_scoped_release released;
          factor = otswap::rsdFactor(z, table, zb, b, nExtrapolated);
        }
        warn_extrapolation(nExtrapolated, z.size(), zb);
        const std::size_t n = factor.size();
        return owned(std::move(factor), {n});
      },
      "redshift"_a.none(), "distances"_a.none(), nb::kw_only(), "bias_redshift"_a.none(),
      "bias"_a.none(),
      "The factor f/(b + 3f/5) at each redshift; warns ExtrapolationWarning once when b(z) "
      "is extrapolated.");

  m.def("rsd_factor_box",
      [] (nb::handle redshift, nb::handle distances, nb::handle bias) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        return otswap::rsdFactorBox(to_double(redshift, "redshift"), table, to_double(bias, "bias"));
      },
      "redshift"_a.none(), "distances"_a.none(), nb::kw_only(), "bias"_a.none(),
      "The factor f/(b + 3f/5) of a box at a single redshift, with a constant bias.");

  m.def("shift_along_line_of_sight",
      [] (nb::handle positions, nb::handle shift, nb::handle axis) {
        const std::vector<double> p = rows3(positions, "positions");
        const std::vector<double> s = column(shift, "shift");
        std::vector<double> out = axis.is_none()
          ? otswap::shiftAlongLineOfSight(p, s)
          : otswap::shiftAlongLineOfSight(p, s, to_unsigned(axis, "axis"));
        const std::size_t n = out.size() / 3;
        return owned(std::move(out), {n, 3});
      },
      "positions"_a.none(), "shift"_a.none(), nb::kw_only(), "axis"_a = nb::none(),
      "Move each position by its shift along its line of sight, or along an axis.");

  m.def("real_space_lightcone",
      [] (nb::handle result, nb::handle tracersSky, nb::handle distances, nb::handle biasRedshift,
          nb::handle bias, nb::handle sigma, nb::handle angleUnit, nb::handle weightByRealizations) {
        const otswap::Result& r = instance<otswap::Result>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        const std::vector<double> given = rows3(tracersSky, "tracers_sky");
        std::vector<double> ts = given;
        if (degrees) sky_to_radians(ts, "tracers_sky");
        const std::vector<double> zb = column(biasRedshift, "bias_redshift");
        const std::vector<double> b = column(bias, "bias");
        const double s = to_double(sigma, "sigma");
        const bool weight = to_bool(weightByRealizations, "weight_by_realizations");
        std::size_t nExtrapolated = 0;
        otswap::RealSpaceCatalog catalog;
        {
          nb::gil_scoped_release released;
          catalog = otswap::realSpaceLightcone(r, ts, table, zb, b, s, weight, nExtrapolated);
        }
        // Right ascension and declination come back as given, in the unit
        // of the call, rather than converted back from radians.
        for (std::size_t i = 0; i < catalog.nObjects; ++i)
          if (!std::isnan(catalog.positions[3*i+2])) {
            catalog.positions[3*i]   = given[3*i];
            catalog.positions[3*i+1] = given[3*i+1];
          }
        std::size_t corrected = 0;
        for (std::size_t i = 0; i < r.nObjects; ++i)
          if ((r.outsideRedshiftCut.empty() || r.outsideRedshiftCut[i] == 0) &&
              (r.outsideMask.empty() || r.outsideMask[i] == 0)) ++corrected;
        warn_extrapolation(nExtrapolated, corrected, zb);
        return catalog;
      },
      "result"_a.none(), "tracers_sky"_a.none(), nb::kw_only(), "distances"_a.none(),
      "bias_redshift"_a.none(), "bias"_a.none(), "sigma"_a.none(), "angle_unit"_a.none(),
      "weight_by_realizations"_a.none() = false,
      "Move a lightcone catalogue from redshift space to real space.");

  m.def("real_space_box",
      [] (nb::handle result, nb::handle tracers, nb::handle axis, nb::handle redshift,
          nb::handle distances, nb::handle bias, nb::handle sigma, nb::handle weightByRealizations) {
        const otswap::Result& r = instance<otswap::Result>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const std::vector<double> t = rows3(tracers, "tracers");
        const unsigned a = to_unsigned(axis, "axis");
        const double z = to_double(redshift, "redshift");
        const double b = to_double(bias, "bias");
        const double s = to_double(sigma, "sigma");
        const bool weight = to_bool(weightByRealizations, "weight_by_realizations");
        nb::gil_scoped_release released;
        return otswap::realSpaceBox(r, t, a, z, table, b, s, weight);
      },
      "result"_a.none(), "tracers"_a.none(), nb::kw_only(), "axis"_a.none(), "redshift"_a.none(),
      "distances"_a.none(), "bias"_a.none(), "sigma"_a.none(),
      "weight_by_realizations"_a.none() = false,
      "Move a box catalogue from redshift space to real space.");
}
