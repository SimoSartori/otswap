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
 *  python/otswap/__init__.pyi, and of include/otswap/io.h, as the
 *  submodule otswap.io specified by python/otswap/io.pyi.
 *
 *  Every argument is taken as a Python object and converted here, so that
 *  any invalid argument raises otswap.Error; every argument is declared
 *  .none(), so that None reaches the conversions and raises otswap.Error
 *  like any other invalid value. Arrays are converted with
 *  numpy.asarray to C-contiguous float64 and their shapes checked; the C++
 *  library then works on its own copy. Right ascension and declination are
 *  converted from degrees here, so the library only ever sees radians; a
 *  Result records the angle unit of its call, and returns its sky arrays in
 *  it, converted once and kept. The other Result arrays are read-only views
 *  on the C++ result, kept alive by the Python object that holds it. The
 *  report of a call is written by the library into a string, and printed
 *  from Python, to sys.stdout, which a notebook shows, rather than by the
 *  library to the process's stderr; the extrapolation of b(z) is an
 *  ExtrapolationWarning instead of a line.
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

static_assert(sizeof(unsigned) == 4, "valid_realizations is exposed as uint32");
static_assert(sizeof(bool) == 1 && sizeof(std::uint8_t) == 1, "valid is exposed as bool");
static_assert(sizeof(std::size_t) == 8, "MpsProfile.count is exposed as uint64");
static_assert(sizeof(otswap::CorrectionStatus) == 1, "status is exposed as uint8");

namespace {

  /// pi/180 as a double: the factor numpy.deg2rad and otswap::skyToRadians
  /// multiply by, so that an angle converted here and one converted with
  /// either are the same double.
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

  /// The argument as a C-contiguous float64 numpy array of the same shape.
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

  /// An (N, 3) array, copied into the flat layout the library takes.
  std::vector<double> rows3 (nb::handle value, const std::string& name)
  {
    const CArray a = as_float64(value, name);
    if (a.ndim() != 2 || a.shape(1) != 3)
      throw otswap::Error(name + " must have shape (N, 3); it has shape " + shape_of(a));
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  /// A one-dimensional array.
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

  /// A one-dimensional array of counts: each entry an integer in
  /// [0, 4294967295].
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

  /// An instance of a bound class; the reference stays valid while the
  /// caller holds the argument.
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

  otswap::Verbosity to_verbosity (nb::handle value)
  {
    if (nb::isinstance<nb::str>(value)) {
      const std::string v = nb::str(value).c_str();
      if (v == "silent") return otswap::Verbosity::Silent;
      if (v == "normal") return otswap::Verbosity::Normal;
      if (v == "detailed") return otswap::Verbosity::Detailed;
    }
    throw otswap::Error("verbosity must be \"silent\", \"normal\" or \"detailed\"; got " +
                        repr_of(value));
  }

  /// Print each line of a report with builtins.print, leaving out those of
  /// the extrapolation of b(z), which ExtrapolationWarning carries. Needs
  /// the GIL.
  void print_report (const std::string& text)
  {
    const nb::object print = nb::module_::import_("builtins").attr("print");
    const std::string extrapolation = "otswap: b(z) extrapolated";
    std::size_t start = 0;
    for (std::size_t end; (end = text.find('\n', start)) != std::string::npos; start = end + 1)
      if (text.compare(start, extrapolation.size(), extrapolation) != 0)
        print(nb::str(text.data() + start, end - start));
  }

  /// The redshift cut from None or a (min, max) pair.
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

  /// A path: a str or an os.PathLike.
  std::string to_path (nb::handle value, const std::string& name)
  {
    nb::object path;
    try {
      path = nb::module_::import_("os").attr("fspath")(value);
    }
    catch (nb::python_error&) {
      throw otswap::Error(name + " must be a str or an os.PathLike; got " + repr_of(value));
    }
    if (!nb::isinstance<nb::str>(path))
      throw otswap::Error(name + " must be a str or an os.PathLike naming a str; got " + repr_of(value));
    return nb::str(path).c_str();
  }

  char to_char (nb::handle value, const std::string& name)
  {
    const std::string s = to_string(value, name);
    if (s.size() != 1)
      throw otswap::Error(name + " must be a single character; got " + repr_of(value));
    return s[0];
  }

  /// A sequence of column names; an int is a 0-based ASCII index.
  std::vector<std::string> column_names (nb::handle value, const std::string& name)
  {
    if (nb::isinstance<nb::str>(value))
      throw otswap::Error(name + " must be a sequence of names; got the str " + repr_of(value));
    std::vector<std::string> names;
    try {
      for (nb::handle item : nb::iter(value)) {
        if (nb::isinstance<nb::str>(item)) names.push_back(nb::str(item).c_str());
        else if (PyLong_Check(item.ptr()) && !PyBool_Check(item.ptr()))
          names.push_back(std::to_string(to_unsigned(item, name + " index")));
        else throw otswap::Error(name + " holds " + repr_of(item) + "; each entry must be a str "
                                 "or a non-negative int");
      }
    }
    catch (nb::python_error&) {
      throw otswap::Error(name + " must be a sequence of names; got " + repr_of(value));
    }
    return names;
  }

  /// A one-dimensional array of 64-bit integers: integer or boolean values
  /// taken as they are, floating ones only when each is an integer in range.
  std::vector<std::int64_t> int64_column (nb::handle value, const std::string& name)
  {
    const nb::module_ np = nb::module_::import_("numpy");
    nb::object raw;
    try {
      raw = np.attr("asarray")(value);
    }
    catch (nb::python_error& e) {
      throw otswap::Error(name + " cannot be converted to an array: " + nb::str(e.value()).c_str());
    }
    const std::string kind = nb::str(raw.attr("dtype").attr("kind")).c_str();
    if (nb::cast<std::size_t>(raw.attr("ndim")) != 1)
      throw otswap::Error(name + " must be one-dimensional; it has shape " +
                          repr_of(nb::object(raw.attr("shape"))));
    if (kind == "f") {
      const std::vector<double> c = column(raw, name);
      std::vector<std::int64_t> out(c.size());
      for (std::size_t i = 0; i < c.size(); ++i) {
        if (!(c[i] >= -9223372036854775808. && c[i] < 9223372036854775808.) || c[i] != std::floor(c[i]))
          throw otswap::Error(name + " holds " + std::to_string(c[i]) + " at entry " +
                              std::to_string(i) + "; every entry must be a 64-bit integer");
        out[i] = (std::int64_t)c[i];
      }
      return out;
    }
    if (kind != "i" && kind != "u" && kind != "b")
      throw otswap::Error(name + " must hold integers; its dtype is " + repr_of(raw.attr("dtype")));
    if (kind == "u" && nb::cast<std::size_t>(raw.attr("size")) > 0 &&
        nb::cast<unsigned long long>(raw.attr("max")()) > 9223372036854775807ULL)
      throw otswap::Error(name + " holds a value above 9223372036854775807, the largest 64-bit integer");
    using IArray = nb::ndarray<const std::int64_t, nb::c_contig, nb::device::cpu>;
    const IArray a = nb::cast<IArray>(np.attr("ascontiguousarray")(raw, "dtype"_a = np.attr("int64")));
    return std::vector<std::int64_t>(a.data(), a.data() + a.size());
  }

  /// The groups named in a sequence of str, in its order; None for the
  /// default, an empty list.
  template <typename G>
  std::vector<G> group_list (nb::handle value, const std::vector<std::pair<std::string, G>>& names)
  {
    std::vector<G> out;
    if (value.is_none()) return out;
    std::string known;
    for (const auto& n : names) known += (known.empty() ? "" : ", ") + n.first;
    if (nb::isinstance<nb::str>(value))
      throw otswap::Error("groups must be a sequence of names; got the str " + repr_of(value));
    try {
      for (nb::handle item : nb::iter(value)) {
        const std::string name = nb::isinstance<nb::str>(item) ? nb::str(item).c_str() : "";
        bool found = false;
        for (const auto& n : names)
          if (n.first == name) { out.push_back(n.second); found = true; }
        if (!found)
          throw otswap::Error("unknown group " + repr_of(item) + "; the groups are " + known);
      }
    }
    catch (nb::python_error&) {
      throw otswap::Error("groups must be a sequence of names; got " + repr_of(value));
    }
    if (out.empty())
      throw otswap::Error("groups is empty; give None for the default groups");
    return out;
  }

  // --------------------------------------------------------- returned arrays

  std::vector<std::size_t> shape_vector (const CArray& a)
  {
    std::vector<std::size_t> s(a.ndim());
    for (std::size_t i = 0; i < a.ndim(); ++i) s[i] = a.shape(i);
    return s;
  }

  /// A new numpy array that owns values.
  template <typename T>
  nb::ndarray<nb::numpy, T> owned (std::vector<T>&& values, const std::vector<std::size_t>& shape)
  {
    auto* held = new std::vector<T>(std::move(values));
    nb::capsule owner(held, [] (void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T>(held->data(), shape.size(), shape.data(), owner);
  }

  /// A DistanceTable lookup applied to every element, keeping the shape.
  template <typename F>
  nb::ndarray<nb::numpy, double> elementwise (nb::handle value, const std::string& name, F lookup)
  {
    const CArray a = as_float64(value, name);
    std::vector<double> out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) out[i] = lookup(a.data()[i]);
    return owned(std::move(out), shape_vector(a));
  }

  /// A read-only view on an array of a Result, which keeps the Python object
  /// holding the Result alive.
  template <typename T>
  nb::ndarray<nb::numpy, const T> view (nb::handle owner, const T* data,
                                        std::initializer_list<std::size_t> shape)
  {
    return nb::ndarray<nb::numpy, const T>(data, shape, owner);
  }

  // ------------------------------------------------------------- warnings

  /// otswap.ExtrapolationWarning, created with the module.
  PyObject* extrapolation_warning = nullptr;

  /// Raise one ExtrapolationWarning for a call that extrapolated b(z) at
  /// count of total redshifts. Needs the GIL.
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

  /// Check a bias table with the library's own rules and messages: the
  /// factor over no redshift checks the table and computes nothing.
  void check_bias_table (const otswap::BiasTable& table)
  {
    static const otswap::DistanceTable withGrowth({0., 1.}, {0., 1.}, {1., 1.});
    otswap::rsdFactor({}, withGrowth, table);
  }

  /// The Result of the Python package: the C++ result, the angle unit of
  /// the call that made it, and its sky arrays in that unit, converted on
  /// first access and kept until a filter changes them.
  struct PyResult : otswap::Result {
    bool degrees = false;
    nb::object tracersSkyInUnit, lagrangianSkyInUnit;

    PyResult () = default;
    PyResult (otswap::Result&& r, const bool inDegrees) : otswap::Result(std::move(r)), degrees(inDegrees) {}
  };

  using ResultHandle = nb::handle_t<PyResult>;
  /// The RealSpaceCatalog of the Python package: the C++ catalogue, the
  /// angle unit of the result it corrects, and its sky array in that unit,
  /// converted on first access and kept.
  struct PyCatalog : otswap::RealSpaceCatalog {
    bool degrees = false;
    nb::object skyInUnit;

    PyCatalog () = default;
    PyCatalog (otswap::RealSpaceCatalog&& c, const bool inDegrees)
      : otswap::RealSpaceCatalog(std::move(c)), degrees(inDegrees) {}
  };

  using CatalogHandle = nb::handle_t<PyCatalog>;

  const PyResult& result_of (ResultHandle self)
  {
    return nb::cast<const PyResult&>(self);
  }

  /// A sky array of a result in the unit of its call: a view in radians, a
  /// read-only copy in degrees, made once and kept in cache.
  nb::object sky_in_unit (nb::handle owner, const std::vector<double>& sky, const bool degrees,
                          nb::object& cache)
  {
    if (sky.empty()) return nb::none();
    const std::size_t n = sky.size() / 3;
    if (!degrees) return nb::cast(view<double>(owner, sky.data(), {n, 3}));
    if (!cache.is_valid()) {
      nb::object converted = nb::cast(owned(otswap::skyToDegrees(sky), {n, 3}));
      converted.attr("flags").attr("writeable") = false;
      cache = converted;
    }
    return cache;
  }

  /// A float64 array of the given shape, copied into the flat layout.
  std::vector<double> shaped (nb::handle value, const std::string& name,
                              const std::vector<std::size_t>& shape)
  {
    const CArray a = as_float64(value, name);
    bool same = a.ndim() == shape.size();
    for (std::size_t i = 0; same && i < shape.size(); ++i) same = a.shape(i) == shape[i];
    if (!same) {
      std::string expected = "(";
      for (std::size_t i = 0; i < shape.size(); ++i)
        expected += (i ? ", " : "") + std::to_string(shape[i]);
      throw otswap::Error(name + " must have shape " + expected + (shape.size() == 1 ? ",)" : ")") +
                          "; it has shape " + shape_of(a));
    }
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  /// 0 or 1 flags of the given shape, from bool or numbers.
  std::vector<std::uint8_t> flags (nb::handle value, const std::string& name,
                                   const std::vector<std::size_t>& shape)
  {
    const std::vector<double> v = shaped(value, name, shape);
    std::vector<std::uint8_t> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (v[i] != 0. && v[i] != 1.)
        throw otswap::Error(name + " holds " + std::to_string(v[i]) + " at entry " + std::to_string(i) +
                            "; every entry must be 0 or 1");
      out[i] = (std::uint8_t)v[i];
    }
    return out;
  }

}

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
        if (!c.redshiftCut) return nb::none();
        return nb::make_tuple(c.redshiftCut->min, c.redshiftCut->max);
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
        if (!c.maxUnobservedPixelsCrossed) return nb::none();
        return nb::int_(*c.maxUnobservedPixelsCrossed);
      }, "The threshold of the latest mask filter applied to the result, or None when none was.")
    .def_ro("displacements", &otswap::SelectionCounts::displacements)
    .def_ro("displacements_crossing_mask", &otswap::SelectionCounts::displacementsCrossingMask)
    .def("__repr__", [] (const otswap::SelectionCounts& c) {
        std::string text = c.message();
        if (text.empty()) return nb::str("SelectionCounts(no selection)");
        text.pop_back();
        return nb::str(text.c_str());
      });

  // -------------------------------------------------------------- MpsProfile

  nb::class_<otswap::MpsProfile>(m, "MpsProfile",
      "The mean particle separation of a lightcone as a function of redshift. Read-only.")
    .def_prop_ro("redshift", [] (nb::handle_t<otswap::MpsProfile> self) {
        const auto& p = nb::cast<const otswap::MpsProfile&>(self);
        return view<double>(self, p.redshift.data(), {p.redshift.size()});
      }, "Bin centres, shape (n_bins,).")
    .def_prop_ro("mps", [] (nb::handle_t<otswap::MpsProfile> self) {
        const auto& p = nb::cast<const otswap::MpsProfile&>(self);
        return view<double>(self, p.mps.data(), {p.mps.size()});
      }, "Mean particle separation at each centre, Mpc/h, shape (n_bins,).")
    .def_prop_ro("count", [] (nb::handle_t<otswap::MpsProfile> self) {
        const auto& p = nb::cast<const otswap::MpsProfile&>(self);
        return view<std::uint64_t>(self, reinterpret_cast<const std::uint64_t*>(p.count.data()),
                                   {p.count.size()});
      }, "Tracers in each bin, shape (n_bins,).")
    .def_ro("redshift_min", &otswap::MpsProfile::redshiftMin)
    .def_ro("redshift_max", &otswap::MpsProfile::redshiftMax)
    .def_ro("representative", &otswap::MpsProfile::representative)
    .def("at", [] (const otswap::MpsProfile& p, nb::handle z) {
          return elementwise(z, "z", [&p] (const double v) { return p.at(v); });
        }, "z"_a.none(), "The mps at each redshift, linear between the nodes and beyond them.");

  // ------------------------------------------------------------------ Result

  nb::class_<PyResult>(m, "Result", "Displacement field produced by a reconstruction.")
    .def_static("from_arrays",
        [] (nb::handle displacement, nb::handle matchedRandom, nb::handle valid, nb::handle tracers,
            nb::handle tracersSky, nb::handle angleUnit, nb::handle distances,
            nb::handle outsideRedshiftCut, nb::handle outsideMask, nb::handle seed) {
          const CArray d = as_float64(displacement, "displacement");
          if (d.ndim() != 3 || d.shape(2) != 3 || d.shape(0) == 0 || d.shape(1) == 0)
            throw otswap::Error("displacement must have shape (n_realizations, n_objects, 3), neither "
                                "0; it has shape " + shape_of(d));
          const std::size_t nRec = d.shape(0), n = d.shape(1);
          otswap::Result r;
          r.nObjects = n;
          r.nRealizations = (unsigned)nRec;
          r.displacement.assign(d.data(), d.data() + d.size());
          r.matchedRandom = shaped(matchedRandom, "matched_random", {nRec, n, 3});
          r.valid = flags(valid, "valid", {nRec, n});
          r.tracers = shaped(tracers, "tracers", {n, 3});
          r.outsideRedshiftCut = outsideRedshiftCut.is_none()
            ? std::vector<std::uint8_t>(n, 0) : flags(outsideRedshiftCut, "outside_redshift_cut", {n});
          r.outsideMask = outsideMask.is_none()
            ? std::vector<std::uint8_t>(n, 0) : flags(outsideMask, "outside_mask", {n});
          r.config.nRealizations = (unsigned)nRec;
          r.config.seed = to_unsigned(seed, "seed");
          bool degrees = false;
          if (!tracersSky.is_none()) {
            degrees = in_degrees(angleUnit);
            r.geometry = otswap::Geometry::Lightcone;
            r.tracersSky = shaped(tracersSky, "tracers_sky", {n, 3});
            if (degrees) sky_to_radians(r.tracersSky, "tracers_sky");
            r.distances = std::make_shared<const otswap::DistanceTable>(
              instance<otswap::DistanceTable>(distances, "distances", "DistanceTable"));
          }
          else if (!angleUnit.is_none() || !distances.is_none())
            throw otswap::Error("angle_unit and distances describe tracers_sky; give them with it, "
                                "for a lightcone result, or not at all");
          otswap::recomputeMeans(r);
          return PyResult(std::move(r), degrees);
        },
        "displacement"_a.none(), "matched_random"_a.none(), "valid"_a.none(), "tracers"_a.none(),
        nb::kw_only(), "tracers_sky"_a = nb::none(), "angle_unit"_a = nb::none(),
        "distances"_a = nb::none(), "outside_redshift_cut"_a = nb::none(),
        "outside_mask"_a = nb::none(), "seed"_a.none() = 0,
        "Build a result from its arrays, as read back from a file; the means are recomputed.")
    .def_prop_ro("n_objects", [] (ResultHandle self) { return result_of(self).nObjects; })
    .def_prop_ro("n_realizations", [] (ResultHandle self) { return result_of(self).nRealizations; })
    .def_prop_ro("geometry", [] (ResultHandle self) {
        return nb::str(result_of(self).geometry == otswap::Geometry::Box ? "box" : "lightcone");
      }, "\"box\" or \"lightcone\".")
    .def_prop_ro("angle_unit", [] (ResultHandle self) -> nb::object {
        const PyResult& r = result_of(self);
        if (r.geometry == otswap::Geometry::Box) return nb::none();
        return nb::str(r.degrees ? "deg" : "rad");
      }, "The angle_unit of the sky arrays, that of the call; None for a box.")
    .def_prop_ro("seed", [] (ResultHandle self) { return result_of(self).config.seed; },
                 "The seed used: the drawn one when 0 was given.")
    .def_prop_ro("displacement", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<double>(self, r.displacement.data(), {r.nRealizations, r.nObjects, 3});
      }, "Displacements, shape (n_realizations, n_objects, 3).")
    .def_prop_ro("matched_random", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<double>(self, r.matchedRandom.data(), {r.nRealizations, r.nObjects, 3});
      }, "Cartesian position of the random matched to each tracer, shape "
         "(n_realizations, n_objects, 3).")
    .def_prop_ro("valid", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.valid.data()),
                          {r.nRealizations, r.nObjects});
      }, "Whether each displacement is valid, shape (n_realizations, n_objects).")
    .def_prop_ro("valid_realizations", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(r.validRealizations.data()),
                                   {r.nObjects});
      }, "Number of valid realizations per tracer, shape (n_objects,).")
    .def_prop_ro("mean_displacement", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<double>(self, r.meanDisplacement.data(), {r.nObjects, 3});
      }, "Displacement averaged over the valid realizations, shape (n_objects, 3).")
    .def_prop_ro("tracers", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<double>(self, r.tracers.data(), {r.nObjects, 3});
      }, "Cartesian position of each tracer, shape (n_objects, 3).")
    .def_prop_ro("tracers_sky", [] (ResultHandle self) {
        PyResult& r = nb::cast<PyResult&>(self);
        return sky_in_unit(self, r.tracersSky, r.degrees, r.tracersSkyInUnit);
      }, "Sky coordinates of each tracer in angle_unit, shape (n_objects, 3); None for a box.")
    .def_prop_ro("lagrangian", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<double>(self, r.lagrangian.data(), {r.nObjects, 3});
      }, "Mean Lagrangian position, tracers + mean_displacement, shape (n_objects, 3).")
    .def_prop_ro("filtered_nside", [] (ResultHandle self) { return result_of(self).filteredNside; },
                 "NSIDE of the mask the result was filtered against, 0 if none.")
    .def_prop_ro("outside_redshift_cut", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.outsideRedshiftCut.data()),
                          {r.nObjects});
      }, "Whether each tracer lay outside the redshift cut, shape (n_objects,).")
    .def_prop_ro("outside_mask", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.outsideMask.data()), {r.nObjects});
      }, "Whether each tracer fell on an unobserved pixel of the mask, shape (n_objects,).")
    .def_prop_ro("selection", [] (const PyResult& r) -> const otswap::SelectionCounts& {
        return r.selection;
      }, nb::rv_policy::reference_internal, "What the reconstruction left out, and why.")
    .def_prop_ro("lagrangian_sky", [] (ResultHandle self) {
        PyResult& r = nb::cast<PyResult&>(self);
        return sky_in_unit(self, r.lagrangianSky, r.degrees, r.lagrangianSkyInUnit);
      }, "Sky coordinates of each tracer's mean Lagrangian position in angle_unit, shape "
         "(n_objects, 3); None for a box.")
    .def_prop_ro("mps", [] (ResultHandle self) -> nb::object {
        const PyResult& r = result_of(self);
        if (std::isnan(r.mps)) return nb::none();
        return nb::float_(r.mps);
      }, "Box: the mean particle separation used, Mpc/h; None for a lightcone.")
    .def_prop_ro("mps_profile", [] (const PyResult& r) -> const otswap::MpsProfile* {
        return r.mpsProfile.redshift.empty() ? nullptr : &r.mpsProfile;
      }, nb::rv_policy::reference_internal, "Lightcone: the mps(z) profile; None for a box.")
    .def_prop_ro("elapsed_seconds", [] (ResultHandle self) { return result_of(self).elapsedSeconds; },
                 "Wall time of the call that made the result, in seconds.");

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
        if (degrees) sky = otswap::skyToDegrees(std::move(sky));
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
          nb::handle convergence, nb::handle seed, nb::handle cellSize, nb::handle verbosity) {
        const std::vector<double> t = rows3(tracers, "tracers");
        std::vector<double> r;
        if (!randoms.is_none()) {
          r = rows3(randoms, "randoms");
          if (r.empty())
            throw otswap::Error("randoms is empty; omit it to have the randoms drawn "
                                "in the bounding box of the tracers");
        }
        const bool computed = mps.is_none();
        const double separation = computed ? 0. : to_double(mps, "mps");
        otswap::Config config = make_config(nRealizations, convergence, seed, cellSize);
        config.verbosity = to_verbosity(verbosity);
        std::ostringstream report;
        config.log = &report;
        otswap::Result result;
        {
          nb::gil_scoped_release released;
          result = computed ? otswap::reconstructBox(t, r, config)
                            : otswap::reconstructBox(t, r, separation, config);
        }
        print_report(report.str());
        return PyResult(std::move(result), false);
      },
      "tracers"_a.none(), "randoms"_a = nb::none(), nb::kw_only(), "mps"_a = nb::none(),
      "n_realizations"_a.none() = 1, "convergence"_a.none() = 1.e-3, "seed"_a.none() = 0,
      "cell_size"_a.none() = 4., "verbosity"_a.none() = "normal",
      "Reconstruct in box geometry, with a constant mean particle separation, given or that of "
      "the tracers' bounding box.");

  m.def("reconstruct_lightcone",
      [] (nb::handle tracersSky, nb::handle randomsSky, nb::handle skyAreaDeg2, nb::handle mask,
          nb::handle nBins, nb::handle distances, nb::handle angleUnit, nb::handle tracers,
          nb::handle randoms, nb::handle nRealizations, nb::handle convergence, nb::handle seed,
          nb::handle cellSize, nb::handle redshiftCut, nb::handle rejectCrossings,
          nb::handle maxUnobservedPixelsCrossed, nb::handle verbosity) {
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
        config.verbosity = to_verbosity(verbosity);
        std::ostringstream report;
        config.log = &report;
        const otswap::RedshiftCut cut = to_cut(redshiftCut);
        std::vector<double> t, r;
        if (!tracers.is_none()) {
          t = rows3(tracers, "tracers");
          r = rows3(randoms, "randoms");
        }

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
        print_report(report.str());
        return PyResult(std::move(result), degrees);
      },
      "tracers_sky"_a.none(), "randoms_sky"_a.none(), nb::kw_only(), "sky_area_deg2"_a = nb::none(),
      "mask"_a = nb::none(), "n_bins"_a.none(), "distances"_a.none(), "angle_unit"_a.none(),
      "tracers"_a = nb::none(), "randoms"_a = nb::none(), "n_realizations"_a.none() = 1,
      "convergence"_a.none() = 1.e-3, "seed"_a.none() = 0, "cell_size"_a.none() = 4.,
      "redshift_cut"_a = nb::none(), "reject_crossings"_a.none() = true,
      "max_unobserved_pixels_crossed"_a.none() = 0, "verbosity"_a.none() = "normal",
      "Reconstruct in lightcone geometry.");

  m.def("reject_mask_crossings",
      [] (nb::handle result, nb::handle mask, nb::handle maxUnobservedPixelsCrossed) {
        PyResult& r = instance<PyResult>(result, "result", "Result");
        const otswap::Mask& k = instance<otswap::Mask>(mask, "mask", "Mask");
        const unsigned limit = to_unsigned(maxUnobservedPixelsCrossed, "max_unobserved_pixels_crossed");
        {
          nb::gil_scoped_release released;
          otswap::rejectMaskCrossings(r, k, limit);
        }
        r.lagrangianSkyInUnit = nb::object();
      },
      "result"_a.none(), "mask"_a.none(), "max_unobserved_pixels_crossed"_a.none() = 0,
      "Mark as invalid the displacements whose path crosses more than "
      "max_unobserved_pixels_crossed distinct unobserved pixels of the mask. Updates result in "
      "place.");

  m.def("recompute_means",
      [] (nb::handle result) {
        PyResult& r = instance<PyResult>(result, "result", "Result");
        otswap::recomputeMeans(r);
        r.lagrangianSkyInUnit = nb::object();
      },
      "result"_a.none(),
      "Recompute valid_realizations, mean_displacement, lagrangian and lagrangian_sky from "
      "displacement and valid, as the reconstruction computes them. Updates result in place.");

  // --------------------------------------------------------------- BiasTable

  nb::class_<otswap::BiasTable>(m, "BiasTable",
      "A b(z) table: the linear bias at redshift nodes, interpolated linearly between them and "
      "extrapolated linearly beyond them. Immutable.")
    .def("__init__", [] (otswap::BiasTable* self, nb::handle redshift, nb::handle bias) {
          otswap::BiasTable table;
          table.redshift = column(redshift, "redshift");
          table.bias = column(bias, "bias");
          check_bias_table(table);
          new (self) otswap::BiasTable(std::move(table));
        }, "redshift"_a.none(), "bias"_a.none(),
        "At least two nodes; redshifts finite and strictly increasing; bias finite and positive.")
    .def_prop_ro("redshift", [] (nb::handle_t<otswap::BiasTable> self) {
        const auto& t = nb::cast<const otswap::BiasTable&>(self);
        return view<double>(self, t.redshift.data(), {t.redshift.size()});
      }, "Redshifts of the nodes, shape (n,).")
    .def_prop_ro("bias", [] (nb::handle_t<otswap::BiasTable> self) {
        const auto& t = nb::cast<const otswap::BiasTable&>(self);
        return view<double>(self, t.bias.data(), {t.bias.size()});
      }, "Bias at each node, shape (n,).");

  // ----------------------------------------------- Redshift-space correction

  nb::class_<PyCatalog>(m, "RealSpaceCatalog",
      "A catalogue moved to real space; row i describes input tracer i.")
    .def_prop_ro("n_objects", [] (const PyCatalog& c) { return c.nObjects; })
    .def_prop_ro("geometry", [] (const PyCatalog& c) {
        return nb::str(c.geometry == otswap::Geometry::Box ? "box" : "lightcone");
      }, "\"box\" or \"lightcone\".")
    .def_prop_ro("angle_unit", [] (const PyCatalog& c) -> nb::object {
        if (c.geometry == otswap::Geometry::Box) return nb::none();
        return nb::str(c.degrees ? "deg" : "rad");
      }, "The angle_unit of sky, that of the result corrected; None for a box.")
    .def_prop_ro("sky", [] (CatalogHandle self) {
        PyCatalog& c = nb::cast<PyCatalog&>(self);
        return sky_in_unit(self, c.sky, c.degrees, c.skyInUnit);
      }, "Lightcone: right ascension and declination of each tracer and its corrected redshift, "
         "in angle_unit, shape (n_objects, 3); None for a box.")
    .def_prop_ro("cartesian", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<double>(self, c.cartesian.data(), {c.nObjects, 3});
      }, "Corrected Cartesian positions, Mpc/h, shape (n_objects, 3); NaN rows for the tracers "
         "not moved.")
    .def_prop_ro("shift", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<double>(self, c.shift.data(), {c.nObjects});
      }, "Shift applied along the line of sight, Mpc/h, shape (n_objects,); NaN where none was.")
    .def_prop_ro("factor", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<double>(self, c.factor.data(), {c.nObjects});
      }, "The factor f/(b + 3f/5) at each tracer, shape (n_objects,); NaN for those left out.")
    .def_prop_ro("status", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<std::uint8_t>(self, reinterpret_cast<const std::uint8_t*>(c.status.data()),
                                  {c.nObjects});
      }, "What the correction did with each tracer, shape (n_objects,): 0 corrected, 1 moved by "
         "its neighbours, 2 no valid neighbour, 3 left out.")
    .def_prop_ro("valid_realizations", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(c.validRealizations.data()),
                                   {c.nObjects});
      }, "Valid realizations of each tracer, as in the result, shape (n_objects,).")
    .def_prop_ro("n_neighbours", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(c.nNeighbours.data()),
                                   {c.nObjects});
      }, "Valid tracers averaged for each tracer, shape (n_objects,).")
    .def_prop_ro("n_realizations_averaged", [] (CatalogHandle self) {
        const auto& c = nb::cast<const PyCatalog&>(self);
        return view<std::uint32_t>(self, reinterpret_cast<const std::uint32_t*>(c.nRealizationsAveraged.data()),
                                   {c.nObjects});
      }, "Sum of the valid realizations of those tracers, shape (n_objects,).")
    .def_prop_ro("n_extrapolated", [] (const PyCatalog& c) { return c.nExtrapolated; },
                 "Tracers at whose redshift b(z) was extrapolated; 0 for a box.")
    .def_prop_ro("sigma", [] (const PyCatalog& c) { return c.sigma; }, "The width of the average, Mpc/h.")
    .def_prop_ro("weight_by_realizations", [] (const PyCatalog& c) { return c.weightByRealizations; })
    .def_prop_ro("axis", [] (const PyCatalog& c) -> nb::object {
        if (c.geometry != otswap::Geometry::Box) return nb::none();
        return nb::int_(c.axis);
      }, "Box: the line-of-sight axis; None for a lightcone.")
    .def_prop_ro("box_redshift", [] (const PyCatalog& c) -> nb::object {
        if (c.geometry != otswap::Geometry::Box) return nb::none();
        return nb::float_(c.boxRedshift);
      }, "Box: the redshift of the box; None for a lightcone.")
    .def_prop_ro("box_bias", [] (const PyCatalog& c) -> nb::object {
        if (c.geometry != otswap::Geometry::Box) return nb::none();
        return nb::float_(c.boxBias);
      }, "Box: the bias of the tracers; None for a lightcone.")
    .def_prop_ro("elapsed_seconds", [] (const PyCatalog& c) { return c.elapsedSeconds; },
                 "Wall time of the call that made the catalogue, in seconds.");

  m.def("radial_projection",
      [] (nb::handle displacement, nb::handle positions) {
        const std::vector<double> d = rows3(displacement, "displacement");
        const std::vector<double> p = rows3(positions, "positions");
        std::vector<double> out = otswap::radialProjection(d, p);
        const std::size_t n = out.size();
        return owned(std::move(out), {n});
      },
      "displacement"_a.none(), "positions"_a.none(),
      "Component of each displacement along the line of sight of its position, shape (N,).");

  m.def("axis_projection",
      [] (nb::handle displacement, nb::handle axis) {
        const std::vector<double> d = rows3(displacement, "displacement");
        std::vector<double> out = otswap::axisProjection(d, to_unsigned(axis, "axis"));
        const std::size_t n = out.size();
        return owned(std::move(out), {n});
      },
      "displacement"_a.none(), "axis"_a.none(),
      "Component of each displacement along a Cartesian axis, shape (N,).");

  m.def("neighbour_average",
      [] (nb::handle positions, nb::handle values, nb::handle validRealizations, nb::handle sigma,
          nb::handle weightByRealizations) {
        const std::vector<double> p = rows3(positions, "positions");
        const std::vector<double> v = column(values, "values");
        const std::vector<unsigned> r = counts(validRealizations, "valid_realizations");
        const double s = to_double(sigma, "sigma");
        const bool weight = to_bool(weightByRealizations, "weight_by_realizations");
        otswap::NeighbourAverage average;
        {
          nb::gil_scoped_release released;
          average = otswap::neighbourAverage(p, v, r, s, weight);
        }
        const std::size_t n = average.values.size();
        return nb::module_::import_("otswap").attr("NeighbourAverage")(
          owned(std::move(average.values), {n}), owned(std::move(average.nNeighbours), {n}),
          owned(std::move(average.nRealizationsAveraged), {n}));
      },
      "positions"_a.none(), "values"_a.none(), "valid_realizations"_a.none(), nb::kw_only(),
      "sigma"_a.none(), "weight_by_realizations"_a.none() = false,
      "Gaussian average of values over the valid neighbours of each object.");

  m.def("rsd_factor",
      [] (nb::handle redshift, nb::handle distances, nb::handle bias) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const otswap::BiasTable& b = instance<otswap::BiasTable>(bias, "bias", "BiasTable");
        const std::vector<double> z = column(redshift, "redshift");
        std::size_t nExtrapolated = 0;
        std::vector<double> factor;
        {
          nb::gil_scoped_release released;
          factor = otswap::rsdFactor(z, table, b, &nExtrapolated);
        }
        warn_extrapolation(nExtrapolated, z.size(), b.redshift);
        const std::size_t n = factor.size();
        return owned(std::move(factor), {n});
      },
      "redshift"_a.none(), "distances"_a.none(), nb::kw_only(), "bias"_a.none(),
      "The factor f/(b + 3f/5) at each redshift; warns ExtrapolationWarning once when b(z) "
      "is extrapolated.");

  m.def("rsd_factor_box",
      [] (nb::handle redshift, nb::handle distances, nb::handle bias) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        return otswap::rsdFactorBox(to_double(redshift, "redshift"), table, to_double(bias, "bias"));
      },
      "redshift"_a.none(), "distances"_a.none(), nb::kw_only(), "bias"_a.none(),
      "The factor f/(b + 3f/5) of a box at a single redshift, with a constant bias.");

  m.def("shift_radially",
      [] (nb::handle positions, nb::handle shift) {
        const std::vector<double> p = rows3(positions, "positions");
        const std::vector<double> s = column(shift, "shift");
        std::vector<double> out = otswap::shiftRadially(p, s);
        const std::size_t n = out.size() / 3;
        return owned(std::move(out), {n, 3});
      },
      "positions"_a.none(), "shift"_a.none(),
      "Move each position by its shift along its own line of sight.");

  m.def("shift_along_axis",
      [] (nb::handle positions, nb::handle shift, nb::handle axis) {
        const std::vector<double> p = rows3(positions, "positions");
        const std::vector<double> s = column(shift, "shift");
        std::vector<double> out = otswap::shiftAlongAxis(p, s, to_unsigned(axis, "axis"));
        const std::size_t n = out.size() / 3;
        return owned(std::move(out), {n, 3});
      },
      "positions"_a.none(), "shift"_a.none(), "axis"_a.none(),
      "Move each position by its shift along a Cartesian axis, without periodic wrapping.");

  m.def("real_space_lightcone",
      [] (nb::handle result, nb::handle distances, nb::handle bias, nb::handle sigma,
          nb::handle weightByRealizations, nb::handle verbosity) {
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const otswap::BiasTable& b = instance<otswap::BiasTable>(bias, "bias", "BiasTable");
        const double s = to_double(sigma, "sigma");
        otswap::CorrectionConfig config;
        config.weightByRealizations = to_bool(weightByRealizations, "weight_by_realizations");
        config.verbosity = to_verbosity(verbosity);
        std::ostringstream report;
        config.log = &report;
        otswap::RealSpaceCatalog catalog;
        {
          nb::gil_scoped_release released;
          catalog = otswap::realSpaceLightcone(r, table, b, s, config);
        }
        print_report(report.str());
        std::size_t kept = 0;
        for (const otswap::CorrectionStatus status : catalog.status)
          kept += status != otswap::CorrectionStatus::LeftOut;
        warn_extrapolation(catalog.nExtrapolated, kept, b.redshift);
        return PyCatalog(std::move(catalog), r.degrees);
      },
      "result"_a.none(), nb::kw_only(), "distances"_a.none(), "bias"_a.none(), "sigma"_a.none(),
      "weight_by_realizations"_a.none() = false, "verbosity"_a.none() = "normal",
      "Move a lightcone catalogue from redshift space to real space.");

  m.def("real_space_box",
      [] (nb::handle result, nb::handle axis, nb::handle redshift, nb::handle distances,
          nb::handle bias, nb::handle sigma, nb::handle weightByRealizations, nb::handle verbosity) {
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const unsigned a = to_unsigned(axis, "axis");
        const double z = to_double(redshift, "redshift");
        const double b = to_double(bias, "bias");
        const double s = to_double(sigma, "sigma");
        otswap::CorrectionConfig config;
        config.weightByRealizations = to_bool(weightByRealizations, "weight_by_realizations");
        config.verbosity = to_verbosity(verbosity);
        std::ostringstream report;
        config.log = &report;
        otswap::RealSpaceCatalog catalog;
        {
          nb::gil_scoped_release released;
          catalog = otswap::realSpaceBox(r, a, z, table, b, s, config);
        }
        print_report(report.str());
        return PyCatalog(std::move(catalog), false);
      },
      "result"_a.none(), nb::kw_only(), "axis"_a.none(), "redshift"_a.none(), "distances"_a.none(),
      "bias"_a.none(), "sigma"_a.none(), "weight_by_realizations"_a.none() = false,
      "verbosity"_a.none() = "normal",
      "Move a box catalogue from redshift space to real space.");

  // ---------------------------------------------------------------------- io

  nb::module_ io = m.def_submodule("io", "Reading and writing ASCII and FITS tables.");

  nb::class_<otswap::io::Column>(io, "Column", "One column of a table to write. Read-only.")
    .def("__init__",
        [] (otswap::io::Column* self, nb::handle name, nb::handle data, nb::handle type,
            nb::handle description, nb::handle unit) {
          otswap::io::Column c;
          c.name = to_string(name, "name");
          const std::string t = to_string(type, "type");
          if (t != "D" && t != "J" && t != "K")
            throw otswap::Error("type must be \"D\", \"J\" or \"K\"; got " + repr_of(type));
          c.type = t[0];
          c.description = to_string(description, "description");
          c.unit = to_string(unit, "unit");
          if (c.type == 'K') c.integers = int64_column(data, "data");
          else               c.data = column(data, "data");
          new (self) otswap::io::Column(std::move(c));
        },
        "name"_a.none(), "data"_a.none(), nb::kw_only(), "type"_a.none() = "D",
        "description"_a.none() = "", "unit"_a.none() = "")
    .def_prop_ro("name", [] (const otswap::io::Column& c) { return nb::str(c.name.c_str()); })
    .def_prop_ro("type", [] (const otswap::io::Column& c) { return nb::str(std::string(1, c.type).c_str()); })
    .def_prop_ro("description", [] (const otswap::io::Column& c) { return nb::str(c.description.c_str()); })
    .def_prop_ro("unit", [] (const otswap::io::Column& c) { return nb::str(c.unit.c_str()); })
    .def_prop_ro("data", [] (nb::handle_t<otswap::io::Column> self) -> nb::object {
        const auto& c = nb::cast<const otswap::io::Column&>(self);
        if (c.type == 'K') return nb::cast(view<std::int64_t>(self, c.integers.data(), {c.integers.size()}));
        return nb::cast(view<double>(self, c.data.data(), {c.data.size()}));
      }, "The values: float64, or int64 for a 'K' column.");

  nb::class_<otswap::io::Table>(io, "Table", "A table read from file. Read-only.")
    .def_prop_ro("n_rows", [] (const otswap::io::Table& t) { return t.nRows; })
    .def_prop_ro("n_columns", [] (const otswap::io::Table& t) { return t.nColumns; })
    .def_prop_ro("n_integer_columns", [] (const otswap::io::Table& t) { return t.nIntegerColumns; })
    .def_prop_ro("values", [] (nb::handle_t<otswap::io::Table> self) {
        const auto& t = nb::cast<const otswap::io::Table&>(self);
        return view<double>(self, t.values.data(), {t.nRows, t.nColumns});
      }, "The columns read as float64, shape (n_rows, n_columns).")
    .def_prop_ro("integers", [] (nb::handle_t<otswap::io::Table> self) {
        const auto& t = nb::cast<const otswap::io::Table&>(self);
        return view<std::int64_t>(self, t.integers.data(), {t.nRows, t.nIntegerColumns});
      }, "The columns read as int64, shape (n_rows, n_integer_columns).");

  io.def("read",
      [] (nb::handle file, nb::handle columns, nb::handle integerColumns, nb::handle delimiter,
          nb::handle comment) {
        const std::string path = to_path(file, "file");
        const std::vector<std::string> c = column_names(columns, "columns");
        const std::vector<std::string> i = column_names(integerColumns, "integer_columns");
        const char d = to_char(delimiter, "delimiter"), k = to_char(comment, "comment");
        nb::gil_scoped_release released;
        return otswap::io::read(path, c, i, d, k);
      },
      "file"_a.none(), "columns"_a.none(), nb::kw_only(), "integer_columns"_a.none() = nb::tuple(),
      "delimiter"_a.none() = " ", "comment"_a.none() = "#",
      "Read selected columns of an ASCII or FITS table.");

  io.def("write",
      [] (nb::handle file, nb::handle columns, nb::handle precision, nb::handle keywords) {
        const std::string path = to_path(file, "file");
        std::vector<const otswap::io::Column*> given;
        std::vector<otswap::io::Column> header;
        try {
          for (nb::handle item : nb::iter(columns)) {
            const otswap::io::Column& c = instance<otswap::io::Column>(item, "each column", "io.Column");
            given.push_back(&c);
            header.push_back({c.name, c.type, c.description, c.unit, {}, {}});
          }
        }
        catch (nb::python_error&) {
          throw otswap::Error("columns must be a sequence of otswap.io.Column; got " + repr_of(columns));
        }
        otswap::io::WriteOptions options;
        const unsigned p = to_unsigned(precision, "precision");
        options.precision = p > 1000u ? 1000 : (int)p;
        if (!keywords.is_none()) {
          try {
            for (nb::handle item : nb::iter(keywords)) {
              if (!nb::isinstance<nb::tuple>(item) || nb::len(item) != 3)
                throw otswap::Error("each keyword must be a (name, value, comment) tuple; got " +
                                    repr_of(item));
              const nb::tuple k = nb::borrow<nb::tuple>(item);
              options.keywords.push_back({to_string(k[0], "a keyword name"),
                                          to_string(k[1], "a keyword value"),
                                          to_string(k[2], "a keyword comment")});
            }
          }
          catch (nb::python_error&) {
            throw otswap::Error("keywords must be a sequence of (name, value, comment) tuples; got " +
                                repr_of(keywords));
          }
        }
        if (given.empty()) throw otswap::Error("no column was given for: " + path);
        auto length = [] (const otswap::io::Column* c) {
          return c->type == 'K' ? c->integers.size() : c->data.size();
        };
        for (const auto* c : given)
          if (length(c) != length(given[0]))
            throw otswap::Error("the column '" + c->name + "' holds " + std::to_string(length(c)) +
                                " entries, against " + std::to_string(length(given[0])) + " in '" +
                                given[0]->name + "'");
        nb::gil_scoped_release released;
        otswap::io::write(path, header, length(given[0]),
                          [&given] (const std::size_t row, std::vector<double>& values,
                                    std::vector<std::int64_t>& integers) {
                            for (std::size_t c = 0; c < given.size(); ++c) {
                              if (given[c]->type == 'K') integers[c] = given[c]->integers[row];
                              else                       values[c] = given[c]->data[row];
                            }
                          }, options);
      },
      "file"_a.none(), "columns"_a.none(), nb::kw_only(), "precision"_a.none() = 9,
      "keywords"_a = nb::none(),
      "Write a table, ASCII or FITS by the extension of file. An existing file is overwritten.");

  io.def("write_displacements",
      [] (nb::handle file, nb::handle result, nb::handle groups) {
        const std::string path = to_path(file, "file");
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        static const std::vector<std::pair<std::string, otswap::io::DisplacementGroup>> names {
          {"index", otswap::io::DisplacementGroup::Index},
          {"tracer_sky", otswap::io::DisplacementGroup::TracerSky},
          {"lagrangian_sky", otswap::io::DisplacementGroup::LagrangianSky},
          {"tracer", otswap::io::DisplacementGroup::Tracer},
          {"lagrangian", otswap::io::DisplacementGroup::Lagrangian},
          {"displacement", otswap::io::DisplacementGroup::Displacement},
          {"valid_realizations", otswap::io::DisplacementGroup::ValidRealizations},
          {"selection", otswap::io::DisplacementGroup::Selection}};
        const std::vector<otswap::io::DisplacementGroup> chosen = group_list(groups, names);
        nb::gil_scoped_release released;
        otswap::io::writeDisplacements(path, r, chosen);
      },
      "file"_a.none(), "result"_a.none(), nb::kw_only(), "groups"_a = nb::none(),
      "Write the mean displacement field of a result, one row per tracer.");

  io.def("write_displacement_field",
      [] (nb::handle file, nb::handle result) {
        const std::string path = to_path(file, "file");
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        nb::gil_scoped_release released;
        otswap::io::writeDisplacementField(path, r);
      },
      "file"_a.none(), "result"_a.none(),
      "Write every realization of a result, losslessly, one row per realization and tracer.");

  io.def("write_real_space_catalog",
      [] (nb::handle file, nb::handle catalog, nb::handle groups) {
        const std::string path = to_path(file, "file");
        const PyCatalog& c = instance<PyCatalog>(catalog, "catalog", "RealSpaceCatalog");
        static const std::vector<std::pair<std::string, otswap::io::CatalogGroup>> names {
          {"index", otswap::io::CatalogGroup::Index},
          {"sky", otswap::io::CatalogGroup::Sky},
          {"cartesian", otswap::io::CatalogGroup::Cartesian},
          {"valid_realizations", otswap::io::CatalogGroup::ValidRealizations},
          {"neighbours", otswap::io::CatalogGroup::Neighbours},
          {"shift", otswap::io::CatalogGroup::Shift},
          {"status", otswap::io::CatalogGroup::Status}};
        const std::vector<otswap::io::CatalogGroup> chosen = group_list(groups, names);
        nb::gil_scoped_release released;
        otswap::io::writeRealSpaceCatalog(path, c, chosen);
      },
      "file"_a.none(), "catalog"_a.none(), nb::kw_only(), "groups"_a = nb::none(),
      "Write a catalogue moved to real space, one row per tracer.");

  io.def("write_mps_profile",
      [] (nb::handle file, nb::handle profile) {
        const std::string path = to_path(file, "file");
        const otswap::MpsProfile& p = instance<otswap::MpsProfile>(profile, "profile", "MpsProfile");
        nb::gil_scoped_release released;
        otswap::io::writeMpsProfile(path, p);
      },
      "file"_a.none(), "profile"_a.none(),
      "Write a mean particle separation profile, one row per bin.");

  io.def("read_bias_table",
      [] (nb::handle file, nb::handle delimiter, nb::handle comment) {
        const std::string path = to_path(file, "file");
        const char d = to_char(delimiter, "delimiter"), k = to_char(comment, "comment");
        return otswap::io::readBiasTable(path, d, k);
      },
      "file"_a.none(), nb::kw_only(), "delimiter"_a.none() = " ", "comment"_a.none() = "#",
      "Read a b(z) table: FITS columns REDSHIFT and BIAS, or ASCII columns 0 and 1.");
}
