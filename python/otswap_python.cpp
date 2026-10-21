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
 *  library then works on its own copy. Config, CorrectionConfig,
 *  RedshiftCut and BoxLimits are the C++ structs themselves, so that
 *  another nanobind module can take and return them; their setters convert
 *  each value to its C++ type, and the library checks it at the call.
 *  Right ascension and declination are
 *  converted from degrees here, so the library only ever sees radians; a
 *  Result records the angle unit of its call, and returns its sky arrays in
 *  it, converted once and kept. The other Result arrays are read-only views
 *  on the C++ result, kept alive by the Python object that holds it, but
 *  for ids, a new copy on each access, since assigning ids replaces them. The
 *  report of a call is written by the library into a string, and printed
 *  from Python, to sys.stdout, which a notebook shows, rather than by the
 *  library to the process's stderr; the extrapolation of b(z) is an
 *  ExtrapolationWarning instead of a line.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>

#include "otswap/OT.h"

#include "internal.h"

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

  nb::str verbosity_name (const otswap::Verbosity v)
  {
    if (v == otswap::Verbosity::Silent) return nb::str("silent");
    if (v == otswap::Verbosity::Detailed) return nb::str("detailed");
    return nb::str("normal");
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

  /// The redshift cut from None, an otswap.RedshiftCut or a (min, max)
  /// pair.
  otswap::RedshiftCut to_cut (nb::handle value)
  {
    otswap::RedshiftCut cut;
    if (value.is_none()) return cut;
    if (nb::isinstance<otswap::RedshiftCut>(value)) return nb::cast<otswap::RedshiftCut>(value);
    const std::vector<double> c = column(value, "redshift_cut");
    if (c.size() != 2)
      throw otswap::Error("redshift_cut must be an otswap.RedshiftCut or a (min, max) pair; it holds " +
                          std::to_string(c.size()) + " values");
    cut.min = c[0];
    cut.max = c[1];
    return cut;
  }

  /// A copy of the otswap.Config argument of a call, its report written to
  /// report; the argument itself is left alone.
  otswap::Config config_of (nb::handle value, std::ostream& report)
  {
    otswap::Config config = instance<otswap::Config>(value, "config", "Config");
    config.log = &report;
    return config;
  }

  /// A copy of the otswap.CorrectionConfig argument of a call, as
  /// config_of.
  otswap::CorrectionConfig correction_of (nb::handle value, std::ostream& report)
  {
    otswap::CorrectionConfig config = instance<otswap::CorrectionConfig>(value, "correction",
                                                                         "CorrectionConfig");
    config.log = &report;
    return config;
  }

  /// +inf, the default of RedshiftCut's max, created with the module: the
  /// constructor tells it from an inf given by the caller by identity.
  PyObject* cut_max_default = nullptr;

  /// A NotImplemented for the comparison operators.
  nb::object not_implemented ()
  {
    return nb::borrow(nb::handle(Py_NotImplemented));
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

  /// The IDs to store for n objects, checked as the library checks them:
  /// none for None, an empty array or exactly 0..n-1, else the values.
  std::vector<std::int64_t> ids_value (nb::handle value, const std::size_t n)
  {
    if (value.is_none()) return {};
    std::vector<std::int64_t> ids = int64_column(value, "ids");
    otswap::internal::check_ids(ids, n);
    for (std::size_t i = 0; i < ids.size(); ++i)
      if (ids[i] != (std::int64_t)i) return ids;
    return {};
  }

  /// A new read-only array of the IDs of n objects: a copy of ids, or 0..n-1
  /// when it is empty.
  nb::ndarray<nb::numpy, const std::int64_t> ids_array (const std::vector<std::int64_t>& ids,
                                                         const std::size_t n)
  {
    auto* held = new std::vector<std::int64_t>(ids);
    if (held->empty()) {
      held->resize(n);
      std::iota(held->begin(), held->end(), std::int64_t(0));
    }
    nb::capsule owner(held, [] (void* p) noexcept { delete static_cast<std::vector<std::int64_t>*>(p); });
    return nb::ndarray<nb::numpy, const std::int64_t>(held->data(), {n}, owner);
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
  /// first access and kept until a filter changes them. Its config.log is
  /// &std::clog, never the stream a call wrote its report into.
  struct PyResult : otswap::Result {
    bool degrees = false;
    nb::object tracersSkyInUnit, lagrangianSkyInUnit;

    PyResult () = default;
    PyResult (otswap::Result&& r, const bool inDegrees) : otswap::Result(std::move(r)), degrees(inDegrees)
    {
      config.log = &std::clog;
    }
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

  /// A corner of a box: three numbers, x, y, z.
  std::array<double, 3> corner (nb::handle value, const std::string& name)
  {
    const std::vector<double> c = column(value, name);
    if (c.size() != 3)
      throw otswap::Error(name + " must hold three numbers, x, y, z; it holds " +
                          std::to_string(c.size()));
    return {c[0], c[1], c[2]};
  }

  /// The box limits from None, an otswap.BoxLimits, or a pair of corners
  /// ((x0, y0, z0), (x1, y1, z1)). The values are checked by the library.
  std::optional<otswap::BoxLimits> to_box (nb::handle value)
  {
    if (value.is_none()) return std::nullopt;
    if (nb::isinstance<otswap::BoxLimits>(value)) return nb::cast<otswap::BoxLimits>(value);
    const CArray a = as_float64(value, "box");
    if (a.ndim() != 2 || a.shape(0) != 2 || a.shape(1) != 3)
      throw otswap::Error("box must be an otswap.BoxLimits or a pair of corners ((x0, y0, z0), "
                          "(x1, y1, z1)); it has shape " + shape_of(a));
    otswap::BoxLimits box;
    for (std::size_t k = 0; k < 3; ++k) {
      box.lower[k] = a.data()[k];
      box.upper[k] = a.data()[3 + k];
    }
    return box;
  }

  nb::object corner_tuple (const std::array<double, 3>& c)
  {
    return nb::make_tuple(c[0], c[1], c[2]);
  }

  nb::object box_or_none (const std::optional<otswap::BoxLimits>& box)
  {
    if (!box) return nb::none();
    return nb::cast(*box);
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

  // ------------------------------------------------------------------ Config

  const otswap::Config configDefaults;

  nb::class_<otswap::Config>(m, "Config",
      "Free parameters of the reconstruction, as otswap::Config. Its fields are checked by the "
      "call that reads them, not here.")
    .def("__init__",
        [] (otswap::Config* self, nb::handle nRealizations, nb::handle convergence, nb::handle seed,
            nb::handle cellSize, nb::handle rejectCrossings, nb::handle maxUnobservedPixelsCrossed,
            nb::handle verbosity) {
          otswap::Config config;
          config.nRealizations = to_unsigned(nRealizations, "n_realizations");
          config.convergence = to_double(convergence, "convergence");
          config.seed = to_unsigned(seed, "seed");
          config.cellSize = to_double(cellSize, "cell_size");
          config.rejectCrossings = to_bool(rejectCrossings, "reject_crossings");
          config.maxUnobservedPixelsCrossed =
            to_unsigned(maxUnobservedPixelsCrossed, "max_unobserved_pixels_crossed");
          config.verbosity = to_verbosity(verbosity);
          new (self) otswap::Config(config);
        },
        nb::kw_only(), "n_realizations"_a.none() = configDefaults.nRealizations,
        "convergence"_a.none() = configDefaults.convergence, "seed"_a.none() = configDefaults.seed,
        "cell_size"_a.none() = configDefaults.cellSize,
        "reject_crossings"_a.none() = configDefaults.rejectCrossings,
        "max_unobserved_pixels_crossed"_a.none() = configDefaults.maxUnobservedPixelsCrossed,
        "verbosity"_a.none() = verbosity_name(configDefaults.verbosity),
        "Every field a keyword, with the defaults of otswap::Config.")
    .def_prop_rw("n_realizations",
        [] (const otswap::Config& c) { return c.nRealizations; },
        [] (otswap::Config& c, nb::handle v) { c.nRealizations = to_unsigned(v, "n_realizations"); },
        nb::for_setter(nb::arg("value").none()), "Independent reconstructions to run.")
    .def_prop_rw("convergence",
        [] (const otswap::Config& c) { return c.convergence; },
        [] (otswap::Config& c, nb::handle v) { c.convergence = to_double(v, "convergence"); },
        nb::for_setter(nb::arg("value").none()),
        "The sweeps stop once the fraction of successful swaps in a sweep no longer exceeds it.")
    .def_prop_rw("seed",
        [] (const otswap::Config& c) { return c.seed; },
        [] (otswap::Config& c, nb::handle v) { c.seed = to_unsigned(v, "seed"); },
        nb::for_setter(nb::arg("value").none()), "Seed of the generator; 0 draws one.")
    .def_prop_rw("cell_size",
        [] (const otswap::Config& c) { return c.cellSize; },
        [] (otswap::Config& c, nb::handle v) { c.cellSize = to_double(v, "cell_size"); },
        nb::for_setter(nb::arg("value").none()),
        "Grid cell size, in units of the representative mean particle separation. Affects speed "
        "only.")
    .def_prop_rw("reject_crossings",
        [] (const otswap::Config& c) { return c.rejectCrossings; },
        [] (otswap::Config& c, nb::handle v) { c.rejectCrossings = to_bool(v, "reject_crossings"); },
        nb::for_setter(nb::arg("value").none()),
        "With a mask, filter the result with reject_mask_crossings before returning it.")
    .def_prop_rw("max_unobserved_pixels_crossed",
        [] (const otswap::Config& c) { return c.maxUnobservedPixelsCrossed; },
        [] (otswap::Config& c, nb::handle v) {
          c.maxUnobservedPixelsCrossed = to_unsigned(v, "max_unobserved_pixels_crossed");
        },
        nb::for_setter(nb::arg("value").none()), "The threshold of that filter.")
    .def_prop_rw("verbosity",
        [] (const otswap::Config& c) { return verbosity_name(c.verbosity); },
        [] (otswap::Config& c, nb::handle v) { c.verbosity = to_verbosity(v); },
        nb::for_setter(nb::arg("value").none()), "What the call reports.")
    .def("__repr__", [] (const otswap::Config& c) {
        const std::string text =
          "Config(n_realizations=" + std::to_string(c.nRealizations) +
          ", convergence=" + repr_of(nb::float_(c.convergence)) + ", seed=" + std::to_string(c.seed) +
          ", cell_size=" + repr_of(nb::float_(c.cellSize)) +
          ", reject_crossings=" + (c.rejectCrossings ? "True" : "False") +
          ", max_unobserved_pixels_crossed=" + std::to_string(c.maxUnobservedPixelsCrossed) +
          ", verbosity=" + repr_of(verbosity_name(c.verbosity)) + ")";
        return nb::str(text.c_str());
      })
    .def("__eq__", [] (const otswap::Config& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<otswap::Config>(other)) return not_implemented();
        const otswap::Config& b = nb::cast<const otswap::Config&>(other);
        return nb::bool_(a.nRealizations == b.nRealizations && a.convergence == b.convergence &&
                         a.seed == b.seed && a.cellSize == b.cellSize &&
                         a.rejectCrossings == b.rejectCrossings &&
                         a.maxUnobservedPixelsCrossed == b.maxUnobservedPixelsCrossed &&
                         a.verbosity == b.verbosity);
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  const otswap::CorrectionConfig correctionDefaults;

  nb::class_<otswap::CorrectionConfig>(m, "CorrectionConfig",
      "Options of real_space_lightcone and real_space_box, as otswap::CorrectionConfig.")
    .def("__init__",
        [] (otswap::CorrectionConfig* self, nb::handle weightByRealizations, nb::handle verbosity) {
          otswap::CorrectionConfig config;
          config.weightByRealizations = to_bool(weightByRealizations, "weight_by_realizations");
          config.verbosity = to_verbosity(verbosity);
          new (self) otswap::CorrectionConfig(config);
        },
        nb::kw_only(), "weight_by_realizations"_a.none() = correctionDefaults.weightByRealizations,
        "verbosity"_a.none() = verbosity_name(correctionDefaults.verbosity),
        "Every field a keyword, with the defaults of otswap::CorrectionConfig.")
    .def_prop_rw("weight_by_realizations",
        [] (const otswap::CorrectionConfig& c) { return c.weightByRealizations; },
        [] (otswap::CorrectionConfig& c, nb::handle v) {
          c.weightByRealizations = to_bool(v, "weight_by_realizations");
        },
        nb::for_setter(nb::arg("value").none()),
        "Weight each neighbour by its count of valid realizations as well.")
    .def_prop_rw("verbosity",
        [] (const otswap::CorrectionConfig& c) { return verbosity_name(c.verbosity); },
        [] (otswap::CorrectionConfig& c, nb::handle v) { c.verbosity = to_verbosity(v); },
        nb::for_setter(nb::arg("value").none()), "What the call reports.")
    .def("__repr__", [] (const otswap::CorrectionConfig& c) {
        const std::string text =
          std::string("CorrectionConfig(weight_by_realizations=") +
          (c.weightByRealizations ? "True" : "False") +
          ", verbosity=" + repr_of(verbosity_name(c.verbosity)) + ")";
        return nb::str(text.c_str());
      })
    .def("__eq__", [] (const otswap::CorrectionConfig& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<otswap::CorrectionConfig>(other)) return not_implemented();
        const otswap::CorrectionConfig& b = nb::cast<const otswap::CorrectionConfig&>(other);
        return nb::bool_(a.weightByRealizations == b.weightByRealizations && a.verbosity == b.verbosity);
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  // ------------------------------------------------------------- RedshiftCut

  cut_max_default = PyFloat_FromDouble(otswap::RedshiftCut().max);
  if (cut_max_default == nullptr) throw nb::python_error();

  nb::class_<otswap::RedshiftCut>(m, "RedshiftCut",
      "A closed redshift range for reconstruct_lightcone and measure_mps_profile, as "
      "otswap::RedshiftCut. The default, (-inf, +inf), cuts nothing. Wherever a RedshiftCut is "
      "taken, a (min, max) pair is accepted too.")
    .def("__init__", [] (otswap::RedshiftCut* self, nb::handle min, nb::handle max) {
          otswap::RedshiftCut cut;
          // A pair alone: the implicit conversion, or RedshiftCut((min, max)).
          std::pair<double, double> bounds;
          if (max.ptr() == cut_max_default && nb::try_cast(min, bounds)) {
            cut.min = bounds.first;
            cut.max = bounds.second;
          }
          else {
            cut.min = to_double(min, "min");
            cut.max = to_double(max, "max");
          }
          new (self) otswap::RedshiftCut(cut);
        }, "min"_a.none() = otswap::RedshiftCut().min, "max"_a.none() = nb::handle(cut_max_default),
        "The range [min, max]; its bounds are checked by the call that reads it.")
    .def_prop_rw("min",
        [] (const otswap::RedshiftCut& c) { return c.min; },
        [] (otswap::RedshiftCut& c, nb::handle v) { c.min = to_double(v, "min"); },
        nb::for_setter(nb::arg("value").none()), "The lowest redshift kept.")
    .def_prop_rw("max",
        [] (const otswap::RedshiftCut& c) { return c.max; },
        [] (otswap::RedshiftCut& c, nb::handle v) { c.max = to_double(v, "max"); },
        nb::for_setter(nb::arg("value").none()), "The highest redshift kept.")
    .def("__repr__", [] (const otswap::RedshiftCut& c) {
        const std::string text = "RedshiftCut(min=" + repr_of(nb::float_(c.min)) + ", max=" +
                                 repr_of(nb::float_(c.max)) + ")";
        return nb::str(text.c_str());
      })
    .def("__eq__", [] (const otswap::RedshiftCut& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<otswap::RedshiftCut>(other)) return not_implemented();
        const otswap::RedshiftCut& b = nb::cast<const otswap::RedshiftCut&>(other);
        return nb::bool_(a.min == b.min && a.max == b.max);
      }, "other"_a.none())
    .attr("__hash__") = nb::none();
  nb::implicitly_convertible<std::pair<double, double>, otswap::RedshiftCut>();

  // --------------------------------------------------------------- BoxLimits

  nb::class_<otswap::BoxLimits>(m, "BoxLimits",
      "The closed limits of a box for reconstruct_box, in Mpc/h. Read-only. Wherever a BoxLimits "
      "is taken, a pair of corners ((x0, y0, z0), (x1, y1, z1)) is accepted too.")
    .def("__init__", [] (otswap::BoxLimits* self, nb::handle lower, nb::handle upper) {
          otswap::BoxLimits box;
          box.lower = corner(lower, "lower");
          box.upper = corner(upper, "upper");
          new (self) otswap::BoxLimits(box);
        }, "lower"_a.none(), "upper"_a.none(),
        "The box from its lower and upper corners, (x, y, z) each.")
    .def("__init__",
        [] (otswap::BoxLimits* self, const std::pair<std::array<double, 3>, std::array<double, 3>>& corners) {
          otswap::BoxLimits box;
          box.lower = corners.first;
          box.upper = corners.second;
          new (self) otswap::BoxLimits(box);
        }, "corners"_a, "The box from a pair of corners, the implicit conversion.")
    .def_prop_ro("lower", [] (const otswap::BoxLimits& b) { return corner_tuple(b.lower); },
                 "x, y, z of the lower corner.")
    .def_prop_ro("upper", [] (const otswap::BoxLimits& b) { return corner_tuple(b.upper); },
                 "x, y, z of the upper corner.")
    .def("__repr__", [] (const otswap::BoxLimits& b) {
        const std::string text = "BoxLimits(lower=" + repr_of(corner_tuple(b.lower)) + ", upper=" +
                                 repr_of(corner_tuple(b.upper)) + ")";
        return nb::str(text.c_str());
      })
    .def("__eq__", [] (const otswap::BoxLimits& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<otswap::BoxLimits>(other)) return not_implemented();
        const otswap::BoxLimits& b = nb::cast<const otswap::BoxLimits&>(other);
        return nb::bool_(a.lower == b.lower && a.upper == b.upper);
      }, "other"_a.none())
    .attr("__hash__") = nb::none();
  nb::implicitly_convertible<std::pair<std::array<double, 3>, std::array<double, 3>>, otswap::BoxLimits>();

  // --------------------------------------------------------- SelectionCounts

  nb::class_<otswap::SelectionCounts>(m, "SelectionCounts",
      "What a reconstruction left out, and why. Read-only.")
    .def_prop_ro("redshift_cut", [] (const otswap::SelectionCounts& c) -> nb::object {
        if (!c.redshiftCut) return nb::none();
        return nb::cast(*c.redshiftCut);
      }, "The redshift cut applied, a new RedshiftCut on each access, or None when no cut was "
         "applied.")
    .def_ro("mask_applied", &otswap::SelectionCounts::maskApplied)
    .def_prop_ro("box_limits", [] (const otswap::SelectionCounts& c) { return box_or_none(c.boxLimits); },
                 "The box limits given to reconstruct_box, or None when none were.")
    .def_ro("tracers", &otswap::SelectionCounts::tracers)
    .def_ro("tracers_outside_redshift_cut", &otswap::SelectionCounts::tracersOutsideRedshiftCut)
    .def_ro("tracers_outside_mask", &otswap::SelectionCounts::tracersOutsideMask)
    .def_ro("tracers_outside_both", &otswap::SelectionCounts::tracersOutsideBoth)
    .def_ro("tracers_outside_box", &otswap::SelectionCounts::tracersOutsideBox)
    .def_ro("randoms", &otswap::SelectionCounts::randoms)
    .def_ro("randoms_outside_redshift_cut", &otswap::SelectionCounts::randomsOutsideRedshiftCut)
    .def_ro("randoms_outside_mask", &otswap::SelectionCounts::randomsOutsideMask)
    .def_ro("randoms_outside_both", &otswap::SelectionCounts::randomsOutsideBoth)
    .def_ro("randoms_outside_box", &otswap::SelectionCounts::randomsOutsideBox)
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
            nb::handle outsideRedshiftCut, nb::handle outsideMask, nb::handle outsideBox,
            nb::handle ids, nb::handle config) {
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
          if (!outsideBox.is_none() && !tracersSky.is_none())
            throw otswap::Error("outside_box describes a box result, and tracers_sky a lightcone "
                                "one; give outside_box without tracers_sky");
          r.outsideBox = outsideBox.is_none()
            ? std::vector<std::uint8_t>(n, 0) : flags(outsideBox, "outside_box", {n});
          r.ids = ids_value(ids, n);
          if (config.is_none()) r.config.nRealizations = (unsigned)nRec;
          else {
            r.config = instance<otswap::Config>(config, "config", "Config");
            if (r.config.nRealizations != nRec)
              throw otswap::Error("config.n_realizations is " + std::to_string(r.config.nRealizations) +
                                  " and the arrays hold " + std::to_string(nRec) +
                                  " realizations; they must be equal");
          }
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
        "outside_mask"_a = nb::none(), "outside_box"_a = nb::none(), "ids"_a = nb::none(),
        "config"_a = nb::none(),
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
    .def_prop_ro("config", [] (const PyResult& r) {
        otswap::Config config = r.config;
        config.log = &std::clog;
        return config;
      }, "The parameters of the run: a new copy on each access, so changing it changes nothing in "
         "the result. seed is the one used: the drawn one when 0 was given.")
    .def_prop_ro("distances", [] (const PyResult& r) -> nb::object {
        if (!r.distances) return nb::none();
        return nb::cast(otswap::DistanceTable(*r.distances));
      }, "The distance table of a lightcone result, a new copy on each access; None for a box.")
    .def_prop_rw("ids",
        [] (const PyResult& r) { return ids_array(r.ids, r.nObjects); },
        [] (PyResult& r, nb::handle value) { r.ids = ids_value(value, r.nObjects); },
        nb::for_getter(nb::rv_policy::reference), nb::for_setter(nb::arg("value").none()),
        "The ID of each object, shape (n_objects,), a new read-only copy on each access; "
        "0..n_objects-1 when none were set. Settable: one 64-bit integer per object, no value twice, "
        "or None.")
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
    .def_prop_ro("outside_box", [] (ResultHandle self) {
        const PyResult& r = result_of(self);
        return view<bool>(self, reinterpret_cast<const bool*>(r.outsideBox.data()), {r.nObjects});
      }, "Whether each tracer lay outside the box limits, shape (n_objects,).")
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
    .def_prop_ro("box", [] (ResultHandle self) { return box_or_none(result_of(self).box); },
                 "Box: the box of the reconstruction, a copy; None for a lightcone.")
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
    .def_prop_ro("redshift", [] (nb::handle_t<otswap::DistanceTable> self) {
        const auto& t = nb::cast<const otswap::DistanceTable&>(self);
        return view<double>(self, t.redshift().data(), {t.redshift().size()});
      }, "The sampled redshifts, shape (n_samples,); a read-only view.")
    .def_prop_ro("distance", [] (nb::handle_t<otswap::DistanceTable> self) {
        const auto& t = nb::cast<const otswap::DistanceTable&>(self);
        return view<double>(self, t.distance().data(), {t.distance().size()});
      }, "The comoving distance at each sampled redshift, Mpc/h; a read-only view.")
    .def_prop_ro("growth_rate", [] (nb::handle_t<otswap::DistanceTable> self) -> nb::object {
        const auto& t = nb::cast<const otswap::DistanceTable&>(self);
        if (!t.hasGrowthRate()) return nb::none();
        return nb::cast(view<double>(self, t.growthRate().data(), {t.growthRate().size()}));
      }, "The growth rate at each sampled redshift; None when the table carries none.")
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
      [] (nb::handle tracers, nb::handle randoms, nb::handle mps, nb::handle box, nb::handle configArg) {
        const std::vector<double> t = rows3(tracers, "tracers");
        std::vector<double> r;
        if (!randoms.is_none()) {
          r = rows3(randoms, "randoms");
          if (r.empty())
            throw otswap::Error("randoms is empty; omit it to have the randoms drawn in the box "
                                "limits, or in the bounding box of the tracers");
        }
        const bool computed = mps.is_none();
        const double separation = computed ? 0. : to_double(mps, "mps");
        const std::optional<otswap::BoxLimits> limits = to_box(box);
        std::ostringstream report;
        const otswap::Config config = config_of(configArg, report);
        otswap::Result result;
        {
          nb::gil_scoped_release released;
          result = computed ? otswap::reconstructBox(t, r, config, limits)
                            : otswap::reconstructBox(t, r, separation, config, limits);
        }
        print_report(report.str());
        return PyResult(std::move(result), false);
      },
      "tracers"_a.none(), "randoms"_a = nb::none(), nb::kw_only(), "mps"_a = nb::none(),
      "box"_a = nb::none(), "config"_a.none() = otswap::Config(),
      "Reconstruct in box geometry, with a constant mean particle separation, given or that of "
      "the tracers' bounding box.");

  m.def("reconstruct_lightcone",
      [] (nb::handle tracersSky, nb::handle randomsSky, nb::handle skyAreaDeg2, nb::handle mask,
          nb::handle nBins, nb::handle distances, nb::handle angleUnit, nb::handle tracers,
          nb::handle randoms, nb::handle redshiftCut, nb::handle configArg) {
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
        std::ostringstream report;
        const otswap::Config config = config_of(configArg, report);
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
      "tracers"_a = nb::none(), "randoms"_a = nb::none(), "redshift_cut"_a = nb::none(),
      "config"_a.none() = otswap::Config(),
      "Reconstruct in lightcone geometry.");

  m.def("measure_mps_profile",
      [] (nb::handle tracersSky, nb::handle skyAreaDeg2, nb::handle mask, nb::handle nBins,
          nb::handle distances, nb::handle redshiftCut, nb::handle angleUnit, nb::handle verbosity) {
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const bool degrees = in_degrees(angleUnit);
        std::vector<double> ts = rows3(tracersSky, "tracers_sky");
        if (degrees) sky_to_radians(ts, "tracers_sky");
        if (!skyAreaDeg2.is_none() && !mask.is_none())
          throw otswap::Error("give sky_area_deg2 or mask, not both; with a mask the sky area is "
                              "the mask's");
        if (skyAreaDeg2.is_none() && mask.is_none())
          throw otswap::Error("give sky_area_deg2 or mask; one is required");
        const otswap::Mask* footprint = mask.is_none() ? nullptr
                                       : &instance<otswap::Mask>(mask, "mask", "Mask");
        const double area = footprint ? 0. : to_double(skyAreaDeg2, "sky_area_deg2");
        const unsigned bins = to_unsigned(nBins, "n_bins");
        const otswap::Verbosity level = to_verbosity(verbosity);
        std::ostringstream report;
        const otswap::RedshiftCut cut = to_cut(redshiftCut);

        otswap::MpsProfile profile;
        {
          nb::gil_scoped_release released;
          profile = footprint ? otswap::measureMpsProfile(ts, *footprint, bins, table, cut, level, &report)
                              : otswap::measureMpsProfile(ts, area, bins, table, cut, level, &report);
        }
        print_report(report.str());
        return profile;
      },
      "tracers_sky"_a.none(), nb::kw_only(), "sky_area_deg2"_a = nb::none(), "mask"_a = nb::none(),
      "n_bins"_a.none(), "distances"_a.none(), "redshift_cut"_a = nb::none(), "angle_unit"_a.none(),
      "verbosity"_a.none() = "normal",
      "Measure the mps(z) profile of a lightcone without reconstructing.");

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
    .def_prop_ro("ids", [] (const PyCatalog& c) { return ids_array(c.ids, c.nObjects); },
                 nb::rv_policy::reference,
                 "The IDs of the objects, copied from the Result corrected, shape (n_objects,), a new "
                 "read-only copy on each access; 0..n_objects-1 when it had none.")
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
          nb::handle correction) {
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const otswap::BiasTable& b = instance<otswap::BiasTable>(bias, "bias", "BiasTable");
        const double s = to_double(sigma, "sigma");
        std::ostringstream report;
        const otswap::CorrectionConfig config = correction_of(correction, report);
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
      "correction"_a.none() = otswap::CorrectionConfig(),
      "Move a lightcone catalogue from redshift space to real space.");

  m.def("real_space_box",
      [] (nb::handle result, nb::handle axis, nb::handle redshift, nb::handle distances,
          nb::handle bias, nb::handle sigma, nb::handle correction) {
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        const otswap::DistanceTable& table = instance<otswap::DistanceTable>(distances, "distances", "DistanceTable");
        const unsigned a = to_unsigned(axis, "axis");
        const double z = to_double(redshift, "redshift");
        const double b = to_double(bias, "bias");
        const double s = to_double(sigma, "sigma");
        std::ostringstream report;
        const otswap::CorrectionConfig config = correction_of(correction, report);
        otswap::RealSpaceCatalog catalog;
        {
          nb::gil_scoped_release released;
          catalog = otswap::realSpaceBox(r, a, z, table, b, s, config);
        }
        print_report(report.str());
        return PyCatalog(std::move(catalog), false);
      },
      "result"_a.none(), nb::kw_only(), "axis"_a.none(), "redshift"_a.none(), "distances"_a.none(),
      "bias"_a.none(), "sigma"_a.none(), "correction"_a.none() = otswap::CorrectionConfig(),
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

  io.def("write_displacement_mean",
      [] (nb::handle file, nb::handle result, nb::handle groups) {
        const std::string path = to_path(file, "file");
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        static const std::vector<std::pair<std::string, otswap::io::DisplacementGroup>> names {
          {"id", otswap::io::DisplacementGroup::Id},
          {"tracer_sky", otswap::io::DisplacementGroup::TracerSky},
          {"lagrangian_sky", otswap::io::DisplacementGroup::LagrangianSky},
          {"tracer", otswap::io::DisplacementGroup::Tracer},
          {"lagrangian", otswap::io::DisplacementGroup::Lagrangian},
          {"displacement", otswap::io::DisplacementGroup::Displacement},
          {"valid_realizations", otswap::io::DisplacementGroup::ValidRealizations},
          {"selection", otswap::io::DisplacementGroup::Selection},
          {"outside_box", otswap::io::DisplacementGroup::OutsideBox}};
        const std::vector<otswap::io::DisplacementGroup> chosen = group_list(groups, names);
        nb::gil_scoped_release released;
        otswap::io::writeDisplacementMean(path, r, chosen);
      },
      "file"_a.none(), "result"_a.none(), nb::kw_only(), "groups"_a = nb::none(),
      "Write the mean displacement field of a result, one row per tracer.");

  io.def("write_displacement_complete",
      [] (nb::handle file, nb::handle result, nb::handle lightconeTables) {
        const std::string path = to_path(file, "file");
        const PyResult& r = instance<PyResult>(result, "result", "Result");
        const bool tables = to_bool(lightconeTables, "lightcone_tables");
        nb::gil_scoped_release released;
        otswap::io::writeDisplacementComplete(path, r, tables);
      },
      "file"_a.none(), "result"_a.none(), nb::kw_only(), "lightcone_tables"_a.none() = true,
      "Write every realization of a result, losslessly, one row per realization and tracer, with "
      "the mps(z) profile and the distance table of a lightcone.");

  io.def("read_displacement_complete",
      [] (nb::handle file, nb::handle angleUnit) {
        const std::string path = to_path(file, "file");
        const bool degrees = in_degrees(angleUnit);
        otswap::Result result;
        {
          nb::gil_scoped_release released;
          result = otswap::io::readDisplacementComplete(path);
        }
        return PyResult(std::move(result), degrees);
      },
      "file"_a.none(), nb::kw_only(), "angle_unit"_a.none(),
      "Read a file written by write_displacement_complete back into a Result.");

  io.def("write_real_space_catalog",
      [] (nb::handle file, nb::handle catalog, nb::handle groups) {
        const std::string path = to_path(file, "file");
        const PyCatalog& c = instance<PyCatalog>(catalog, "catalog", "RealSpaceCatalog");
        static const std::vector<std::pair<std::string, otswap::io::CatalogGroup>> names {
          {"id", otswap::io::CatalogGroup::Id},
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
