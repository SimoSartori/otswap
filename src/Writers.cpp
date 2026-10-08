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
 *  @file src/Writers.cpp
 *
 *  @brief The writers of a Result, a RealSpaceCatalog and an MpsProfile.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

#include "internal.h"

namespace {

  using otswap::io::Column;
  using otswap::io::Keyword;

  /// Significant digits of the writers that are not lossless, and of the
  /// lossless one.
  constexpr int kPrecision = 9;
  constexpr int kLossless = 17;

  /// A set of columns and the function that fills their entries of a row,
  /// from the given offset of the two buffers.
  struct Group {
    std::vector<Column> columns;
    std::function<void(std::size_t row, double* values, std::int64_t* integers)> fill;
  };

  Column column (const std::string& name, const char type, const std::string& description,
                 const std::string& unit = "")
  {
    return Column{name, type, description, unit, {}, {}};
  }

  /// The shortest decimal of up to 17 significant digits that reads back as
  /// value.
  std::string number (const double value)
  {
    for (int digits = 15; ; ++digits) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      out << std::setprecision(digits) << value;
      if (digits == kLossless || std::strtod(out.str().c_str(), nullptr) == value) return out.str();
    }
  }

  template <typename T>
  std::string integer (const T value)
  {
    return std::to_string(value);
  }

  /// The size a field must have, checked before any row is written.
  void check_size (const std::size_t size, const std::size_t expected, const std::string& field)
  {
    if (size != expected)
      throw otswap::Error("cannot write " + field + ": it holds " + std::to_string(size) +
                          " entries where " + std::to_string(expected) + " are expected");
  }

  const char* geometry_name (const otswap::Geometry g)
  {
    return g == otswap::Geometry::Box ? "box" : "lightcone";
  }

  /// The keywords every product starts with.
  std::vector<Keyword> common (const std::string& product, const otswap::Geometry geometry,
                               const std::size_t nObjects)
  {
    return {{"PRODUCT", product, "what the file holds"},
            {"OTSWAPV", OTSWAP_VERSION, "otswap version"},
            {"GEOMETRY", geometry_name(geometry), "geometry of the reconstruction"},
            {"NOBJECTS", integer(nObjects), "rows: one per input object"}};
  }

  /// The keywords of a reconstruction.
  std::vector<Keyword> result_keywords (const std::string& product, const otswap::Result& r)
  {
    std::vector<Keyword> k = common(product, r.geometry, r.nObjects);
    k.push_back({"NREC", integer(r.nRealizations), "realizations"});
    k.push_back({"SEED", integer(r.config.seed), "seed of the random streams"});
    k.push_back({"CONVERG", number(r.config.convergence), "convergence threshold of the swap loop"});
    if (r.geometry == otswap::Geometry::Box && !std::isnan(r.mps))
      k.push_back({"MPS", number(r.mps), "mean particle separation, Mpc/h"});
    const otswap::SelectionCounts& s = r.selection;
    if (s.redshiftCut && std::isfinite(s.redshiftCut->min))
      k.push_back({"ZCUTMIN", number(s.redshiftCut->min), "lower bound of the redshift cut"});
    if (s.redshiftCut && std::isfinite(s.redshiftCut->max))
      k.push_back({"ZCUTMAX", number(s.redshiftCut->max), "upper bound of the redshift cut"});
    if (s.maskApplied) k.push_back({"MASK", "1", "objects selected by a mask"});
    if (r.filteredNside != 0)
      k.push_back({"FILTNSID", integer(r.filteredNside), "NSIDE of the mask filter"});
    if (s.maxUnobservedPixelsCrossed)
      k.push_back({"MAXPIXCR", integer(*s.maxUnobservedPixelsCrossed),
                   "unobserved pixels a displacement may cross"});
    return k;
  }

  /// A 3-column group over a flat [object][3] array; angles converted to
  /// degrees when the array is a sky one.
  Group triple (const std::vector<double>& a, const std::string& prefix, const char* const suffix[3],
                const std::string& description, const char* const unit[3], const bool sky)
  {
    Group g;
    for (int c = 0; c < 3; ++c)
      g.columns.push_back(column(prefix + suffix[c], 'D', description, unit[c]));
    const std::vector<double>* source = &a;
    g.fill = [source, sky] (const std::size_t row, double* values, std::int64_t*) {
      const double* p = &(*source)[3*row];
      if (!sky) {
        values[0] = p[0]; values[1] = p[1]; values[2] = p[2];
        return;
      }
      const std::vector<double> degrees =
        otswap::skyToDegrees({otswap::internal::fold_ra(p[0]), p[1], p[2]});
      values[0] = degrees[0]; values[1] = degrees[1]; values[2] = degrees[2];
    };
    return g;
  }

  const char* const kXyz[3] = {"X", "Y", "Z"};
  const char* const kSky[3] = {"RA", "Dec", "Red"};
  const char* const kMpc[3] = {"Mpc/h", "Mpc/h", "Mpc/h"};
  const char* const kDeg[3] = {"deg", "deg", ""};

  Group index_group ()
  {
    return {{column("index", 'K', "row of the tracer in the input catalogue, from 0")},
            [] (const std::size_t row, double*, std::int64_t* integers) {
              integers[0] = (std::int64_t)row;
            }};
  }

  /// Write the groups, nRows rows.
  void write_groups (const std::string& file, const std::vector<Group>& groups,
                     const std::size_t nRows, const otswap::io::WriteOptions& options)
  {
    std::vector<Column> columns;
    std::vector<std::size_t> offset;
    for (const Group& g : groups) {
      offset.push_back(columns.size());
      columns.insert(columns.end(), g.columns.begin(), g.columns.end());
    }
    otswap::io::write(file, columns, nRows,
                      [&groups, &offset] (const std::size_t row, std::vector<double>& values,
                                          std::vector<std::int64_t>& integers) {
                        for (std::size_t k = 0; k < groups.size(); ++k)
                          groups[k].fill(row, values.data() + offset[k], integers.data() + offset[k]);
                      }, options);
  }

  template <typename G>
  void check_groups (const std::vector<G>& groups, const std::vector<G>& lightconeOnly,
                     const bool lightcone, const char* const names[])
  {
    for (std::size_t i = 0; i < groups.size(); ++i) {
      for (std::size_t j = 0; j < i; ++j)
        if (groups[j] == groups[i])
          throw otswap::Error(std::string("the group ") + names[(int)groups[i]] + " is given twice");
      if (!lightcone)
        for (const G only : lightconeOnly)
          if (groups[i] == only)
            throw otswap::Error(std::string("the group ") + names[(int)groups[i]] +
                                " describes a lightcone, and this is a box");
    }
  }

}


// ============================================================================


void otswap::io::writeDisplacements (const std::string& file, const Result& result,
                                     const std::vector<DisplacementGroup>& groups)
{
  using G = DisplacementGroup;
  static const char* const names[] = {"Index", "TracerSky", "LagrangianSky", "Tracer", "Lagrangian",
                                      "Displacement", "ValidRealizations", "Selection"};
  const bool lightcone = result.geometry == Geometry::Lightcone;
  const std::vector<G> chosen = !groups.empty() ? groups
    : lightcone ? std::vector<G>{G::TracerSky, G::LagrangianSky, G::Tracer, G::Lagrangian,
                                 G::Displacement, G::ValidRealizations, G::Selection}
                : std::vector<G>{G::Tracer, G::Lagrangian, G::Displacement, G::ValidRealizations};
  check_groups(chosen, {G::TracerSky, G::LagrangianSky, G::Selection}, lightcone, names);
  internal::check_flags(result);

  const std::size_t n = result.nObjects;
  std::vector<Group> out;
  for (const G g : chosen) {
    switch (g) {
    case G::Index:
      out.push_back(index_group());
      break;
    case G::TracerSky:
      check_size(result.tracersSky.size(), 3 * n, "tracersSky");
      out.push_back(triple(result.tracersSky, "trac", kSky,
                           "observed right ascension, declination and redshift of the tracer",
                           kDeg, true));
      break;
    case G::LagrangianSky:
      check_size(result.lagrangianSky.size(), 3 * n, "lagrangianSky");
      out.push_back(triple(result.lagrangianSky, "lagr", kSky,
                           "sky coordinates of the mean Lagrangian position; redshift NaN beyond "
                           "the distance table", kDeg, true));
      break;
    case G::Tracer:
      check_size(result.tracers.size(), 3 * n, "tracers");
      out.push_back(triple(result.tracers, "trac", kXyz, "comoving position of the tracer", kMpc,
                           false));
      break;
    case G::Lagrangian:
      check_size(result.lagrangian.size(), 3 * n, "lagrangian");
      out.push_back(triple(result.lagrangian, "lagr", kXyz,
                           "mean Lagrangian position, tracer + mean displacement", kMpc, false));
      break;
    case G::Displacement:
      check_size(result.meanDisplacement.size(), 3 * n, "meanDisplacement");
      out.push_back(triple(result.meanDisplacement, "displ", kXyz,
                           "mean displacement over the valid realizations", kMpc, false));
      break;
    case G::ValidRealizations: {
      check_size(result.validRealizations.size(), n, "validRealizations");
      const std::vector<unsigned>* v = &result.validRealizations;
      out.push_back({{column("nValidRec", 'J', "number of valid realizations of the tracer")},
                     [v] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = (double)(*v)[row];
                     }});
      break;
    }
    case G::Selection: {
      const Result* r = &result;
      out.push_back({{column("outsideRedshiftCut", 'J', "1 if left out of the reconstruction by the redshift cut"),
                      column("outsideMask", 'J', "1 if left out of the reconstruction by the mask")},
                     [r] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = r->outsideRedshiftCut.empty() ? 0. : (double)r->outsideRedshiftCut[row];
                       values[1] = r->outsideMask.empty() ? 0. : (double)r->outsideMask[row];
                     }});
      break;
    }
    }
  }

  WriteOptions options;
  options.precision = kPrecision;
  options.keywords = result_keywords("displacements", result);
  write_groups(file, out, n, options);
}


// ============================================================================


void otswap::io::writeDisplacementField (const std::string& file, const Result& result)
{
  const std::size_t n = result.nObjects;
  const std::size_t rows = (std::size_t)result.nRealizations * n;
  check_size(result.tracers.size(), 3 * n, "tracers");
  check_size(result.displacement.size(), 3 * rows, "displacement");
  check_size(result.matchedRandom.size(), 3 * rows, "matchedRandom");
  check_size(result.valid.size(), rows, "valid");
  internal::check_flags(result);

  const Result* r = &result;
  const std::vector<Column> columns {
    column("realization", 'J', "realization, from 0"),
    column("index", 'K', "row of the tracer in the input catalogue, from 0"),
    column("tracX", 'D', "comoving position of the tracer", "Mpc/h"),
    column("tracY", 'D', "comoving position of the tracer", "Mpc/h"),
    column("tracZ", 'D', "comoving position of the tracer", "Mpc/h"),
    column("lagrX", 'D', "the random matched in this realization", "Mpc/h"),
    column("lagrY", 'D', "the random matched in this realization", "Mpc/h"),
    column("lagrZ", 'D', "the random matched in this realization", "Mpc/h"),
    column("displX", 'D', "displacement in this realization, lagr - trac", "Mpc/h"),
    column("displY", 'D', "displacement in this realization, lagr - trac", "Mpc/h"),
    column("displZ", 'D', "displacement in this realization, lagr - trac", "Mpc/h"),
    column("valid", 'J', "1 valid, 0 rejected by a filter or left out"),
    column("outsideRedshiftCut", 'J', "1 if left out of the reconstruction by the redshift cut"),
    column("outsideMask", 'J', "1 if left out of the reconstruction by the mask")};

  WriteOptions options;
  options.precision = kLossless;
  options.keywords = result_keywords("displacement field", result);
  write(file, columns, rows,
        [r, n] (const std::size_t row, std::vector<double>& values, std::vector<std::int64_t>& integers) {
          const std::size_t i = row % n;
          values[0] = (double)(row / n);
          integers[1] = (std::int64_t)i;
          for (std::size_t c = 0; c < 3; ++c) {
            values[2+c] = r->tracers[3*i+c];
            values[5+c] = r->matchedRandom[3*row+c];
            values[8+c] = r->displacement[3*row+c];
          }
          values[11] = (double)r->valid[row];
          values[12] = r->outsideRedshiftCut.empty() ? 0. : (double)r->outsideRedshiftCut[i];
          values[13] = r->outsideMask.empty() ? 0. : (double)r->outsideMask[i];
        }, options);
}


// ============================================================================


void otswap::io::writeRealSpaceCatalog (const std::string& file, const RealSpaceCatalog& catalog,
                                        const std::vector<CatalogGroup>& groups)
{
  using G = CatalogGroup;
  static const char* const names[] = {"Index", "Sky", "Cartesian", "ValidRealizations", "Neighbours",
                                      "Shift", "Status"};
  const bool lightcone = catalog.geometry == Geometry::Lightcone;
  const std::vector<G> chosen = !groups.empty() ? groups
    : lightcone ? std::vector<G>{G::Sky, G::Cartesian, G::ValidRealizations, G::Neighbours, G::Status}
                : std::vector<G>{G::Cartesian, G::ValidRealizations, G::Neighbours, G::Status};
  check_groups(chosen, {G::Sky}, lightcone, names);

  const std::size_t n = catalog.nObjects;
  const RealSpaceCatalog* c = &catalog;
  std::vector<Group> out;
  for (const G g : chosen) {
    switch (g) {
    case G::Index:
      out.push_back(index_group());
      break;
    case G::Sky:
      check_size(catalog.sky.size(), 3 * n, "sky");
      out.push_back(triple(catalog.sky, "trac", kSky,
                           "right ascension and declination of the tracer, corrected redshift",
                           kDeg, true));
      break;
    case G::Cartesian:
      check_size(catalog.cartesian.size(), 3 * n, "cartesian");
      out.push_back(triple(catalog.cartesian, "trac", kXyz, "corrected comoving position", kMpc, false));
      break;
    case G::ValidRealizations:
      check_size(catalog.validRealizations.size(), n, "validRealizations");
      out.push_back({{column("nValidRec", 'J', "number of valid OT realizations of the tracer")},
                     [c] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = (double)c->validRealizations[row];
                     }});
      break;
    case G::Neighbours:
      check_size(catalog.nNeighbours.size(), n, "nNeighbours");
      check_size(catalog.nRealizationsAveraged.size(), n, "nRealizationsAveraged");
      out.push_back({{column("nNeighbours", 'J', "tracers with a valid realization averaged within 3 "
                                                 "sigma, the tracer included if valid"),
                      column("nRealizationsAveraged", 'J', "sum of their nValidRec")},
                     [c] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = (double)c->nNeighbours[row];
                       values[1] = (double)c->nRealizationsAveraged[row];
                     }});
      break;
    case G::Shift:
      check_size(catalog.shift.size(), n, "shift");
      check_size(catalog.factor.size(), n, "factor");
      out.push_back({{column("shift", 'D', "shift applied along the line of sight", "Mpc/h"),
                      column("rsdFactor", 'D', "f/(b + 3f/5)")},
                     [c] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = c->shift[row];
                       values[1] = c->factor[row];
                     }});
      break;
    case G::Status:
      check_size(catalog.status.size(), n, "status");
      out.push_back({{column("status", 'J', "0 corrected, 1 moved by its neighbours (no valid "
                                            "realization), 2 no valid tracer within 3 sigma, 3 left "
                                            "out of the reconstruction")},
                     [c] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = (double)(int)c->status[row];
                     }});
      break;
    }
  }

  WriteOptions options;
  options.precision = kPrecision;
  options.keywords = common("real-space catalogue", catalog.geometry, n);
  options.keywords.push_back({"SIGMA", number(catalog.sigma), "width of the average, Mpc/h"});
  options.keywords.push_back({"WEIGHTED", catalog.weightByRealizations ? "1" : "0",
                              "neighbours weighted by their valid realizations"});
  if (lightcone)
    options.keywords.push_back({"NEXTRAP", integer(catalog.nExtrapolated),
                                "tracers at which b(z) was extrapolated"});
  else {
    options.keywords.push_back({"AXIS", integer(catalog.axis), "line-of-sight axis"});
    options.keywords.push_back({"ZBOX", number(catalog.boxRedshift), "redshift of the box"});
    options.keywords.push_back({"BIAS", number(catalog.boxBias), "linear bias of the tracers"});
  }
  write_groups(file, out, n, options);
}


// ============================================================================


void otswap::io::writeMpsProfile (const std::string& file, const MpsProfile& profile)
{
  const std::size_t n = profile.redshift.size();
  if (n == 0) throw Error("cannot write an empty mean particle separation profile");
  check_size(profile.mps.size(), n, "the mps of the profile");
  check_size(profile.count.size(), n, "the counts of the profile");

  std::vector<Column> columns {
    {"redshift", 'D', "centre of the redshift bin", "", profile.redshift, {}},
    {"MPS", 'D', "mean particle separation in the bin", "Mpc/h", profile.mps, {}},
    {"nTracers", 'J', "tracers in the bin", "", std::vector<double>(n), {}}};
  for (std::size_t i = 0; i < n; ++i) columns[2].data[i] = (double)profile.count[i];

  WriteOptions options;
  options.precision = kPrecision;
  options.keywords = {{"PRODUCT", "mps profile", "what the file holds"},
                      {"OTSWAPV", OTSWAP_VERSION, "otswap version"},
                      {"ZMIN", number(profile.redshiftMin), "lowest redshift binned"},
                      {"ZMAX", number(profile.redshiftMax), "highest redshift binned"},
                      {"MPSREPR", number(profile.representative),
                       "representative mps, the count-weighted median, Mpc/h"}};
  write(file, columns, options);
}
