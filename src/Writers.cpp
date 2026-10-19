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
#include <utility>
#include <vector>

#include "internal.h"

namespace {

  using otswap::io::Column;
  using otswap::io::Keyword;

  /// Significant digits of the writers that are not lossless, and of the
  /// lossless one.
  constexpr int kPrecision = 9;
  constexpr int kLossless = 17;

  /// The EXTNAME of the main table of both displacement writers, and of
  /// the parts of the complete file.
  const char* const kMainExtname = "DISPLACEMENTS";
  const char* const kProfileExtname = "MPSPROF";
  const char* const kTableExtname = "DISTANCE";

  /// Radians to degrees, as skyToDegrees converts.
  constexpr double kRadToDeg = 180. / 3.14159265358979323846;

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

  const char* verbosity_name (const otswap::Verbosity v)
  {
    return v == otswap::Verbosity::Silent ? "silent" : v == otswap::Verbosity::Normal ? "normal" : "detailed";
  }

  /// The keywords of a reconstruction, shared by both displacement
  /// writers.
  std::vector<Keyword> result_keywords (const std::string& product, const otswap::Result& r)
  {
    std::vector<Keyword> k = common(product, r.geometry, r.nObjects);
    k.push_back({"NREC", integer(r.nRealizations), "realizations"});
    k.push_back({"SEED", integer(r.config.seed), "seed of the random streams"});
    k.push_back({"CONVERG", number(r.config.convergence), "convergence threshold of the swap loop"});
    k.push_back({"CELLSIZE", number(r.config.cellSize), "grid cell size, in units of the mps"});
    k.push_back({"REJCROSS", r.config.rejectCrossings ? "1" : "0", "config: reject mask crossings"});
    k.push_back({"CFMAXPIX", integer(r.config.maxUnobservedPixelsCrossed),
                 "config: unobserved pixels a crossing may cross"});
    k.push_back({"VERBOSE", verbosity_name(r.config.verbosity), "what the call reported"});
    if (r.geometry == otswap::Geometry::Box && !std::isnan(r.mps))
      k.push_back({"MPS", number(r.mps), "mean particle separation, Mpc/h"});
    if (r.box) {
      static const char* const axis[3] = {"X", "Y", "Z"};
      for (int c = 0; c < 3; ++c) {
        k.push_back({std::string("BOX") + axis[c] + "LO", number(r.box->lower[c]),
                     std::string("lower ") + (char)('x' + c) + " of the box, Mpc/h"});
        k.push_back({std::string("BOX") + axis[c] + "HI", number(r.box->upper[c]),
                     std::string("upper ") + (char)('x' + c) + " of the box, Mpc/h"});
      }
    }
    const otswap::SelectionCounts& s = r.selection;
    if (s.boxLimits) k.push_back({"BOXLIM", "1", "the box is the limits given"});
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
    const std::pair<const char*, std::pair<std::size_t, const char*>> counts[] = {
      {"NTRAC",   {s.tracers, "tracers given"}},
      {"NTRCUT",  {s.tracersOutsideRedshiftCut, "of them, outside the redshift cut"}},
      {"NTRMASK", {s.tracersOutsideMask, "of them, on an unobserved pixel"}},
      {"NTRBOTH", {s.tracersOutsideBoth, "of them, outside both"}},
      {"NTRBOX",  {s.tracersOutsideBox, "of them, outside the box limits"}},
      {"NRAND",   {s.randoms, "randoms given"}},
      {"NRNCUT",  {s.randomsOutsideRedshiftCut, "of them, outside the redshift cut"}},
      {"NRNMASK", {s.randomsOutsideMask, "of them, on an unobserved pixel"}},
      {"NRNBOTH", {s.randomsOutsideBoth, "of them, outside both"}},
      {"NRNBOX",  {s.randomsOutsideBox, "of them, outside the box limits"}},
      {"NDISPL",  {s.displacements, "valid displacements before the filters"}},
      {"NDISPCR", {s.displacementsCrossingMask, "of them, rejected as crossing the mask"}}};
    for (const auto& c : counts)
      if (c.second.first != 0) k.push_back({c.first, integer(c.second.first), c.second.second});
    return k;
  }

  /// The checks of both displacement writers on the reconstruction as a
  /// whole, which the keywords record.
  void check_result (const otswap::Result& r)
  {
    if (r.nRealizations != r.config.nRealizations)
      throw otswap::Error("the result holds nRealizations = " + std::to_string(r.nRealizations) +
                          " and config.nRealizations = " + std::to_string(r.config.nRealizations) +
                          "; they must be equal");
    if (r.selection.boxLimits && !r.box)
      throw otswap::Error("the result holds selection.boxLimits and no box; box must equal the limits");
    if (r.selection.boxLimits && (r.selection.boxLimits->lower != r.box->lower ||
                                  r.selection.boxLimits->upper != r.box->upper))
      throw otswap::Error("the result holds selection.boxLimits different from its box; they must be "
                          "equal");
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

  /// The id column, with its description, shared by the three writers.
  Column id_column ()
  {
    return column("id", 'K', "tracer ID; its input row, from 0, when no IDs were set");
  }

  /// The ID of object i: ids[i], or i when ids is empty.
  std::int64_t id_of (const std::vector<std::int64_t>& ids, const std::size_t i)
  {
    return ids.empty() ? (std::int64_t)i : ids[i];
  }

  Group id_group (const std::vector<std::int64_t>& ids)
  {
    const std::vector<std::int64_t>* source = &ids;
    return {{id_column()},
            [source] (const std::size_t row, double*, std::int64_t* integers) {
              integers[0] = id_of(*source, row);
            }};
  }

  /// Write the groups, nRows rows; in FITS, in an HDU of that extname when
  /// it is not empty.
  void write_groups (const std::string& file, const std::vector<Group>& groups,
                     const std::size_t nRows, const otswap::io::WriteOptions& options,
                     const std::string& extname = "")
  {
    std::vector<Column> columns;
    std::vector<std::size_t> offset;
    for (const Group& g : groups) {
      offset.push_back(columns.size());
      columns.insert(columns.end(), g.columns.begin(), g.columns.end());
    }
    otswap::internal::write_parts(file, {{extname, "", columns, nRows,
      [&groups, &offset] (const std::size_t row, std::vector<double>& values,
                          std::vector<std::int64_t>& integers) {
        for (std::size_t k = 0; k < groups.size(); ++k)
          groups[k].fill(row, values.data() + offset[k], integers.data() + offset[k]);
      }, options}});
  }

  template <typename G>
  void check_groups (const std::vector<G>& groups, const std::vector<G>& lightconeOnly,
                     const bool lightcone, const char* const names[],
                     const std::vector<G>& boxOnly = {})
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
      if (lightcone)
        for (const G only : boxOnly)
          if (groups[i] == only)
            throw otswap::Error(std::string("the group ") + names[(int)groups[i]] +
                                " describes a box, and this is a lightcone");
    }
  }

}


// ============================================================================


void otswap::io::writeDisplacementMean (const std::string& file, const Result& result,
                                        const std::vector<DisplacementGroup>& groups)
{
  using G = DisplacementGroup;
  static const char* const names[] = {"Id", "TracerSky", "LagrangianSky", "Tracer", "Lagrangian",
                                      "Displacement", "ValidRealizations", "Selection", "OutsideBox"};
  const bool lightcone = result.geometry == Geometry::Lightcone;
  std::vector<G> chosen = groups;
  if (chosen.empty() && lightcone)
    chosen = {G::TracerSky, G::LagrangianSky, G::Tracer, G::Lagrangian, G::Displacement,
              G::ValidRealizations, G::Selection};
  else if (chosen.empty()) {
    chosen = {G::Tracer, G::Lagrangian, G::Displacement, G::ValidRealizations};
    if (result.selection.boxLimits) chosen.push_back(G::OutsideBox);
  }
  if (groups.empty() && !result.ids.empty()) chosen.insert(chosen.begin(), G::Id);
  check_groups(chosen, {G::TracerSky, G::LagrangianSky, G::Selection}, lightcone, names,
               {G::OutsideBox});
  check_result(result);
  internal::check_flags(result);
  internal::check_ids(result.ids, result.nObjects);

  const std::size_t n = result.nObjects;
  std::vector<Group> out;
  for (const G g : chosen) {
    switch (g) {
    case G::Id:
      out.push_back(id_group(result.ids));
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
    case G::OutsideBox: {
      const Result* r = &result;
      out.push_back({{column("outsideBox", 'J', "1 if left out of the reconstruction by the box limits")},
                     [r] (const std::size_t row, double* values, std::int64_t*) {
                       values[0] = r->outsideBox.empty() ? 0. : (double)r->outsideBox[row];
                     }});
      break;
    }
    }
  }

  WriteOptions options;
  options.precision = kPrecision;
  options.keywords = result_keywords("displacements", result);
  write_groups(file, out, n, options, kMainExtname);
}


// ============================================================================


void otswap::io::writeDisplacementComplete (const std::string& file, const Result& result,
                                            const bool lightconeTables)
{
  const std::size_t n = result.nObjects;
  const std::size_t rows = (std::size_t)result.nRealizations * n;
  const bool lightcone = result.geometry == Geometry::Lightcone;
  check_result(result);
  check_size(result.tracers.size(), 3 * n, "tracers");
  if (lightcone) check_size(result.tracersSky.size(), 3 * n, "tracersSky");
  check_size(result.displacement.size(), 3 * rows, "displacement");
  check_size(result.matchedRandom.size(), 3 * rows, "matchedRandom");
  check_size(result.valid.size(), rows, "valid");
  internal::check_flags(result);
  internal::check_ids(result.ids, n);

  const MpsProfile& profile = result.mpsProfile;
  const bool writeProfile = lightconeTables && lightcone && !profile.redshift.empty();
  const bool writeTable = lightconeTables && lightcone && result.distances != nullptr;
  if (writeProfile) {
    check_size(profile.mps.size(), profile.redshift.size(), "the mps of the profile");
    check_size(profile.count.size(), profile.redshift.size(), "the counts of the profile");
  }

  std::vector<Column> columns {column("realization", 'J', "realization, from 0"), id_column()};
  if (lightcone) {
    columns.push_back(column("tracRA", 'D', "observed right ascension of the tracer, as given", "deg"));
    columns.push_back(column("tracDec", 'D', "observed declination of the tracer", "deg"));
    columns.push_back(column("tracRed", 'D', "observed redshift of the tracer"));
  }
  const std::size_t xyz = columns.size();
  for (const char* c : kXyz) columns.push_back(column(std::string("trac") + c, 'D',
                                                      "comoving position of the tracer", "Mpc/h"));
  for (const char* c : kXyz) columns.push_back(column(std::string("lagr") + c, 'D',
                                                      "the random matched in this realization", "Mpc/h"));
  for (const char* c : kXyz)
    columns.push_back(column(std::string("displ") + c, 'D', "displacement in this realization, lagr - trac",
                             "Mpc/h"));
  columns.push_back(column("valid", 'J', "1 valid, 0 rejected by a filter or left out"));
  columns.push_back(column("outsideRedshiftCut", 'J', "1 if left out of the reconstruction by the redshift cut"));
  columns.push_back(column("outsideMask", 'J', "1 if left out of the reconstruction by the mask"));
  columns.push_back(column("outsideBox", 'J', "1 if left out of the reconstruction by the box limits"));

  internal::TablePart primary;
  primary.extname = kMainExtname;
  primary.columns = std::move(columns);
  primary.nRows = rows;
  primary.options.precision = kLossless;
  primary.options.keywords = result_keywords("complete displacements", result);
  if (writeProfile) primary.options.keywords.push_back({"MPSPART", "1", "the mps(z) profile is written"});
  if (writeTable) primary.options.keywords.push_back({"DISTPART", "1", "the distance table is written"});
  const Result* r = &result;
  primary.fillRow = [r, n, lightcone, xyz] (const std::size_t row, std::vector<double>& values,
                                            std::vector<std::int64_t>& integers) {
    const std::size_t i = row % n;
    values[0] = (double)(row / n);
    integers[1] = id_of(r->ids, i);
    if (lightcone) {
      values[2] = r->tracersSky[3*i] * kRadToDeg;
      values[3] = r->tracersSky[3*i+1] * kRadToDeg;
      values[4] = r->tracersSky[3*i+2];
    }
    for (std::size_t c = 0; c < 3; ++c) {
      values[xyz+c] = r->tracers[3*i+c];
      values[xyz+3+c] = r->matchedRandom[3*row+c];
      values[xyz+6+c] = r->displacement[3*row+c];
    }
    values[xyz+9] = (double)r->valid[row];
    values[xyz+10] = r->outsideRedshiftCut.empty() ? 0. : (double)r->outsideRedshiftCut[i];
    values[xyz+11] = r->outsideMask.empty() ? 0. : (double)r->outsideMask[i];
    values[xyz+12] = r->outsideBox.empty() ? 0. : (double)r->outsideBox[i];
  };

  std::vector<internal::TablePart> parts {primary};
  if (writeProfile) {
    internal::TablePart part;
    part.extname = kProfileExtname;
    part.twin = "mps";
    part.columns = {column("redshift", 'D', "centre of the redshift bin"),
                    column("MPS", 'D', "mean particle separation in the bin", "Mpc/h"),
                    column("nTracers", 'J', "tracers in the bin")};
    part.nRows = profile.redshift.size();
    part.options.precision = kLossless;
    part.options.keywords = {{"PRODUCT", "complete displacements mps profile", "what the table holds"},
                             {"OTSWAPV", OTSWAP_VERSION, "otswap version"},
                             {"ZMIN", number(profile.redshiftMin), "lowest redshift binned"},
                             {"ZMAX", number(profile.redshiftMax), "highest redshift binned"},
                             {"MPSREPR", number(profile.representative),
                              "representative mps, the count-weighted median, Mpc/h"}};
    const MpsProfile* p = &profile;
    part.fillRow = [p] (const std::size_t row, std::vector<double>& values, std::vector<std::int64_t>&) {
      values[0] = p->redshift[row];
      values[1] = p->mps[row];
      values[2] = (double)p->count[row];
    };
    parts.push_back(part);
  }

  if (writeTable) {
    const DistanceTable* t = result.distances.get();
    const bool growth = t->hasGrowthRate();
    internal::TablePart part;
    part.extname = kTableExtname;
    part.twin = "distances";
    part.columns = {column("redshift", 'D', "sampled redshift"),
                    column("distance", 'D', "comoving distance", "Mpc/h")};
    if (growth) part.columns.push_back(column("growthRate", 'D', "linear growth rate f"));
    part.nRows = t->redshift().size();
    part.options.precision = kLossless;
    part.options.keywords = {{"PRODUCT", "complete displacements distance table", "what the table holds"},
                             {"OTSWAPV", OTSWAP_VERSION, "otswap version"},
                             {"ZMIN", number(t->minRedshift()), "lowest redshift of the table"},
                             {"ZMAX", number(t->maxRedshift()), "highest redshift of the table"}};
    part.fillRow = [t, growth] (const std::size_t row, std::vector<double>& values,
                                std::vector<std::int64_t>&) {
      values[0] = t->redshift()[row];
      values[1] = t->distance()[row];
      if (growth) values[2] = t->growthRate()[row];
    };
    parts.push_back(part);
  }

  internal::write_parts(file, parts);
}


// ============================================================================


void otswap::io::writeRealSpaceCatalog (const std::string& file, const RealSpaceCatalog& catalog,
                                        const std::vector<CatalogGroup>& groups)
{
  using G = CatalogGroup;
  static const char* const names[] = {"Id", "Sky", "Cartesian", "ValidRealizations", "Neighbours",
                                      "Shift", "Status"};
  const bool lightcone = catalog.geometry == Geometry::Lightcone;
  std::vector<G> chosen = groups;
  if (chosen.empty()) {
    chosen = lightcone ? std::vector<G>{G::Sky, G::Cartesian, G::ValidRealizations, G::Neighbours, G::Status}
                       : std::vector<G>{G::Cartesian, G::ValidRealizations, G::Neighbours, G::Status};
    if (!catalog.ids.empty()) chosen.insert(chosen.begin(), G::Id);
  }
  check_groups(chosen, {G::Sky}, lightcone, names);
  internal::check_ids(catalog.ids, catalog.nObjects);

  const std::size_t n = catalog.nObjects;
  const RealSpaceCatalog* c = &catalog;
  std::vector<Group> out;
  for (const G g : chosen) {
    switch (g) {
    case G::Id:
      out.push_back(id_group(catalog.ids));
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
      out.push_back({{column("status", 'J', "0 corrected, 1 moved by neighbours, 2 no valid neighbour, "
                                            "3 left out")},
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
