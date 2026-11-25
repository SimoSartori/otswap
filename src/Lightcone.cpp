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
 *  @file src/Lightcone.cpp
 *
 *  @brief Lightcone geometry: the mean particle separation profile and
 *  the entry points that use it.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "detmath.h"
#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;

  /// pi/180 and 180/pi as doubles: the factors numpy.deg2rad and
  /// numpy.rad2deg multiply by.
  constexpr double kDegToRad = kPi / 180.;
  constexpr double kRadToDeg = 180. / kPi;

  /// Square degrees on the whole sky, 4*pi*(180/pi)^2.
  const double kFullSkyDeg2 = 4. * kPi * (180. / kPi) * (180. / kPi);

  /// Relative precision demanded of a bin's mean particle separation.
  constexpr double kMpsPrecision = 0.02;

  /// Smallest number of tracers a redshift bin may hold.
  ///
  /// mps = (N/V)^(-1/3), so the relative error on mps is one third of the
  /// relative error on the count, d(mps)/mps = (1/3)(dN/N). The count in a
  /// bin is Poisson, dN = sqrt(N), which gives d(mps)/mps = 1/(3*sqrt(N)).
  /// Requiring that to be at most kMpsPrecision gives
  /// N >= 1/(9*kMpsPrecision^2), which is 277.8 at kMpsPrecision = 0.02, hence 278.
  constexpr unsigned kMinTracersPerBin =
    (unsigned)(1. / (9. * kMpsPrecision * kMpsPrecision)) + 1;

  const double kNaN = std::numeric_limits<double>::quiet_NaN();

  /// Sky coordinates and distance of a Cartesian position, the right
  /// ascension folded into [0, 2 pi) by fold_ra.
  void sky_of (const double x, const double y, const double z,
               double& ra, double& dec, double& distance)
  {
    otswap::internal::to_sky(x, y, z, ra, dec, distance);
    ra = otswap::internal::fold_ra(ra);
  }

  void check_cut (const otswap::RedshiftCut& cut)
  {
    if (std::isnan(cut.min) || std::isnan(cut.max))
      throw otswap::Error("a bound of the redshift cut is NaN");
    if (cut.min > cut.max)
      throw otswap::Error("the redshift cut is [" + std::to_string(cut.min) + ", " +
                          std::to_string(cut.max) + "]; its lower bound exceeds the upper one");
  }

  void check_cut (const otswap::DistanceCut& cut)
  {
    if (std::isnan(cut.min) || std::isnan(cut.max))
      throw otswap::Error("a bound of the distance cut is NaN");
    if (cut.min > cut.max)
      throw otswap::Error("the distance cut is [" + std::to_string(cut.min) + ", " +
                          std::to_string(cut.max) + "] in coordinate units; its lower bound exceeds "
                          "the upper one");
  }

  /// One line of the selection report, for the tracers or the randoms:
  /// outsideCut counts the objects outside the redshift cut or, for a
  /// distance cut, outside the distance cut.
  void selection_line (std::ostream& out, const otswap::SelectionCounts& counts, const char* what,
                       const std::size_t given, const std::size_t outsideCut,
                       const std::size_t outsideMask, const std::size_t outsideBoth)
  {
    const bool cut = counts.redshiftCut || counts.distanceCut;
    out << "otswap: kept " << given - (outsideCut + outsideMask - outsideBoth) << " of " << given
        << " " << what << ": ";
    if (counts.redshiftCut)
      out << outsideCut << " outside the redshift cut [" << counts.redshiftCut->min << ", "
          << counts.redshiftCut->max << "]";
    else if (counts.distanceCut)
      out << outsideCut << " outside the distance cut [" << counts.distanceCut->min << ", "
          << counts.distanceCut->max << "] in coordinate units";
    if (cut && counts.maskApplied) out << ", ";
    if (counts.maskApplied) out << outsideMask << " outside the mask";
    if (cut && counts.maskApplied) out << " (" << outsideBoth << " outside both)";
    out << '\n';
  }

  /// The line of a Detailed report that describes an mps(z) or mps(d)
  /// profile.
  void profile_line (std::ostream& out, const otswap::MpsProfile& p)
  {
    const auto range = std::minmax_element(p.mps.begin(), p.mps.end());
    if (p.distance.empty())
      out << "otswap: mps(z) from " << *range.first << " to " << *range.second << " in coordinate units, in "
          << p.mps.size() << (p.mps.size() == 1 ? " bin" : " bins") << " over ["
          << p.redshiftMin << ", " << p.redshiftMax << "]\n";
    else
      out << "otswap: mps(d) from " << *range.first << " to " << *range.second << " in coordinate units, in "
          << p.mps.size() << (p.mps.size() == 1 ? " bin" : " bins") << " over ["
          << p.distanceMin << ", " << p.distanceMax << "] in coordinate units\n";
  }

  /// The profile of mps_profile and mps_profile_distance: the values, one
  /// per tracer at stride @p stride, binned uniformly over their range; a
  /// bin edge is a redshift turned into a distance by the table, or, with
  /// none, a distance itself.
  otswap::MpsProfile binned_profile (const double* values, const std::size_t stride,
                                     const std::size_t nObjects, const double skyAreaDeg2,
                                     const unsigned nBins, const otswap::DistanceTable* distances)
  {
    using otswap::Error;
    const std::string what = distances != nullptr ? "redshift" : "distance";

    if (nBins == 0)
      throw Error("the number of " + what + " bins is 0; it must be at least 1");

    if (!std::isfinite(skyAreaDeg2) || skyAreaDeg2 <= 0. || skyAreaDeg2 > kFullSkyDeg2)
      throw Error("the sky area is " + std::to_string(skyAreaDeg2) +
                  " square degrees; it must be positive and at most the whole sky, " +
                  std::to_string(kFullSkyDeg2));

    double zMin = values[0], zMax = values[0];
    for (std::size_t i = 0; i < nObjects; ++i) {
      zMin = std::min(zMin, values[stride*i]);
      zMax = std::max(zMax, values[stride*i]);
    }

    if (!(zMax > zMin))
      throw Error("every tracer is at " + what + " " + std::to_string(zMin) +
                  "; the mean particle separation cannot be measured in " + what + " bins");

    if (distances != nullptr && (zMin < distances->minRedshift() || zMax > distances->maxRedshift()))
      throw Error("the tracer redshifts span [" + std::to_string(zMin) + ", " +
                  std::to_string(zMax) + "], outside the distance table's [" +
                  std::to_string(distances->minRedshift()) + ", " +
                  std::to_string(distances->maxRedshift()) + "]");

    const double step = (zMax - zMin) / (double)nBins;

    otswap::MpsProfile profile;
    std::vector<double>& centre = distances != nullptr ? profile.redshift : profile.distance;
    centre.resize(nBins);
    profile.mps.resize(nBins);
    profile.count.assign(nBins, 0);
    (distances != nullptr ? profile.redshiftMin : profile.distanceMin) = zMin;
    (distances != nullptr ? profile.redshiftMax : profile.distanceMax) = zMax;

    for (std::size_t i = 0; i < nObjects; ++i) {
      int bin = (int)((values[stride*i] - zMin) / step);
      if (bin >= (int)nBins) bin = (int)nBins - 1;
      if (bin < 0) bin = 0;
      ++profile.count[(std::size_t)bin];
    }

    const unsigned supported =
      (unsigned)std::max<std::size_t>(1, nObjects / kMinTracersPerBin);

    for (unsigned i = 0; i < nBins; ++i)
      if (profile.count[i] < kMinTracersPerBin)
        throw Error(what + " bin " + std::to_string(i) + " (from 0) of " + std::to_string(nBins) +
                    ", covering [" + std::to_string(zMin + i*step) + ", " +
                    std::to_string(zMin + (i+1)*step) + "], holds " +
                    std::to_string(profile.count[i]) + " tracers, below the " +
                    std::to_string(kMinTracersPerBin) + " needed for its density to be "
                    "meaningful to " + std::to_string(kMpsPrecision * 100.) +
                    " per cent; this catalog of " + std::to_string(nObjects) +
                    " tracers supports at most " + std::to_string(supported) + " bins");

    const double skyFraction = skyAreaDeg2 / kFullSkyDeg2;

    for (unsigned i = 0; i < nBins; ++i) {
      const double lo = zMin + i * step;
      const double hi = (i + 1 == nBins) ? zMax : zMin + (i+1) * step;

      centre[i] = 0.5 * (lo + hi);

      const double dLo = distances != nullptr ? distances->distanceAt(lo) : lo;
      const double dHi = distances != nullptr ? distances->distanceAt(hi) : hi;
      const double volume = skyFraction * (4. * kPi / 3.) * (dHi*dHi*dHi - dLo*dLo*dLo);

      if (!(volume > 0.))
        throw Error(what + " bin " + std::to_string(i) + " (from 0) has a non-positive shell volume");

      profile.mps[i] = otswap::internal::det_pow((double)profile.count[i] / volume, -1./3.);
    }

    profile.representative = otswap::internal::representative(profile);
    return profile;
  }

}


// ============================================================================


otswap::MpsProfile
otswap::internal::mps_profile (const std::vector<double>& tracersSky,
                               const double skyAreaDeg2,
                               const unsigned nBins,
                               const DistanceTable& distances)
{
  return binned_profile(tracersSky.data() + 2, 3, tracersSky.size() / 3, skyAreaDeg2, nBins, &distances);
}


// ============================================================================


otswap::MpsProfile
otswap::internal::mps_profile_distance (const std::vector<double>& distance,
                                        const double skyAreaDeg2,
                                        const unsigned nBins)
{
  return binned_profile(distance.data(), 1, distance.size(), skyAreaDeg2, nBins, nullptr);
}


// ============================================================================


double otswap::MpsProfile::at (const double x) const
{
  const bool inDistance = !distance.empty();
  const std::vector<double>& nodes = inDistance ? distance : redshift;

  if (!inDistance && (redshift.empty() || redshift.size() != mps.size()))
    throw Error("the mean particle separation profile is empty or malformed: it holds " +
                std::to_string(redshift.size()) + " redshifts and " + std::to_string(mps.size()) +
                " values");
  if (inDistance && (!redshift.empty() || distance.size() != mps.size()))
    throw Error("the mean particle separation profile is malformed: it holds " +
                std::to_string(redshift.size()) + " redshifts, " + std::to_string(distance.size()) +
                " distances and " + std::to_string(mps.size()) + " values; it must hold redshifts "
                "or distances, as many as values");

  const double value = internal::profile_at(nodes, mps, x);

  if (!(value > 0.)) {
    const std::string name = inDistance ? "d" : "z";
    throw Error("the mean particle separation profile, extrapolated linearly from its "
                "terminal nodes at " + name + " = " + std::to_string(nodes.front()) + " and " + name +
                " = " + std::to_string(nodes.back()) + ", gives " + std::to_string(value) +
                " at " + name + " = " + std::to_string(x) + "; it must be positive. Use fewer " +
                (inDistance ? "distance" : "redshift") + " bins, so that the terminal bins are less noisy");
  }

  return value;
}


// ============================================================================


double otswap::internal::representative (const MpsProfile& profile)
{
  std::vector<std::size_t> order(profile.mps.size());
  std::iota(order.begin(), order.end(), (std::size_t)0);
  std::sort(order.begin(), order.end(),
            [&profile] (const std::size_t a, const std::size_t b)
            { return profile.mps[a] < profile.mps[b]; });

  std::size_t total = 0;
  for (const std::size_t c : profile.count) total += c;

  const double half = 0.5 * (double)total;
  std::size_t running = 0;

  for (const std::size_t i : order) {
    running += profile.count[i];
    if ((double)running >= half) return profile.mps[i];
  }

  return profile.mps[order.back()];
}


// ============================================================================


std::vector<double> otswap::toCartesian (const std::vector<double>& sky,
                                         const DistanceTable& distances)
{
  const std::size_t nObjects = internal::check_sky(sky, "the sky array", true);

  std::vector<double> cartesian(3 * nObjects);

  for (std::size_t i = 0; i < nObjects; ++i) {
    double distance = 0.;
    try {
      distance = distances.distanceAt(sky[3*i+2]);
    }
    catch (const Error& e) {
      throw Error("object " + std::to_string(i) + " (from 0) of the sky array: " + e.what());
    }
    internal::to_cartesian(sky[3*i], sky[3*i+1], distance,
                           cartesian[3*i], cartesian[3*i+1], cartesian[3*i+2]);
  }

  return cartesian;
}


// ============================================================================


std::vector<double> otswap::toSky (const std::vector<double>& cartesian,
                                   const DistanceTable& distances)
{
  if (cartesian.size() % 3 != 0)
    throw Error("the Cartesian array holds " + std::to_string(cartesian.size()) +
                " entries, which is not a multiple of three");

  const std::size_t nObjects = cartesian.size() / 3;
  std::vector<double> sky(3 * nObjects, kNaN);

  for (std::size_t i = 0; i < nObjects; ++i) {
    const double* p = &cartesian[3*i];
    for (std::size_t c = 0; c < 3; ++c)
      if (std::isinf(p[c]))
        throw Error("object " + std::to_string(i) + " (from 0) of the Cartesian array: component " +
                    std::to_string(c) + " (from 0) is infinite");
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) continue;

    double distance = 0.;
    sky_of(p[0], p[1], p[2], sky[3*i], sky[3*i+1], distance);
    try {
      sky[3*i+2] = distances.redshiftAt(distance);
    }
    catch (const Error& e) {
      throw Error("object " + std::to_string(i) + " (from 0) of the Cartesian array: " + e.what() +
                  "; use a distance table covering a wider redshift range than [" +
                  std::to_string(distances.minRedshift()) + ", " +
                  std::to_string(distances.maxRedshift()) + "]");
    }
  }

  return sky;
}


// ============================================================================


void otswap::internal::sky_of_positions (const std::vector<double>& cartesian,
                                         const DistanceTable* distances,
                                         std::vector<double>& sky,
                                         std::vector<double>* distance)
{
  sky.resize(cartesian.size());
  std::fill(sky.begin(), sky.end(), kNaN);
  if (distance != nullptr) distance->assign(cartesian.size() / 3, kNaN);

  for (std::size_t i = 0; 3*i < cartesian.size(); ++i) {
    const double* p = &cartesian[3*i];
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) continue;
    double d = 0.;
    sky_of(p[0], p[1], p[2], sky[3*i], sky[3*i+1], d);
    if (distance != nullptr) (*distance)[i] = d;
    if (distances == nullptr) continue;
    try {
      sky[3*i+2] = distances->redshiftAt(d);
    }
    catch (const Error&) {}
  }
}


// ============================================================================


void otswap::pointToSky (const double x, const double y, const double z,
                         double& distance, double& ra, double& dec) noexcept
{
  if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
    distance = ra = dec = kNaN;
    return;
  }
  sky_of(x, y, z, ra, dec, distance);
}


// ============================================================================


std::vector<double> otswap::skyToRadians (std::vector<double> skyDegrees)
{
  if (skyDegrees.size() % 3 != 0)
    throw Error("the sky array holds " + std::to_string(skyDegrees.size()) +
                " entries, which is not a multiple of three");

  for (std::size_t i = 0; i < skyDegrees.size(); i += 3)
    if (std::isfinite(skyDegrees[i+1]) && std::fabs(skyDegrees[i+1]) > 90.)
      throw Error("the sky array holds a declination of " + std::to_string(skyDegrees[i+1]) +
                  " degrees at object " + std::to_string(i / 3) + " (from 0), outside [-90, 90]");

  for (std::size_t i = 0; i < skyDegrees.size(); i += 3) {
    skyDegrees[i]   *= kDegToRad;
    skyDegrees[i+1] *= kDegToRad;
  }
  return skyDegrees;
}


// ============================================================================


std::vector<double> otswap::skyToDegrees (std::vector<double> skyRadians)
{
  if (skyRadians.size() % 3 != 0)
    throw Error("the sky array holds " + std::to_string(skyRadians.size()) +
                " entries, which is not a multiple of three");

  for (std::size_t i = 0; i < skyRadians.size(); i += 3) {
    skyRadians[i]   *= kRadToDeg;
    skyRadians[i+1] *= kRadToDeg;
    if (skyRadians[i] >= 360.) skyRadians[i] = 0.;
  }
  return skyRadians;
}


// ============================================================================


namespace {

  /// Every key of the selection kept: the selection of a Cartesian
  /// overload without a mask.
  otswap::internal::Selection keep_all (const std::size_t n)
  {
    otswap::internal::Selection s;
    s.outsideCut.assign(n, 0);
    s.outsideMask.assign(n, 0);
    s.outsideBox.assign(n, 0);
    s.keep.resize(n);
    std::iota(s.keep.begin(), s.keep.end(), (std::size_t)0);
    return s;
  }

  /// Right ascension, declination and distance of each position, by
  /// pointToSky: the array the Cartesian overloads select on, the distance
  /// standing where a sky array has its redshift.
  std::vector<double> sky_and_distance (const std::vector<double>& positions)
  {
    std::vector<double> out(positions.size());
    for (std::size_t i = 0; 3*i < positions.size(); ++i)
      otswap::pointToSky(positions[3*i], positions[3*i+1], positions[3*i+2],
                         out[3*i+2], out[3*i], out[3*i+1]);
    return out;
  }

  /// Whether a cut has a finite bound, and so is applied and reported.
  template <typename Cut>
  bool finite (const Cut& cut)
  {
    return std::isfinite(cut.min) || std::isfinite(cut.max);
  }

  /// The selection of a Cartesian overload: the distance cut and the mask
  /// (null for none) on the right ascension, declination and distance of
  /// each position, by pointToSky; every object kept when neither applies.
  otswap::internal::Selection select_cartesian (const std::vector<double>& positions,
                                                const otswap::DistanceCut& cut,
                                                const otswap::Mask* mask)
  {
    if (mask == nullptr && !finite(cut)) return keep_all(positions.size() / 3);
    return otswap::internal::select(sky_and_distance(positions), otswap::RedshiftCut{cut.min, cut.max}, mask);
  }

  /// The tracers kept of the Cartesian overloads and of their
  /// measureMpsProfile: the array checked, its sky and distance, and the
  /// selection by the distance cut and the mask (null for the sky area).
  struct CartesianTracers {
    std::size_t nObjects = 0;
    std::vector<double> sky;          ///< ra, dec, distance per tracer
    otswap::internal::Selection selection;
    otswap::SelectionCounts counts;
    std::vector<double> keptDistance; ///< distance per tracer kept
  };

  CartesianTracers cartesian_tracers (const std::vector<double>& tracers, const otswap::DistanceCut& cut,
                                      const otswap::Mask* mask)
  {
    CartesianTracers c;
    c.nObjects = otswap::internal::check_coordinates(tracers, "the tracer array");
    check_cut(cut);
    c.sky = sky_and_distance(tracers);
    c.selection = mask == nullptr && !finite(cut)
      ? keep_all(c.nObjects)
      : otswap::internal::select(c.sky, otswap::RedshiftCut{cut.min, cut.max}, mask);
    c.counts = otswap::internal::tracer_counts(otswap::RedshiftCut{}, mask, c.nObjects, c.selection);
    if (finite(cut)) c.counts.distanceCut = cut;
    c.counts.tracersOutsideRedshiftCut = 0;
    c.counts.tracersOutsideDistanceCut = c.selection.nOutsideCut;
    c.keptDistance.resize(c.selection.keep.size());
    for (std::size_t k = 0; k < c.keptDistance.size(); ++k)
      c.keptDistance[k] = c.sky[3*c.selection.keep[k]+2];
    return c;
  }

  /// The Cartesian overloads, the mask null for the sky area. The objects
  /// are selected on the sky directions pointToSky gives; the profile is in
  /// distance; the table, when given, only fills the redshifts of the
  /// result.
  otswap::Result cartesian_lightcone (const std::vector<double>& tracers,
                                      const std::vector<double>& randoms,
                                      const double skyAreaDeg2,
                                      const otswap::Mask* mask,
                                      const unsigned nBins,
                                      const otswap::DistanceTable* distances,
                                      const otswap::Config& config,
                                      const otswap::DistanceCut& cut)
  {
    using namespace otswap;

    const auto start = std::chrono::steady_clock::now();
    const unsigned seed = internal::check_config(config);

    CartesianTracers c = cartesian_tracers(tracers, cut, mask);
    const std::size_t nObjects = c.nObjects;
    const std::size_t nRandoms = internal::check_coordinates(randoms, "the random array");

    if (nObjects > (std::size_t)std::numeric_limits<unsigned>::max())
      throw Error("the catalog holds " + std::to_string(nObjects) + " objects; it must hold at most " +
                  std::to_string(std::numeric_limits<unsigned>::max()));

    const internal::Selection randomSelection = select_cartesian(randoms, cut, mask);

    SelectionCounts counts = c.counts;
    counts.randoms = nRandoms;
    counts.randomsOutsideDistanceCut = randomSelection.nOutsideCut;
    counts.randomsOutsideMask = randomSelection.nOutsideMask;
    counts.randomsOutsideBoth = randomSelection.nOutsideBoth;

    const std::vector<std::size_t>& keepTracers = c.selection.keep;
    const std::vector<std::size_t>& keepRandoms = randomSelection.keep;
    const bool dropTracers = keepTracers.size() != nObjects;
    const bool dropRandoms = keepRandoms.size() != nRandoms;
    if (dropTracers || dropRandoms)
      internal::check_kept(counts, keepTracers.size(), keepRandoms.size(), config.nRealizations);

    const std::size_t nKept = keepTracers.size();
    internal::check_random_supply(keepRandoms.size(), nKept, config.nRealizations);

    const bool reporting = config.log != nullptr && config.verbosity != Verbosity::Silent;
    const bool detailed = reporting && config.verbosity == Verbosity::Detailed;
    const bool rejecting = mask != nullptr && config.rejectCrossings;
    const std::string selected = reporting ? counts.message() : std::string();
    if (reporting) internal::write_lines(*config.log, selected);

    MpsProfile profile = internal::mps_profile_distance(c.keptDistance, skyAreaDeg2, nBins);

    std::ostringstream profileLine;
    profileLine.imbue(std::locale::classic());
    if (detailed) profile_line(profileLine, profile);
    if (!rejecting && detailed) internal::write_lines(*config.log, profileLine.str());

    const std::vector<double> keptTracers = dropTracers ? internal::rows(tracers, keepTracers) : std::vector<double>();
    const std::vector<double> keptRandoms = dropRandoms ? internal::rows(randoms, keepRandoms) : std::vector<double>();
    const std::vector<double>& t = dropTracers ? keptTracers : tracers;
    const std::vector<double>& r = dropRandoms ? keptRandoms : randoms;

    std::vector<double> mps(nKept);
    for (std::size_t i = 0; i < nKept; ++i)
      mps[i] = profile.at(c.keptDistance[i]);

    Result result = internal::reconstruct(t, r, mps, profile.representative, config, seed,
                                          nullptr, nullptr, false);
    if (dropTracers) result = internal::expand(result, c.selection, nObjects);
    result.outsideDistanceCut.swap(result.outsideRedshiftCut);

    result.tracers = tracers;
    result.tracersSky = std::move(c.sky);
    for (std::size_t i = 0; i < nObjects; ++i) {
      double& third = result.tracersSky[3*i+2];
      double z = kNaN;
      if (distances != nullptr) {
        try {
          z = distances->redshiftAt(third);
        }
        catch (const Error&) {}
      }
      third = z;
    }
    result.geometry = Geometry::Lightcone;
    result.mpsProfile = std::move(profile);
    if (distances != nullptr) result.distances = std::make_shared<const DistanceTable>(*distances);
    result.selection = counts;
    internal::summarize(result);

    if (rejecting) {
      rejectMaskCrossings(result, *mask, config.maxUnobservedPixelsCrossed);
      if (reporting) internal::write_lines(*config.log, result.selection.message().substr(selected.size()));
      if (detailed) internal::write_lines(*config.log, profileLine.str());
    }

    result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (reporting)
      internal::write_lines(*config.log, internal::reconstruction_line("reconstructLightcone", result));

    return result;
  }

  /// The sky overloads, the mask null for those with the sky area. The
  /// selections come first, decided on the sky coordinates, so everything
  /// after them sees the kept objects only.
  otswap::Result lightcone (const std::vector<double>& tracersSky,
                            const std::vector<double>& randomsSky,
                            const double skyAreaDeg2,
                            const otswap::Mask* mask,
                            const unsigned nBins,
                            const otswap::DistanceTable& distances,
                            const otswap::Config& config,
                            const otswap::RedshiftCut& cut)
  {
    using namespace otswap;

    const auto start = std::chrono::steady_clock::now();
    const unsigned seed = internal::check_config(config);

    const std::size_t nObjects = internal::check_sky(tracersSky, "the tracer sky array");
    const std::size_t nRandoms = internal::check_sky(randomsSky, "the random sky array");

    if (nObjects > (std::size_t)std::numeric_limits<unsigned>::max())
      throw Error("the catalog holds " + std::to_string(nObjects) + " objects; it must hold at most " +
                  std::to_string(std::numeric_limits<unsigned>::max()));

    check_cut(cut);
    const internal::Selection tracerSelection = internal::select(tracersSky, cut, mask);
    const internal::Selection randomSelection = internal::select(randomsSky, cut, mask);

    SelectionCounts counts = internal::tracer_counts(cut, mask, nObjects, tracerSelection);
    counts.randoms = nRandoms;
    counts.randomsOutsideRedshiftCut = randomSelection.nOutsideCut;
    counts.randomsOutsideMask = randomSelection.nOutsideMask;
    counts.randomsOutsideBoth = randomSelection.nOutsideBoth;

    const std::vector<std::size_t>& keepTracers = tracerSelection.keep;
    const std::vector<std::size_t>& keepRandoms = randomSelection.keep;
    const bool dropTracers = keepTracers.size() != nObjects;
    const bool dropRandoms = keepRandoms.size() != nRandoms;
    if (dropTracers || dropRandoms)
      internal::check_kept(counts, keepTracers.size(), keepRandoms.size(), config.nRealizations);

    const std::vector<double> keptTracersSky = dropTracers ? internal::rows(tracersSky, keepTracers) : std::vector<double>();
    const std::vector<double>& ts = dropTracers ? keptTracersSky : tracersSky;
    const std::size_t nKept = ts.size() / 3;

    internal::check_random_supply(keepRandoms.size(), nKept, config.nRealizations);

    const bool reporting = config.log != nullptr && config.verbosity != Verbosity::Silent;
    const bool detailed = reporting && config.verbosity == Verbosity::Detailed;
    const bool rejecting = mask != nullptr && config.rejectCrossings;
    const std::string selected = reporting ? counts.message() : std::string();
    if (reporting) internal::write_lines(*config.log, selected);

    MpsProfile profile = internal::mps_profile(ts, skyAreaDeg2, nBins, distances);

    std::ostringstream profileLine;
    profileLine.imbue(std::locale::classic());
    if (detailed) profile_line(profileLine, profile);
    if (!rejecting && detailed) internal::write_lines(*config.log, profileLine.str());

    std::vector<double> convertedTracers = toCartesian(ts, distances);
    const std::vector<double> convertedRandoms =
      dropRandoms ? toCartesian(internal::rows(randomsSky, keepRandoms), distances)
                  : toCartesian(randomsSky, distances);
    const std::vector<double>& t = convertedTracers;
    const std::vector<double>& r = convertedRandoms;

    std::vector<double> mps(nKept);
    for (std::size_t i = 0; i < nKept; ++i)
      mps[i] = profile.at(ts[3*i+2]);

    Result result = internal::reconstruct(t, r, mps, profile.representative, config, seed,
                                          nullptr, nullptr, false);
    if (dropTracers) result = internal::expand(result, tracerSelection, nObjects);

    if (!dropTracers) result.tracers = std::move(convertedTracers);
    else {
      result.tracers.assign(3 * nObjects, kNaN);
      std::size_t next = 0;
      for (std::size_t i = 0; i < nObjects; ++i) {
        double* row = &result.tracers[3*i];
        if (next < nKept && keepTracers[next] == i) {
          for (std::size_t c = 0; c < 3; ++c) row[c] = convertedTracers[3*next+c];
          ++next;
          continue;
        }
        try {
          const double distance = distances.distanceAt(tracersSky[3*i+2]);
          internal::to_cartesian(tracersSky[3*i], tracersSky[3*i+1], distance, row[0], row[1], row[2]);
        }
        catch (const Error&) {}
      }
    }
    result.tracersSky = tracersSky;
    result.geometry = Geometry::Lightcone;
    result.mpsProfile = std::move(profile);
    result.distances = std::make_shared<const DistanceTable>(distances);
    result.selection = counts;
    internal::summarize(result);

    if (rejecting) {
      rejectMaskCrossings(result, *mask, config.maxUnobservedPixelsCrossed);
      if (reporting) internal::write_lines(*config.log, result.selection.message().substr(selected.size()));
      if (detailed) internal::write_lines(*config.log, profileLine.str());
    }

    result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (reporting)
      internal::write_lines(*config.log, internal::reconstruction_line("reconstructLightcone", result));

    return result;
  }

  /// The profile of measureMpsProfile: the steps of lightcone that concern
  /// the tracers, in its order, with the mask null for the overload with
  /// the sky area. The kept-tracer minimum applies whenever a selection
  /// does, since the randoms that decide it in lightcone are not given.
  otswap::MpsProfile tracer_profile (const std::vector<double>& tracersSky,
                                     const double skyAreaDeg2,
                                     const otswap::Mask* mask,
                                     const unsigned nBins,
                                     const otswap::DistanceTable& distances,
                                     const otswap::RedshiftCut& cut,
                                     const otswap::Verbosity verbosity,
                                     std::ostream* log)
  {
    using namespace otswap;

    const std::size_t nObjects = internal::check_sky(tracersSky, "the tracer sky array");

    check_cut(cut);
    const internal::Selection selection = internal::select(tracersSky, cut, mask);
    const SelectionCounts counts = internal::tracer_counts(cut, mask, nObjects, selection);

    const std::vector<std::size_t>& keep = selection.keep;
    const bool selected = counts.redshiftCut || counts.maskApplied;
    if (selected) internal::check_kept_tracers(counts, keep.size());

    const bool reporting = log != nullptr && verbosity != Verbosity::Silent;
    if (reporting && selected) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      selection_line(out, counts, "tracers", counts.tracers, counts.tracersOutsideRedshiftCut,
                     counts.tracersOutsideMask, counts.tracersOutsideBoth);
      internal::write_lines(*log, out.str());
    }

    const bool dropTracers = keep.size() != nObjects;
    const std::vector<double> keptSky = dropTracers ? internal::rows(tracersSky, keep) : std::vector<double>();
    MpsProfile result = internal::mps_profile(dropTracers ? keptSky : tracersSky, skyAreaDeg2, nBins,
                                              distances);

    if (reporting) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      if (verbosity == Verbosity::Detailed) profile_line(out, result);
      out << "otswap: measureMpsProfile: " << result.mps.size()
          << (result.mps.size() == 1 ? " bin" : " bins") << " over [" << result.redshiftMin << ", "
          << result.redshiftMax << "]\n";
      internal::write_lines(*log, out.str());
    }

    return result;
  }

  /// The profile of the Cartesian measureMpsProfile: the steps of
  /// cartesian_lightcone that concern the tracers, in its order, the mask
  /// null for the sky area.
  otswap::MpsProfile cartesian_profile (const std::vector<double>& tracers,
                                        const double skyAreaDeg2,
                                        const otswap::Mask* mask,
                                        const unsigned nBins,
                                        const otswap::DistanceCut& cut,
                                        const otswap::Verbosity verbosity,
                                        std::ostream* log)
  {
    using namespace otswap;

    const CartesianTracers c = cartesian_tracers(tracers, cut, mask);
    const bool selected = c.counts.distanceCut || c.counts.maskApplied;
    if (selected) internal::check_kept_tracers(c.counts, c.selection.keep.size());

    const bool reporting = log != nullptr && verbosity != Verbosity::Silent;
    if (reporting && selected) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      selection_line(out, c.counts, "tracers", c.counts.tracers, c.counts.tracersOutsideDistanceCut,
                     c.counts.tracersOutsideMask, c.counts.tracersOutsideBoth);
      internal::write_lines(*log, out.str());
    }

    MpsProfile result = internal::mps_profile_distance(c.keptDistance, skyAreaDeg2, nBins);

    if (reporting) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      if (verbosity == Verbosity::Detailed) profile_line(out, result);
      out << "otswap: measureMpsProfile: " << result.mps.size()
          << (result.mps.size() == 1 ? " bin" : " bins") << " over [" << result.distanceMin << ", "
          << result.distanceMax << "] in coordinate units\n";
      internal::write_lines(*log, out.str());
    }

    return result;
  }

}


namespace {

  /// One line of the selection report of a box with limits, for the
  /// tracers or the randoms.
  void box_line (std::ostream& out, const otswap::BoxLimits& box, const char* what,
                 const std::size_t given, const std::size_t outside)
  {
    out << "otswap: kept " << given - outside << " of " << given << " " << what << ": " << outside
        << " outside the box";
    for (std::size_t k = 0; k < 3; ++k)
      out << (k ? " x [" : " [") << box.lower[k] << ", " << box.upper[k] << "]";
    out << '\n';
  }

  /// One line of the selection report of a periodic box, for the tracers
  /// or the randoms.
  void wrap_line (std::ostream& out, const otswap::BoxLimits& box, const char* what,
                  const std::size_t given, const std::size_t wrapped)
  {
    out << "otswap: wrapped " << wrapped << " of " << given << " " << what
        << " into the periodic box";
    for (std::size_t k = 0; k < 3; ++k)
      out << (k ? " x [" : " [") << box.lower[k] << ", " << box.upper[k] << "]";
    out << '\n';
  }

}


// ============================================================================


std::string otswap::SelectionCounts::message () const
{
  std::ostringstream out;
  out.imbue(std::locale::classic());

  if (boxLimits && periodic) {
    wrap_line(out, *boxLimits, "tracers", tracers, tracersWrapped);
    if (randoms > 0) wrap_line(out, *boxLimits, "randoms", randoms, randomsWrapped);
  }
  else if (boxLimits) {
    box_line(out, *boxLimits, "tracers", tracers, tracersOutsideBox);
    if (randoms > 0) box_line(out, *boxLimits, "randoms", randoms, randomsOutsideBox);
  }

  if (redshiftCut || distanceCut || maskApplied) {
    selection_line(out, *this, "tracers", tracers,
                   redshiftCut ? tracersOutsideRedshiftCut : tracersOutsideDistanceCut, tracersOutsideMask,
                   tracersOutsideBoth);
    selection_line(out, *this, "randoms", randoms,
                   redshiftCut ? randomsOutsideRedshiftCut : randomsOutsideDistanceCut, randomsOutsideMask,
                   randomsOutsideBoth);
  }

  if (maxUnobservedPixelsCrossed)
    out << "otswap: rejected " << displacementsCrossingMask << " of " << displacements
        << " displacements crossing more than " << *maxUnobservedPixelsCrossed
        << " unobserved pixels in a row\n";

  return out.str();
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config,
                                             const RedshiftCut& cut)
{
  return lightcone(tracersSky, randomsSky, skyAreaDeg2, nullptr, nBins, distances, config, cut);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracers,
                                             const std::vector<double>& randoms,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const Config& config,
                                             const DistanceTable* distances,
                                             const DistanceCut& cut)
{
  return cartesian_lightcone(tracers, randoms, skyAreaDeg2, nullptr, nBins, distances, config, cut);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const Mask& mask,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config,
                                             const RedshiftCut& cut)
{
  return lightcone(tracersSky, randomsSky, mask.skyAreaDeg2(), &mask, nBins, distances, config, cut);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracers,
                                             const std::vector<double>& randoms,
                                             const Mask& mask,
                                             const unsigned nBins,
                                             const Config& config,
                                             const DistanceTable* distances,
                                             const DistanceCut& cut)
{
  return cartesian_lightcone(tracers, randoms, mask.skyAreaDeg2(), &mask, nBins, distances, config, cut);
}


// ============================================================================


otswap::MpsProfile otswap::measureMpsProfile (const std::vector<double>& tracersSky,
                                              const double skyAreaDeg2,
                                              const unsigned nBins,
                                              const DistanceTable& distances,
                                              const RedshiftCut& cut,
                                              const Verbosity verbosity,
                                              std::ostream* log)
{
  return tracer_profile(tracersSky, skyAreaDeg2, nullptr, nBins, distances, cut, verbosity, log);
}


// ============================================================================


otswap::MpsProfile otswap::measureMpsProfile (const std::vector<double>& tracersSky,
                                              const Mask& mask,
                                              const unsigned nBins,
                                              const DistanceTable& distances,
                                              const RedshiftCut& cut,
                                              const Verbosity verbosity,
                                              std::ostream* log)
{
  return tracer_profile(tracersSky, mask.skyAreaDeg2(), &mask, nBins, distances, cut, verbosity, log);
}


// ============================================================================


otswap::MpsProfile otswap::measureMpsProfile (const std::vector<double>& tracers,
                                              const double skyAreaDeg2,
                                              const unsigned nBins,
                                              const DistanceCut& cut,
                                              const Verbosity verbosity,
                                              std::ostream* log)
{
  return cartesian_profile(tracers, skyAreaDeg2, nullptr, nBins, cut, verbosity, log);
}


// ============================================================================


otswap::MpsProfile otswap::measureMpsProfile (const std::vector<double>& tracers,
                                              const Mask& mask,
                                              const unsigned nBins,
                                              const DistanceCut& cut,
                                              const Verbosity verbosity,
                                              std::ostream* log)
{
  return cartesian_profile(tracers, mask.skyAreaDeg2(), &mask, nBins, cut, verbosity, log);
}
