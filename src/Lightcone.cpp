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
  /// N >= 1/(9*kMpsPrecision^2), which is 277.8 at eps = 0.02, hence 278.
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

  /// One line of the selection report, for the tracers or the randoms.
  void selection_line (std::ostream& out, const otswap::SelectionCounts& counts, const char* what,
                       const std::size_t given, const std::size_t outsideCut,
                       const std::size_t outsideMask, const std::size_t outsideBoth)
  {
    out << "otswap: kept " << given - (outsideCut + outsideMask - outsideBoth) << " of " << given
        << " " << what << ": ";
    if (counts.redshiftCut)
      out << outsideCut << " outside the redshift cut [" << counts.redshiftCut->min << ", "
          << counts.redshiftCut->max << "]";
    if (counts.redshiftCut && counts.maskApplied) out << ", ";
    if (counts.maskApplied) out << outsideMask << " outside the mask";
    if (counts.redshiftCut && counts.maskApplied) out << " (" << outsideBoth << " outside both)";
    out << '\n';
  }

  /// The line of a Detailed report that describes an mps(z) profile.
  void profile_line (std::ostream& out, const otswap::MpsProfile& p)
  {
    const auto range = std::minmax_element(p.mps.begin(), p.mps.end());
    out << "otswap: mps(z) from " << *range.first << " to " << *range.second << " Mpc/h in "
        << p.mps.size() << (p.mps.size() == 1 ? " bin" : " bins") << " over ["
        << p.redshiftMin << ", " << p.redshiftMax << "]\n";
  }

}


// ============================================================================


otswap::MpsProfile
otswap::internal::mps_profile (const std::vector<double>& tracersSky,
                               const double skyAreaDeg2,
                               const unsigned nBins,
                               const DistanceTable& distances)
{
  const std::size_t nObjects = tracersSky.size() / 3;

  if (nBins == 0)
    throw Error("nBins is zero; at least one redshift bin is required");

  if (!std::isfinite(skyAreaDeg2) || skyAreaDeg2 <= 0. || skyAreaDeg2 > kFullSkyDeg2)
    throw Error("the sky area is " + std::to_string(skyAreaDeg2) +
                " square degrees; it must be positive and at most the whole sky, " +
                std::to_string(kFullSkyDeg2));

  double zMin = tracersSky[2], zMax = tracersSky[2];
  for (std::size_t i = 0; i < nObjects; ++i) {
    zMin = std::min(zMin, tracersSky[3*i+2]);
    zMax = std::max(zMax, tracersSky[3*i+2]);
  }

  if (!(zMax > zMin))
    throw Error("every tracer is at redshift " + std::to_string(zMin) +
                "; the mean particle separation cannot be measured in redshift bins");

  if (zMin < distances.minRedshift() || zMax > distances.maxRedshift())
    throw Error("the tracer redshifts span [" + std::to_string(zMin) + ", " +
                std::to_string(zMax) + "], outside the distance table's [" +
                std::to_string(distances.minRedshift()) + ", " +
                std::to_string(distances.maxRedshift()) + "]");

  const double step = (zMax - zMin) / (double)nBins;

  MpsProfile profile;
  profile.redshift.resize(nBins);
  profile.mps.resize(nBins);
  profile.count.assign(nBins, 0);
  profile.redshiftMin = zMin;
  profile.redshiftMax = zMax;

  for (std::size_t i = 0; i < nObjects; ++i) {
    int bin = (int)((tracersSky[3*i+2] - zMin) / step);
    if (bin >= (int)nBins) bin = (int)nBins - 1;
    if (bin < 0) bin = 0;
    ++profile.count[(std::size_t)bin];
  }

  const unsigned supported =
    (unsigned)std::max<std::size_t>(1, nObjects / kMinTracersPerBin);

  for (unsigned i = 0; i < nBins; ++i)
    if (profile.count[i] < kMinTracersPerBin)
      throw Error("redshift bin " + std::to_string(i) + " of " + std::to_string(nBins) +
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

    profile.redshift[i] = 0.5 * (lo + hi);

    const double dLo = distances.distanceAt(lo);
    const double dHi = distances.distanceAt(hi);
    const double volume = skyFraction * (4. * kPi / 3.) * (dHi*dHi*dHi - dLo*dLo*dLo);

    if (!(volume > 0.))
      throw Error("redshift bin " + std::to_string(i) + " has a non-positive shell volume");

    profile.mps[i] = det_pow((double)profile.count[i] / volume, -1./3.);
  }

  profile.representative = representative(profile);
  return profile;
}


// ============================================================================


double otswap::MpsProfile::at (const double z) const
{
  if (redshift.empty() || redshift.size() != mps.size())
    throw Error("the mean particle separation profile is empty or malformed: it holds " +
                std::to_string(redshift.size()) + " redshifts and " + std::to_string(mps.size()) +
                " values");

  const double value = internal::profile_at(redshift, mps, z);

  if (!(value > 0.))
    throw Error("the mean particle separation profile, extrapolated linearly from its "
                "terminal nodes at z = " + std::to_string(redshift.front()) + " and z = " +
                std::to_string(redshift.back()) + ", gives " + std::to_string(value) +
                " at z = " + std::to_string(z) + "; it must be positive. Use fewer redshift "
                "bins, so that the terminal bins are less noisy");

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
      throw Error("object " + std::to_string(i) + " of the sky array: " + e.what());
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
        throw Error("object " + std::to_string(i) + " of the Cartesian array: component " +
                    std::to_string(c) + " is infinite");
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) continue;

    double distance = 0.;
    sky_of(p[0], p[1], p[2], sky[3*i], sky[3*i+1], distance);
    try {
      sky[3*i+2] = distances.redshiftAt(distance);
    }
    catch (const Error& e) {
      throw Error("object " + std::to_string(i) + " of the Cartesian array: " + e.what() +
                  "; use a distance table covering a wider redshift range than [" +
                  std::to_string(distances.minRedshift()) + ", " +
                  std::to_string(distances.maxRedshift()) + "]");
    }
  }

  return sky;
}


// ============================================================================


void otswap::internal::sky_of_positions (const std::vector<double>& cartesian,
                                         const DistanceTable& distances,
                                         std::vector<double>& sky)
{
  sky.resize(cartesian.size());
  std::fill(sky.begin(), sky.end(), kNaN);

  for (std::size_t i = 0; 3*i < cartesian.size(); ++i) {
    const double* p = &cartesian[3*i];
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) continue;
    double distance = 0.;
    sky_of(p[0], p[1], p[2], sky[3*i], sky[3*i+1], distance);
    try {
      sky[3*i+2] = distances.redshiftAt(distance);
    }
    catch (const Error&) {}
  }
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
                  " degrees at object " + std::to_string(i / 3) + ", outside [-90, 90]");

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

  /// Every lightcone overload: the Cartesian arrays are null for the sky
  /// overloads, and the mask null for those with the sky area. The
  /// selections come first, decided on the sky coordinates, so everything
  /// after them sees the kept objects only; with Cartesian arrays the same
  /// rows are dropped from them, which is why the sky random array must
  /// describe the same randoms.
  otswap::Result lightcone (const std::vector<double>* tracers,
                            const std::vector<double>* randoms,
                            const std::vector<double>& tracersSky,
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

    std::size_t nObjects = 0, nRandoms = 0;
    if (tracers == nullptr) {
      nObjects = internal::check_sky(tracersSky, "the tracer sky array");
      nRandoms = internal::check_sky(randomsSky, "the random sky array");
    }
    else {
      nObjects = internal::check_coordinates(*tracers, "the tracer array");
      nRandoms = internal::check_coordinates(*randoms, "the random array");
      const std::size_t nObjectsSky = internal::check_sky(tracersSky, "the tracer sky array");
      const std::size_t nRandomsSky = internal::check_sky(randomsSky, "the random sky array");

      if (nObjects != nObjectsSky)
        throw Error("the Cartesian tracer array describes " + std::to_string(nObjects) +
                    " objects and the sky one " + std::to_string(nObjectsSky));

      if (nRandoms != nRandomsSky)
        throw Error("the Cartesian random array describes " + std::to_string(nRandoms) +
                    " objects and the sky one " + std::to_string(nRandomsSky));
    }

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

    MpsProfile profile = internal::mps_profile(ts, skyAreaDeg2, nBins, distances);

    std::vector<double> convertedTracers, convertedRandoms, keptTracers, keptRandoms;
    if (tracers == nullptr) {
      convertedTracers = toCartesian(ts, distances);
      convertedRandoms = dropRandoms ? toCartesian(internal::rows(randomsSky, keepRandoms), distances)
                                     : toCartesian(randomsSky, distances);
    }
    else {
      if (dropTracers) keptTracers = internal::rows(*tracers, keepTracers);
      if (dropRandoms) keptRandoms = internal::rows(*randoms, keepRandoms);
    }
    const std::vector<double>& t = tracers == nullptr ? convertedTracers
                                 : dropTracers ? keptTracers : *tracers;
    const std::vector<double>& r = tracers == nullptr ? convertedRandoms
                                 : dropRandoms ? keptRandoms : *randoms;

    std::vector<double> mps(nKept);
    for (std::size_t i = 0; i < nKept; ++i)
      mps[i] = profile.at(ts[3*i+2]);

    Result result = internal::reconstruct(t, r, mps, profile.representative, config, seed);
    if (dropTracers) result = internal::expand(result, tracerSelection, nObjects);

    if (tracers != nullptr) result.tracers = *tracers;
    else if (!dropTracers) result.tracers = std::move(convertedTracers);
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

    if (mask != nullptr && config.rejectCrossings)
      rejectMaskCrossings(result, *mask, config.maxUnobservedPixelsCrossed);

    result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (config.log != nullptr && config.verbosity != Verbosity::Silent) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      out << result.selection.message();
      if (config.verbosity == Verbosity::Detailed) profile_line(out, result.mpsProfile);
      out << internal::reconstruction_line("reconstructLightcone", result);
      *config.log << out.str() << std::flush;
    }

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

    const bool dropTracers = keep.size() != nObjects;
    const std::vector<double> keptSky = dropTracers ? internal::rows(tracersSky, keep) : std::vector<double>();
    MpsProfile result = internal::mps_profile(dropTracers ? keptSky : tracersSky, skyAreaDeg2, nBins,
                                              distances);

    if (log != nullptr && verbosity != Verbosity::Silent) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      if (selected)
        selection_line(out, counts, "tracers", counts.tracers, counts.tracersOutsideRedshiftCut,
                       counts.tracersOutsideMask, counts.tracersOutsideBoth);
      if (verbosity == Verbosity::Detailed) profile_line(out, result);
      out << "otswap: measureMpsProfile: " << result.mps.size()
          << (result.mps.size() == 1 ? " bin" : " bins") << " over [" << result.redshiftMin << ", "
          << result.redshiftMax << "]\n";
      *log << out.str() << std::flush;
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

}


// ============================================================================


std::string otswap::SelectionCounts::message () const
{
  std::ostringstream out;
  out.imbue(std::locale::classic());

  if (boxLimits) {
    box_line(out, *boxLimits, "tracers", tracers, tracersOutsideBox);
    if (randoms > 0) box_line(out, *boxLimits, "randoms", randoms, randomsOutsideBox);
  }

  if (redshiftCut || maskApplied) {
    selection_line(out, *this, "tracers", tracers, tracersOutsideRedshiftCut, tracersOutsideMask,
                   tracersOutsideBoth);
    selection_line(out, *this, "randoms", randoms, randomsOutsideRedshiftCut, randomsOutsideMask,
                   randomsOutsideBoth);
  }

  if (maxUnobservedPixelsCrossed)
    out << "otswap: rejected " << displacementsCrossingMask << " of " << displacements
        << " displacements crossing more than " << *maxUnobservedPixelsCrossed
        << " unobserved pixels\n";

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
  return lightcone(nullptr, nullptr, tracersSky, randomsSky, skyAreaDeg2, nullptr, nBins,
                   distances, config, cut);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracers,
                                             const std::vector<double>& randoms,
                                             const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config,
                                             const RedshiftCut& cut)
{
  return lightcone(&tracers, &randoms, tracersSky, randomsSky, skyAreaDeg2, nullptr, nBins,
                   distances, config, cut);
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
  return lightcone(nullptr, nullptr, tracersSky, randomsSky, mask.skyAreaDeg2(), &mask, nBins,
                   distances, config, cut);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracers,
                                             const std::vector<double>& randoms,
                                             const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const Mask& mask,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config,
                                             const RedshiftCut& cut)
{
  return lightcone(&tracers, &randoms, tracersSky, randomsSky, mask.skyAreaDeg2(), &mask, nBins,
                   distances, config, cut);
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
