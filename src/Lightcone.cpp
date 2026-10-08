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
 *  the two entry points that use it.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "detmath.h"
#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;

  // Square degrees on the whole sky, 4*pi*(180/pi)^2.
  const double kFullSkyDeg2 = 4. * kPi * (180. / kPi) * (180. / kPi);

  // Relative precision demanded of a bin's mean particle separation.
  constexpr double kMpsPrecision = 0.02;

  // Smallest number of tracers a redshift bin may hold.
  //
  // mps = (N/V)^(-1/3), so the relative error on mps is one third of the
  // relative error on the count, d(mps)/mps = (1/3)(dN/N). The count in a
  // bin is Poisson, dN = sqrt(N), which gives d(mps)/mps = 1/(3*sqrt(N)).
  // Requiring that to be at most kMpsPrecision gives
  // N >= 1/(9*kMpsPrecision^2), which is 277.8 at eps = 0.02, hence 278.
  constexpr unsigned kMinTracersPerBin =
    (unsigned)(1. / (9. * kMpsPrecision * kMpsPrecision)) + 1;

  const double kNaN = std::numeric_limits<double>::quiet_NaN();

  /// Sky coordinates and distance of a Cartesian position, the right
  /// ascension in [0, 2 pi): normalize_ra returns 2 pi itself for a
  /// negative angle smaller in magnitude than half its rounding step, and
  /// -0 for -0; both are folded here to +0.
  void sky_of (const double x, const double y, const double z,
               double& ra, double& dec, double& distance)
  {
    otswap::internal::to_sky(x, y, z, ra, dec, distance);
    if (ra == 0. || ra >= 2. * kPi) ra = 0.;
  }

  void check_cut (const otswap::RedshiftCut& cut)
  {
    if (std::isnan(cut.min) || std::isnan(cut.max))
      throw otswap::Error("a bound of the redshift cut is NaN");
    if (cut.min > cut.max)
      throw otswap::Error("the redshift cut is [" + std::to_string(cut.min) + ", " +
                          std::to_string(cut.max) + "]; its lower bound exceeds the upper one");
  }

  /// Both selections, evaluated on every object of a sky array: whether its
  /// redshift lies outside the closed range of the cut, and whether it falls
  /// on an unobserved pixel of the mask, if one is given. The objects kept
  /// are those flagged by neither, in increasing order.
  struct Selection {
    std::vector<std::uint8_t> outsideCut, outsideMask;
    std::vector<std::size_t> keep;
    std::size_t nOutsideCut = 0, nOutsideMask = 0, nOutsideBoth = 0;
  };

  Selection select (const std::vector<double>& sky, const otswap::RedshiftCut& cut,
                    const otswap::Mask* mask)
  {
    const std::size_t n = sky.size() / 3;
    Selection s;
    s.outsideCut.assign(n, 0);
    s.outsideMask.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const bool outsideCut = !(sky[3*i+2] >= cut.min && sky[3*i+2] <= cut.max);
      const bool outsideMask = mask != nullptr && !mask->allows(sky[3*i], sky[3*i+1]);
      s.outsideCut[i] = outsideCut ? 1 : 0;
      s.outsideMask[i] = outsideMask ? 1 : 0;
      s.nOutsideCut += outsideCut ? 1 : 0;
      s.nOutsideMask += outsideMask ? 1 : 0;
      s.nOutsideBoth += (outsideCut && outsideMask) ? 1 : 0;
      if (!outsideCut && !outsideMask) s.keep.push_back(i);
    }
    return s;
  }

  // The rows of a flat 3-column array named by keep, in that order.
  std::vector<double> rows (const std::vector<double>& a, const std::vector<std::size_t>& keep)
  {
    std::vector<double> out(3 * keep.size());
    for (std::size_t k = 0; k < keep.size(); ++k)
      for (std::size_t c = 0; c < 3; ++c) out[3*k+c] = a[3*keep[k]+c];
    return out;
  }

  /// The supply left by the selections, checked with the counts they
  /// dropped, before anything else is computed.
  void check_kept (const otswap::SelectionCounts& counts,
                   const std::size_t keptTracers, const std::size_t keptRandoms,
                   const unsigned nRealizations)
  {
    const std::size_t nTracers = counts.tracers, nRandoms = counts.randoms;
    const std::string cut = "the redshift cut [" + std::to_string(counts.redshiftCut.min) + ", " +
                            std::to_string(counts.redshiftCut.max) + "]";
    const std::string range = !counts.maskApplied ? cut
                            : counts.redshiftCutApplied ? cut + " and the mask" : "the mask";
    if (keptTracers < otswap::internal::min_objects())
      throw otswap::Error(range + " keeps " + std::to_string(keptTracers) + " of " +
                          std::to_string(nTracers) + " tracers (" +
                          std::to_string(nTracers - keptTracers) + " dropped); at least " +
                          std::to_string(otswap::internal::min_objects()) + " are required");
    const std::size_t needed = (std::size_t)nRealizations * keptTracers;
    if (keptRandoms < needed)
      throw otswap::Error(range + " keeps " + std::to_string(keptRandoms) + " of " +
                          std::to_string(nRandoms) + " randoms (" +
                          std::to_string(nRandoms - keptRandoms) + " dropped), for " +
                          std::to_string(keptTracers) + " tracers kept; " +
                          std::to_string(nRealizations) + " realizations need " +
                          std::to_string(needed));
  }

  /// The result of the kept tracers, expanded to one row per input tracer:
  /// a tracer left out has NaN displacement and matchedRandom rows, no valid
  /// realization, and its flags set.
  otswap::Result expand (const otswap::Result& kept, const Selection& selection,
                         const std::size_t nObjects)
  {
    const std::vector<std::size_t>& keep = selection.keep;
    const unsigned nRealizations = kept.nRealizations;

    otswap::Result full;
    full.nObjects = nObjects;
    full.nRealizations = nRealizations;
    full.displacement.assign(3 * (std::size_t)nRealizations * nObjects, kNaN);
    full.matchedRandom.assign(3 * (std::size_t)nRealizations * nObjects, kNaN);
    full.valid.assign((std::size_t)nRealizations * nObjects, 0);
    full.outsideRedshiftCut = selection.outsideCut;
    full.outsideMask = selection.outsideMask;
    full.filteredNside = kept.filteredNside;

    for (std::size_t k = 0; k < keep.size(); ++k) {
      const std::size_t i = keep[k];
      for (unsigned rec = 0; rec < nRealizations; ++rec) {
        const std::size_t from = (std::size_t)rec * keep.size() + k;
        const std::size_t to = (std::size_t)rec * nObjects + i;
        full.valid[to] = kept.valid[from];
        for (std::size_t c = 0; c < 3; ++c) {
          full.displacement[3*to+c] = kept.displacement[3*from+c];
          full.matchedRandom[3*to+c] = kept.matchedRandom[3*from+c];
        }
      }
    }

    otswap::internal::summarize(full);
    return full;
  }

}


// ============================================================================


otswap::internal::MpsProfile
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
  profile.count.assign(nBins, 0u);

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

  return profile;
}


// ============================================================================


double otswap::internal::mps_at (const MpsProfile& profile, const double z)
{
  const double value = profile_at(profile.redshift, profile.mps, z);

  if (!(value > 0.))
    throw Error("the mean particle separation profile, extrapolated linearly from its "
                "terminal nodes at z = " + std::to_string(profile.redshift.front()) + " and z = " +
                std::to_string(profile.redshift.back()) + ", gives " + std::to_string(value) +
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
  for (const unsigned c : profile.count) total += c;

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
    const Selection tracerSelection = select(tracersSky, cut, mask);
    const Selection randomSelection = select(randomsSky, cut, mask);

    SelectionCounts counts;
    counts.redshiftCutApplied = std::isfinite(cut.min) || std::isfinite(cut.max);
    counts.redshiftCut = cut;
    counts.maskApplied = mask != nullptr;
    counts.tracers = nObjects;
    counts.tracersOutsideRedshiftCut = tracerSelection.nOutsideCut;
    counts.tracersOutsideMask = tracerSelection.nOutsideMask;
    counts.tracersOutsideBoth = tracerSelection.nOutsideBoth;
    counts.randoms = nRandoms;
    counts.randomsOutsideRedshiftCut = randomSelection.nOutsideCut;
    counts.randomsOutsideMask = randomSelection.nOutsideMask;
    counts.randomsOutsideBoth = randomSelection.nOutsideBoth;

    const std::vector<std::size_t>& keepTracers = tracerSelection.keep;
    const std::vector<std::size_t>& keepRandoms = randomSelection.keep;
    const bool dropTracers = keepTracers.size() != nObjects;
    const bool dropRandoms = keepRandoms.size() != nRandoms;
    if (dropTracers || dropRandoms)
      check_kept(counts, keepTracers.size(), keepRandoms.size(), config.nRealizations);

    const std::vector<double> keptTracersSky = dropTracers ? rows(tracersSky, keepTracers) : std::vector<double>();
    const std::vector<double>& ts = dropTracers ? keptTracersSky : tracersSky;
    const std::size_t nKept = ts.size() / 3;

    internal::check_random_supply(keepRandoms.size(), nKept, config.nRealizations);

    const internal::MpsProfile profile = internal::mps_profile(ts, skyAreaDeg2, nBins, distances);

    std::vector<double> convertedTracers, convertedRandoms, keptTracers, keptRandoms;
    if (tracers == nullptr) {
      convertedTracers = toCartesian(ts, distances);
      convertedRandoms = dropRandoms ? toCartesian(rows(randomsSky, keepRandoms), distances)
                                     : toCartesian(randomsSky, distances);
    }
    else {
      if (dropTracers) keptTracers = rows(*tracers, keepTracers);
      if (dropRandoms) keptRandoms = rows(*randoms, keepRandoms);
    }
    const std::vector<double>& t = tracers == nullptr ? convertedTracers
                                 : dropTracers ? keptTracers : *tracers;
    const std::vector<double>& r = tracers == nullptr ? convertedRandoms
                                 : dropRandoms ? keptRandoms : *randoms;

    std::vector<double> mps(nKept);
    for (std::size_t i = 0; i < nKept; ++i)
      mps[i] = internal::mps_at(profile, ts[3*i+2]);

    Result result = internal::reconstruct(t, r, mps, internal::representative(profile), config, seed);
    if (dropTracers) result = expand(result, tracerSelection, nObjects);

    if (mask != nullptr && config.rejectCrossings) {
      const auto valid = [&result] {
        return (std::size_t)std::count(result.valid.begin(), result.valid.end(), 1);
      };
      counts.crossingsRejected = true;
      counts.maxUnobservedPixelsCrossed = config.maxUnobservedPixelsCrossed;
      counts.displacements = valid();
      rejectMaskCrossings(result, *mask, config.maxUnobservedPixelsCrossed);
      counts.displacementsCrossingMask = counts.displacements - valid();
    }

    result.lagrangianSky.assign(3 * nObjects, kNaN);
    for (std::size_t k = 0; k < nKept; ++k) {
      const std::size_t i = dropTracers ? keepTracers[k] : k;
      if (result.validRealizations[i] == 0) continue;
      double* sky = &result.lagrangianSky[3*i];
      double distance = 0.;
      sky_of(t[3*k] + result.meanDisplacement[3*i], t[3*k+1] + result.meanDisplacement[3*i+1],
             t[3*k+2] + result.meanDisplacement[3*i+2], sky[0], sky[1], distance);
      try {
        sky[2] = distances.redshiftAt(distance);
      }
      catch (const Error&) {}
    }

    result.selection = counts;

    if (config.verbose) {
      const std::string message = counts.message();
      if (!message.empty()) std::clog << message << std::flush;
    }

    return result;
  }

}


// ============================================================================


std::string otswap::SelectionCounts::message () const
{
  std::ostringstream out;
  out.imbue(std::locale::classic());

  auto line = [&] (const char* what, const std::size_t given, const std::size_t outsideCut,
                   const std::size_t outsideMask, const std::size_t outsideBoth) {
    out << "otswap: kept " << given - (outsideCut + outsideMask - outsideBoth) << " of " << given
        << " " << what << ": ";
    if (redshiftCutApplied)
      out << outsideCut << " outside the redshift cut [" << redshiftCut.min << ", "
          << redshiftCut.max << "]";
    if (redshiftCutApplied && maskApplied) out << ", ";
    if (maskApplied) out << outsideMask << " outside the mask";
    if (redshiftCutApplied && maskApplied) out << " (" << outsideBoth << " outside both)";
    out << '\n';
  };

  if (redshiftCutApplied || maskApplied) {
    line("tracers", tracers, tracersOutsideRedshiftCut, tracersOutsideMask, tracersOutsideBoth);
    line("randoms", randoms, randomsOutsideRedshiftCut, randomsOutsideMask, randomsOutsideBoth);
  }

  if (crossingsRejected)
    out << "otswap: rejected " << displacementsCrossingMask << " of " << displacements
        << " displacements crossing more than " << maxUnobservedPixelsCrossed
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
