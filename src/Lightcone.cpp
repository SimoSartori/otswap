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
#include <limits>
#include <numeric>
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

  void check_cut (const otswap::RedshiftCut& cut)
  {
    if (std::isnan(cut.min) || std::isnan(cut.max))
      throw otswap::Error("a bound of the redshift cut is NaN");
    if (cut.min > cut.max)
      throw otswap::Error("the redshift cut is [" + std::to_string(cut.min) + ", " +
                          std::to_string(cut.max) + "]; its lower bound exceeds the upper one");
  }

  // Indices, increasing, of the objects of a sky array whose redshift lies
  // in the closed range of the cut.
  std::vector<std::size_t> inside (const std::vector<double>& sky, const otswap::RedshiftCut& cut)
  {
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < sky.size() / 3; ++i)
      if (sky[3*i+2] >= cut.min && sky[3*i+2] <= cut.max) keep.push_back(i);
    return keep;
  }

  // The rows of a flat 3-column array named by keep, in that order.
  std::vector<double> rows (const std::vector<double>& a, const std::vector<std::size_t>& keep)
  {
    std::vector<double> out(3 * keep.size());
    for (std::size_t k = 0; k < keep.size(); ++k)
      for (std::size_t c = 0; c < 3; ++c) out[3*k+c] = a[3*keep[k]+c];
    return out;
  }

  // The supply left by the cut, checked with the counts it dropped, before
  // anything else is computed.
  void check_kept (const otswap::RedshiftCut& cut,
                   const std::size_t keptTracers, const std::size_t nTracers,
                   const std::size_t keptRandoms, const std::size_t nRandoms,
                   const unsigned nRealizations)
  {
    const std::string range = "the redshift cut [" + std::to_string(cut.min) + ", " +
                              std::to_string(cut.max) + "]";
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

  // The result of the kept tracers, expanded to one row per input tracer:
  // a cut tracer has NaN displacement and matchedRandom rows, no valid
  // realization, and its flag set.
  otswap::Result expand (const otswap::Result& kept, const std::vector<std::size_t>& keep,
                         const std::size_t nObjects)
  {
    const unsigned nRealizations = kept.nRealizations;

    otswap::Result full;
    full.nObjects = nObjects;
    full.nRealizations = nRealizations;
    full.displacement.assign(3 * (std::size_t)nRealizations * nObjects, kNaN);
    full.matchedRandom.assign(3 * (std::size_t)nRealizations * nObjects, kNaN);
    full.valid.assign((std::size_t)nRealizations * nObjects, 0);
    full.outsideRedshiftCut.assign(nObjects, 1);
    full.filteredNside = kept.filteredNside;

    for (std::size_t k = 0; k < keep.size(); ++k) {
      const std::size_t i = keep[k];
      full.outsideRedshiftCut[i] = 0;
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
  const std::size_t nObjects = internal::check_coordinates(sky, "the sky array", true);

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


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config,
                                             const RedshiftCut& cut)
{
  const unsigned seed = internal::check_config(config);

  const std::size_t nObjects = internal::check_coordinates(tracersSky, "the tracer sky array");
  const std::size_t nRandoms = internal::check_coordinates(randomsSky, "the random sky array");

  // The cut comes first: everything below sees the kept objects only.
  check_cut(cut);
  const std::vector<std::size_t> keepTracers = inside(tracersSky, cut);
  const std::vector<std::size_t> keepRandoms = inside(randomsSky, cut);
  const bool cutTracers = keepTracers.size() != nObjects;
  const bool cutRandoms = keepRandoms.size() != nRandoms;
  if (cutTracers || cutRandoms)
    check_kept(cut, keepTracers.size(), nObjects, keepRandoms.size(), nRandoms, config.nRealizations);

  const std::vector<double> keptTracersSky = cutTracers ? rows(tracersSky, keepTracers) : std::vector<double>();
  const std::vector<double> keptRandomsSky = cutRandoms ? rows(randomsSky, keepRandoms) : std::vector<double>();
  const std::vector<double>& ts = cutTracers ? keptTracersSky : tracersSky;
  const std::vector<double>& rs = cutRandoms ? keptRandomsSky : randomsSky;
  const std::size_t nKept = ts.size() / 3;

  internal::check_random_supply(rs.size() / 3, nKept, config.nRealizations);

  const internal::MpsProfile profile = internal::mps_profile(ts, skyAreaDeg2, nBins, distances);

  const std::vector<double> tracers = toCartesian(ts, distances);
  const std::vector<double> randoms = toCartesian(rs, distances);

  std::vector<double> mps(nKept);
  for (std::size_t i = 0; i < nKept; ++i)
    mps[i] = internal::mps_at(profile, ts[3*i+2]);

  const Result result = internal::reconstruct(tracers, randoms, mps,
                                              internal::representative(profile), config, seed);
  return cutTracers ? expand(result, keepTracers, nObjects) : result;
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
  const unsigned seed = internal::check_config(config);

  const std::size_t nObjects = internal::check_coordinates(tracers, "the tracer array");
  const std::size_t nRandoms = internal::check_coordinates(randoms, "the random array");
  const std::size_t nObjectsSky = internal::check_coordinates(tracersSky, "the tracer sky array");

  // randomsSky is not read by this overload: mps(z) is measured from the
  // tracers, and Result carries no sky columns. It is validated all the
  // same, and its object count must equal the Cartesian random array's,
  // because the two arrays are required to describe the same randoms; a
  // caller whose arrays are out of step gets an error here.
  const std::size_t nRandomsSky = internal::check_coordinates(randomsSky, "the random sky array");

  if (nObjects != nObjectsSky)
    throw Error("the Cartesian tracer array describes " + std::to_string(nObjects) +
                " objects and the sky one " + std::to_string(nObjectsSky));

  if (nRandoms != nRandomsSky)
    throw Error("the Cartesian random array describes " + std::to_string(nRandoms) +
                " objects and the sky one " + std::to_string(nRandomsSky));

  // The cut comes first, decided on the sky redshifts; the same rows are
  // dropped from the Cartesian arrays.
  check_cut(cut);
  const std::vector<std::size_t> keepTracers = inside(tracersSky, cut);
  const std::vector<std::size_t> keepRandoms = inside(randomsSky, cut);
  const bool cutTracers = keepTracers.size() != nObjects;
  const bool cutRandoms = keepRandoms.size() != nRandoms;
  if (cutTracers || cutRandoms)
    check_kept(cut, keepTracers.size(), nObjects, keepRandoms.size(), nRandoms, config.nRealizations);

  const std::vector<double> keptTracers = cutTracers ? rows(tracers, keepTracers) : std::vector<double>();
  const std::vector<double> keptTracersSky = cutTracers ? rows(tracersSky, keepTracers) : std::vector<double>();
  const std::vector<double> keptRandoms = cutRandoms ? rows(randoms, keepRandoms) : std::vector<double>();
  const std::vector<double>& t = cutTracers ? keptTracers : tracers;
  const std::vector<double>& ts = cutTracers ? keptTracersSky : tracersSky;
  const std::vector<double>& r = cutRandoms ? keptRandoms : randoms;
  const std::size_t nKept = t.size() / 3;

  internal::check_random_supply(r.size() / 3, nKept, config.nRealizations);

  const internal::MpsProfile profile = internal::mps_profile(ts, skyAreaDeg2, nBins, distances);

  std::vector<double> mps(nKept);
  for (std::size_t i = 0; i < nKept; ++i)
    mps[i] = internal::mps_at(profile, ts[3*i+2]);

  const Result result = internal::reconstruct(t, r, mps,
                                              internal::representative(profile), config, seed);
  return cutTracers ? expand(result, keepTracers, nObjects) : result;
}
