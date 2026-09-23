/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
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
#include <numeric>
#include <string>
#include <vector>

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

    profile.mps[i] = std::pow((double)profile.count[i] / volume, -1./3.);
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


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config)
{
  const unsigned seed = internal::check_config(config);

  const std::size_t nObjects = internal::check_coordinates(tracersSky, "the tracer sky array");
  const std::size_t nRandoms = internal::check_coordinates(randomsSky, "the random sky array");

  internal::check_random_supply(nRandoms, nObjects, config.nRealizations);

  const internal::MpsProfile profile =
    internal::mps_profile(tracersSky, skyAreaDeg2, nBins, distances);

  std::vector<double> tracers(3 * nObjects), randoms(3 * nRandoms);

  for (std::size_t i = 0; i < nObjects; ++i)
    internal::to_cartesian(tracersSky[3*i], tracersSky[3*i+1],
                           distances.distanceAt(tracersSky[3*i+2]),
                           tracers[3*i], tracers[3*i+1], tracers[3*i+2]);

  for (std::size_t i = 0; i < nRandoms; ++i)
    internal::to_cartesian(randomsSky[3*i], randomsSky[3*i+1],
                           distances.distanceAt(randomsSky[3*i+2]),
                           randoms[3*i], randoms[3*i+1], randoms[3*i+2]);

  std::vector<double> mps(nObjects);
  for (std::size_t i = 0; i < nObjects; ++i)
    mps[i] = internal::mps_at(profile, tracersSky[3*i+2]);

  return internal::reconstruct(tracers, randoms, mps,
                               internal::representative(profile), config, seed);
}


// ============================================================================


otswap::Result otswap::reconstructLightcone (const std::vector<double>& tracers,
                                             const std::vector<double>& randoms,
                                             const std::vector<double>& tracersSky,
                                             const std::vector<double>& randomsSky,
                                             const double skyAreaDeg2,
                                             const unsigned nBins,
                                             const DistanceTable& distances,
                                             const Config& config)
{
  const unsigned seed = internal::check_config(config);

  const std::size_t nObjects = internal::check_coordinates(tracers, "the tracer array");
  const std::size_t nRandoms = internal::check_coordinates(randoms, "the random array");
  const std::size_t nObjectsSky = internal::check_coordinates(tracersSky, "the tracer sky array");

  // randomsSky is not read by this overload: mps(z) is measured from the
  // tracers, and Result carries no sky columns. It is still checked, and
  // the check stays, because the two random arrays are required to
  // describe the same objects. A caller that has them out of step has a
  // real error, and this is the only place that can see it; the pair is
  // otherwise consumed by the sky overload, where a mismatch would be
  // silent and would corrupt the conversion.
  const std::size_t nRandomsSky = internal::check_coordinates(randomsSky, "the random sky array");

  if (nObjects != nObjectsSky)
    throw Error("the Cartesian tracer array describes " + std::to_string(nObjects) +
                " objects and the sky one " + std::to_string(nObjectsSky));

  if (nRandoms != nRandomsSky)
    throw Error("the Cartesian random array describes " + std::to_string(nRandoms) +
                " objects and the sky one " + std::to_string(nRandomsSky));

  internal::check_random_supply(nRandoms, nObjects, config.nRealizations);

  const internal::MpsProfile profile =
    internal::mps_profile(tracersSky, skyAreaDeg2, nBins, distances);

  std::vector<double> mps(nObjects);
  for (std::size_t i = 0; i < nObjects; ++i)
    mps[i] = internal::mps_at(profile, tracersSky[3*i+2]);

  return internal::reconstruct(tracers, randoms, mps,
                               internal::representative(profile), config, seed);
}
