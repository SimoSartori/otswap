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
 *  @file src/DistanceTable.cpp
 *
 *  @brief The sampled redshift to comoving distance relation.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <limits>

#include "internal.h"

namespace {

  // Speed of light in km/s. With H0 = 100 h km/s/Mpc, the Hubble distance
  // c/H0 is (c/100) / h Mpc, which is c/100 Mpc/h. The table is in Mpc/h,
  // so h cancels out of the comoving distance and only enters if a caller
  // converts to Mpc. h is validated regardless: a non-positive h is a
  // caller error whatever the unit.
  constexpr double kSpeedOfLight = 299792.458;

  constexpr double kHubbleDistance = kSpeedOfLight / 100.;

  // The growth equation is started deep in matter domination, where
  // D is proportional to a to well below the accuracy anything here
  // needs. The start is also kept well below the first sampled point, so
  // that a table reaching to high redshift still has a run-up.
  constexpr double kGrowthStartRedshift = 1000.;
  constexpr double kGrowthStartMargin = 100.;

  // Runge-Kutta steps for the run-up, and per sampled interval after it.
  constexpr unsigned kGrowthLeadSteps = 4000;
  constexpr unsigned kGrowthSteps = 8;

  // Sub-intervals of Simpson's rule per sampled step. The integrand
  // 1/E(z) is smooth and monotonic, so the error falls as the fourth
  // power of the step; four sub-intervals on a step of a few times 1e-3
  // in z puts the truncation error far below the linear interpolation
  // error the table itself carries.
  constexpr unsigned kSimpsonSteps = 4;

  double one_over_E (const double z, const double OmegaM,
                     const double w0, const double wa)
  {
    const double opz = 1. + z;
    const double matter = OmegaM * opz * opz * opz;

    const double exponent = 3. * (1. + w0 + wa) * std::log1p(z) - 3. * wa * z / opz;
    if (!std::isfinite(exponent))
      throw otswap::Error("the dark-energy density exponent is not finite at z = " +
                          std::to_string(z));

    const double darkEnergy = (1. - OmegaM) * std::exp(exponent);
    const double E2 = matter + darkEnergy;

    if (!std::isfinite(E2) || E2 <= 0.)
      throw otswap::Error("the given parameters make E^2(z) non-positive or non-finite at z = " +
                          std::to_string(z));

    return 1. / std::sqrt(E2);
  }


  // The dark-energy density relative to its value today, in the scale
  // factor. It is the CPL form written with 1+z = 1/a.
  double f_DE (const double a, const double w0, const double wa)
  {
    return std::pow(a, -3. * (1. + w0 + wa)) * std::exp(-3. * wa * (1. - a));
  }


  double E2_of_a (const double a, const double OmegaM, const double w0, const double wa)
  {
    return OmegaM / (a*a*a) + (1. - OmegaM) * f_DE(a, w0, wa);
  }


  double dE2_da (const double a, const double OmegaM, const double w0, const double wa)
  {
    const double dfDE = f_DE(a, w0, wa) * (-3. * (1. + w0 + wa) / a + 3. * wa);
    return -3. * OmegaM / (a*a*a*a) + (1. - OmegaM) * dfDE;
  }


  // The linear growth equation, as a first-order system in the scale
  // factor:
  //
  //   D'' + (3/a + E'/E) D' - (3/2) OmegaM / (a^5 E^2) D = 0
  //
  // with ' = d/da and E'/E = (E^2)' / (2 E^2). The growth rate follows
  // from the solution as f = dlnD/dlna = a D'/D, so the normalisation of
  // D never enters and no fitting form is involved.
  void growth_derivatives (const double a, const double D, const double dD,
                           const double OmegaM, const double w0, const double wa,
                           double& outD, double& outdD)
  {
    const double E2 = E2_of_a(a, OmegaM, w0, wa);

    if (!std::isfinite(E2) || E2 <= 0.)
      throw otswap::Error("the given parameters make E^2 non-positive or non-finite at a = " +
                          std::to_string(a));

    const double damping = 3./a + dE2_da(a, OmegaM, w0, wa) / (2. * E2);
    const double source = 1.5 * OmegaM / (a*a*a*a*a * E2);

    outD = dD;
    outdD = -damping * dD + source * D;
  }


  void growth_advance (double& a, double& D, double& dD, const double target,
                       const unsigned steps, const double OmegaM,
                       const double w0, const double wa)
  {
    const double h = (target - a) / (double)steps;

    for (unsigned n = 0; n < steps; ++n) {
      double k1D, k1dD, k2D, k2dD, k3D, k3dD, k4D, k4dD;

      growth_derivatives(a, D, dD, OmegaM, w0, wa, k1D, k1dD);
      growth_derivatives(a + 0.5*h, D + 0.5*h*k1D, dD + 0.5*h*k1dD,
                         OmegaM, w0, wa, k2D, k2dD);
      growth_derivatives(a + 0.5*h, D + 0.5*h*k2D, dD + 0.5*h*k2dD,
                         OmegaM, w0, wa, k3D, k3dD);
      growth_derivatives(a + h, D + h*k3D, dD + h*k3dD, OmegaM, w0, wa, k4D, k4dD);

      D  += (h/6.) * (k1D  + 2.*k2D  + 2.*k3D  + k4D);
      dD += (h/6.) * (k1dD + 2.*k2dD + 2.*k3dD + k4dD);
      a  += h;
    }

    a = target;
  }

}


// ============================================================================


otswap::DistanceTable::DistanceTable (const double OmegaM, const double h,
                                      const double w0, const double wa,
                                      const double zMin, const double zMax,
                                      const unsigned nSamples)
{
  if (!std::isfinite(OmegaM) || OmegaM <= 0. || OmegaM > 1.)
    throw Error("OmegaM is " + std::to_string(OmegaM) + "; it must lie in (0, 1]");

  if (!std::isfinite(h) || h <= 0.)
    throw Error("h is " + std::to_string(h) + "; it must be positive");

  if (!std::isfinite(w0) || !std::isfinite(wa))
    throw Error("w0 and wa must both be finite");

  if (!std::isfinite(zMin) || zMin < 0.)
    throw Error("zMin is " + std::to_string(zMin) + "; it must be finite and non-negative");

  if (!std::isfinite(zMax) || zMax <= zMin)
    throw Error("zMax is " + std::to_string(zMax) + "; it must be finite and exceed zMin " +
                std::to_string(zMin));

  if (nSamples < 2)
    throw Error("nSamples is " + std::to_string(nSamples) + "; at least two are required");

  m_redshift.resize(nSamples);
  m_distance.resize(nSamples);

  for (unsigned i = 0; i < nSamples; ++i)
    m_redshift[i] = zMin + (zMax - zMin) * (double)i / (double)(nSamples - 1);
  m_redshift.front() = zMin;
  m_redshift.back() = zMax;

  double integral = 0.;

  if (zMin > 0.) {
    const unsigned lead = kSimpsonSteps * nSamples;
    const double step = zMin / (double)lead;
    double sum = one_over_E(0., OmegaM, w0, wa) + one_over_E(zMin, OmegaM, w0, wa);
    for (unsigned k = 1; k < lead; ++k)
      sum += (k % 2 == 1 ? 4. : 2.) * one_over_E(k * step, OmegaM, w0, wa);
    integral = sum * step / 3.;
  }

  m_distance[0] = kHubbleDistance * integral;

  for (unsigned i = 1; i < nSamples; ++i) {
    const double lo = m_redshift[i-1], hi = m_redshift[i];
    const double step = (hi - lo) / (double)kSimpsonSteps;
    double sum = one_over_E(lo, OmegaM, w0, wa) + one_over_E(hi, OmegaM, w0, wa);
    for (unsigned k = 1; k < kSimpsonSteps; ++k)
      sum += (k % 2 == 1 ? 4. : 2.) * one_over_E(lo + k * step, OmegaM, w0, wa);
    integral += sum * step / 3.;
    m_distance[i] = kHubbleDistance * integral;
  }

  for (unsigned i = 1; i < nSamples; ++i)
    if (!(m_distance[i] > m_distance[i-1]))
      throw Error("the sampled comoving distance is not strictly increasing at z = " +
                  std::to_string(m_redshift[i]) + "; the sampling is too fine for double precision");

  m_growthRate.resize(nSamples);

  // The samples ascend in redshift, so they descend in scale factor. The
  // integration runs the other way, from the earliest time forward, and
  // fills the array from its last entry back.
  const double aFirst = 1. / (1. + m_redshift.back());
  const double aStart = std::min(1. / (1. + kGrowthStartRedshift),
                                 aFirst / kGrowthStartMargin);

  double a = aStart, D = aStart, dD = 1.;

  growth_advance(a, D, dD, aFirst, kGrowthLeadSteps, OmegaM, w0, wa);
  m_growthRate[nSamples-1] = a * dD / D;

  for (unsigned i = nSamples - 1; i > 0; --i) {
    const double target = 1. / (1. + m_redshift[i-1]);
    growth_advance(a, D, dD, target, kGrowthSteps, OmegaM, w0, wa);
    m_growthRate[i-1] = a * dD / D;
  }

  for (unsigned i = 0; i < nSamples; ++i)
    if (!std::isfinite(m_growthRate[i]) || m_growthRate[i] <= 0.)
      throw Error("the growth rate came out non-positive or non-finite at z = " +
                  std::to_string(m_redshift[i]));
}


// ============================================================================


otswap::DistanceTable::DistanceTable (std::vector<double> redshift,
                                      std::vector<double> distance,
                                      std::vector<double> growthRate)
  : m_redshift(std::move(redshift)), m_distance(std::move(distance)),
    m_growthRate(std::move(growthRate))
{
  if (m_redshift.size() != m_distance.size())
    throw Error("the redshift and distance arrays differ in length: " +
                std::to_string(m_redshift.size()) + " against " +
                std::to_string(m_distance.size()));

  if (m_redshift.size() < 2)
    throw Error("the table holds " + std::to_string(m_redshift.size()) +
                " points; at least two are required");

  for (std::size_t i = 0; i < m_redshift.size(); ++i) {
    if (!std::isfinite(m_redshift[i]) || !std::isfinite(m_distance[i]))
      throw Error("the table holds a non-finite value at entry " + std::to_string(i));

    if (i > 0 && m_redshift[i] <= m_redshift[i-1])
      throw Error("the redshift array is not strictly increasing at entry " +
                  std::to_string(i));

    if (i > 0 && m_distance[i] <= m_distance[i-1])
      throw Error("the distance array is not strictly increasing at entry " +
                  std::to_string(i));
  }

  if (m_growthRate.empty()) return;

  if (m_growthRate.size() != m_redshift.size())
    throw Error("the growth-rate array holds " + std::to_string(m_growthRate.size()) +
                " entries against " + std::to_string(m_redshift.size()) +
                " redshifts; it must either match or be empty");

  // The growth rate need not be monotonic, so only finiteness and sign
  // are checked.
  for (std::size_t i = 0; i < m_growthRate.size(); ++i)
    if (!std::isfinite(m_growthRate[i]) || m_growthRate[i] <= 0.)
      throw Error("the growth-rate array holds a non-positive or non-finite value at entry " +
                  std::to_string(i));
}


// ============================================================================


double otswap::DistanceTable::distanceAt (const double z) const
{
  if (!std::isfinite(z))
    throw Error("the requested redshift is not finite");

  if (z < m_redshift.front() || z > m_redshift.back())
    throw Error("the redshift " + std::to_string(z) + " is outside the tabulated range [" +
                std::to_string(m_redshift.front()) + ", " +
                std::to_string(m_redshift.back()) + "]");

  if (z == m_redshift.front()) return m_distance.front();
  if (z == m_redshift.back()) return m_distance.back();

  const std::size_t i =
    (std::size_t)(std::lower_bound(m_redshift.begin(), m_redshift.end(), z) - m_redshift.begin());
  const double f = (z - m_redshift[i-1]) / (m_redshift[i] - m_redshift[i-1]);
  return m_distance[i-1] + f * (m_distance[i] - m_distance[i-1]);
}


// ============================================================================


double otswap::DistanceTable::redshiftAt (const double distance) const
{
  if (!std::isfinite(distance))
    throw Error("the requested comoving distance is not finite");

  if (distance < m_distance.front() || distance > m_distance.back())
    throw Error("the comoving distance " + std::to_string(distance) +
                " is outside the tabulated range [" + std::to_string(m_distance.front()) +
                ", " + std::to_string(m_distance.back()) + "]");

  if (distance == m_distance.front()) return m_redshift.front();
  if (distance == m_distance.back()) return m_redshift.back();

  const std::size_t i =
    (std::size_t)(std::lower_bound(m_distance.begin(), m_distance.end(), distance) - m_distance.begin());
  const double f = (distance - m_distance[i-1]) / (m_distance[i] - m_distance[i-1]);
  return m_redshift[i-1] + f * (m_redshift[i] - m_redshift[i-1]);
}


// ============================================================================


double otswap::DistanceTable::growthRateAt (const double z) const
{
  if (m_growthRate.empty())
    throw Error("this distance table carries no growth rate: it was built from arrays "
                "without one");

  if (!std::isfinite(z))
    throw Error("the requested redshift is not finite");

  if (z < m_redshift.front() || z > m_redshift.back())
    throw Error("the redshift " + std::to_string(z) + " is outside the tabulated range [" +
                std::to_string(m_redshift.front()) + ", " +
                std::to_string(m_redshift.back()) + "]");

  return internal::profile_at(m_redshift, m_growthRate, z);
}


// ============================================================================


bool otswap::DistanceTable::hasGrowthRate () const
{
  return !m_growthRate.empty();
}


// ============================================================================


double otswap::DistanceTable::minRedshift () const
{
  return m_redshift.front();
}


// ============================================================================


double otswap::DistanceTable::maxRedshift () const
{
  return m_redshift.back();
}
