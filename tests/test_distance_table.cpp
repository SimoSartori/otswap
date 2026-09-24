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
 *  @file tests/test_distance_table.cpp
 *
 *  @brief The sampled redshift to comoving distance relation.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "otswap/OT.h"
#include "check.h"

using namespace otswap;

namespace {

  // Flat LCDM has a closed-form growth factor,
  //   D(a) proportional to E(a) * integral_0^a da' / (a' E(a'))^3,
  // which gives an independent reference for the integrated growth rate.
  double E_of_a (const double a, const double OmegaM)
  { return std::sqrt(OmegaM/(a*a*a) + (1.-OmegaM)); }

  double D_lcdm (const double a, const double OmegaM)
  {
    const int steps = 200000;
    const double h = a/steps;
    double sum = 0.;
    for (int k = 1; k <= steps; ++k) {
      const double x = (k-0.5)*h;
      const double t = x*E_of_a(x, OmegaM);
      sum += 1./(t*t*t);
    }
    return E_of_a(a, OmegaM)*sum*h;
  }

  double f_lcdm (const double a, const double OmegaM)
  {
    const double d = 1.e-5;
    return (std::log(D_lcdm(a*(1.+d), OmegaM)) - std::log(D_lcdm(a*(1.-d), OmegaM)))
         / (std::log(a*(1.+d)) - std::log(a*(1.-d)));
  }

}

int main ()
{
  group("the generated table is monotonic and covers its range");
  {
    const DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 2000);
    check_close(t.minRedshift(), 0., 1.e-15, "the table starts where it was asked to");
    check_close(t.maxRedshift(), 2., 1.e-15, "the table ends where it was asked to");
    check_close(t.distanceAt(0.), 0., 1.e-12, "the comoving distance vanishes at z = 0");

    double previous = -1.;
    for (int i = 0; i <= 200; ++i) {
      const double z = 2. * i / 200.;
      const double d = t.distanceAt(z);
      check(d > previous, "the comoving distance increases with redshift");
      previous = d;
    }
  }

  group("redshiftAt inverts distanceAt to the interpolation tolerance");
  {
    const DistanceTable t(0.31, 0.674, -1., 0., 0., 3., 4000);

    // Linear interpolation between samples of a smooth function has an
    // error of order (step^2/8) * |f''|; over 4000 samples on [0, 3] the
    // round trip is far tighter than this.
    for (int i = 1; i < 300; ++i) {
      const double z = 3. * i / 300.;
      check_close(t.redshiftAt(t.distanceAt(z)), z, 1.e-9, "z -> D_C -> z");
    }

    for (int i = 1; i < 300; ++i) {
      const double d = t.distanceAt(3.) * i / 300.;
      check_close(t.distanceAt(t.redshiftAt(d)), d, 1.e-6, "D_C -> z -> D_C");
    }
  }

  group("lookups outside the sampled range raise, and name the value");
  {
    const DistanceTable t(0.3, 0.7, -1., 0., 0.1, 1., 500);
    check_throws([&] { t.distanceAt(0.05); }, "below the first sampled redshift");
    check_throws([&] { t.distanceAt(1.5); }, "above the last sampled redshift");
    check_throws([&] { t.redshiftAt(-1.); }, "below the first sampled distance");
    check_throws([&] { t.redshiftAt(1.e9); }, "above the last sampled distance");
    check_throws([&] { t.distanceAt(std::nan("")); }, "a non-finite redshift");

    t.distanceAt(0.1);
    t.distanceAt(1.);
  }

  group("out-of-range parameters raise");
  {
    check_throws([] { DistanceTable(0., 0.7, -1., 0., 0., 1.); }, "OmegaM of zero");
    check_throws([] { DistanceTable(1.5, 0.7, -1., 0., 0., 1.); }, "OmegaM above one");
    check_throws([] { DistanceTable(0.3, 0., -1., 0., 0., 1.); }, "a non-positive h");
    check_throws([] { DistanceTable(0.3, 0.7, -1., 0., -0.1, 1.); }, "a negative zMin");
    check_throws([] { DistanceTable(0.3, 0.7, -1., 0., 1., 1.); }, "zMax equal to zMin");
    check_throws([] { DistanceTable(0.3, 0.7, -1., 0., 0., 1., 1); }, "a single sample");
  }

  group("a caller-supplied table is validated");
  {
    check_throws([] { DistanceTable({0., 1.}, {0., 1., 2.}); }, "arrays of different length");
    check_throws([] { DistanceTable({0.}, {0.}); }, "a table of one point");
    check_throws([] { DistanceTable({0., 1., 0.5}, {0., 1., 2.}); },
                 "a redshift array that is not increasing");
    check_throws([] { DistanceTable({0., 1., 2.}, {0., 2., 1.}); },
                 "a distance array that is not increasing");

    const DistanceTable t({0., 1., 2.}, {0., 100., 180.});
    check_close(t.distanceAt(0.5), 50., 1.e-12, "the supplied table interpolates linearly");
    check_close(t.redshiftAt(140.), 1.5, 1.e-12, "and inverts linearly");
  }

  group("the dark-energy equation of state changes the distance");
  {
    const DistanceTable lcdm(0.3, 0.7, -1., 0., 0., 1., 1000);
    const DistanceTable cpl(0.3, 0.7, -0.9, 0.2, 0., 1., 1000);
    check(lcdm.distanceAt(1.) != cpl.distanceAt(1.),
          "w0 and wa are not ignored");
  }

  group("the growth rate is exact in an Einstein-de Sitter universe");
  {
    // With OmegaM = 1, D is exactly proportional to a, so f is exactly 1
    // at every redshift. Nothing is fitted, so this has to come out.
    const DistanceTable t(1.0, 0.7, -1., 0., 0., 3., 600);
    check(t.hasGrowthRate(), "the generated table carries a growth rate");

    for (double z : {0., 0.5, 1., 2., 3.})
      check_close(t.growthRateAt(z), 1., 1.e-10, "f is one at every redshift");
  }

  group("the growth rate matches the closed-form LCDM solution");
  {
    const DistanceTable t(0.3, 0.7, -1., 0., 0., 3., 2000);

    for (double z : {0., 0.5, 1., 2., 3.}) {
      const double a = 1./(1.+z);
      const double reference = f_lcdm(a, 0.3);
      check_close(t.growthRateAt(z), reference, 1.e-5 * reference,
                  "f agrees with the closed form");
    }

    // It must also be a real solution of the equation, not a power of
    // Omega_m: f rises towards one at early times and falls below it as
    // dark energy takes over.
    check(t.growthRateAt(3.) > t.growthRateAt(0.), "f decreases towards the present");
    check(t.growthRateAt(3.) < 1., "and stays below one");
  }

  group("the growth rate responds to w0 and wa, which a fitted power would not");
  {
    const DistanceTable lcdm(0.3, 0.7, -1., 0., 0., 2., 1500);
    const DistanceTable cpl (0.3, 0.7, -0.8, 0.5, 0., 2., 1500);

    // Omega_m(z) is identical in the two, so any form that depended on it
    // alone would give the same f. The integrated solution does not.
    const double separation =
      std::fabs(lcdm.growthRateAt(0.5) - cpl.growthRateAt(0.5)) / lcdm.growthRateAt(0.5);
    check(separation > 0.05,
          "a different equation of state moves f by more than five per cent");
  }

  group("a table built without a growth rate says so, and refuses to supply one");
  {
    const DistanceTable t({0., 1., 2.}, {0., 100., 180.});
    check(!t.hasGrowthRate(), "hasGrowthRate is false");
    check_throws([&] { t.growthRateAt(1.); }, "growthRateAt raises");

    t.distanceAt(1.);
    t.redshiftAt(100.);
  }

  group("a supplied growth rate is interpolated and validated");
  {
    const DistanceTable t({0., 1., 2.}, {0., 100., 180.}, {0.5, 0.8, 0.95});
    check(t.hasGrowthRate(), "hasGrowthRate is true");
    check_close(t.growthRateAt(0.5), 0.65, 1.e-12, "f interpolates linearly between nodes");
    check_close(t.growthRateAt(1.), 0.8, 1.e-12, "and is exact on a node");
    check_throws([&] { t.growthRateAt(2.5); }, "outside the sampled range it raises");

    check_throws([] { DistanceTable({0., 1., 2.}, {0., 100., 180.}, {0.5, 0.8}); },
                 "a growth-rate array of the wrong length raises");
    check_throws([] { DistanceTable({0., 1., 2.}, {0., 100., 180.}, {0.5, 0., 0.95}); },
                 "a non-positive growth rate raises");
    check_throws([] { DistanceTable({0., 1., 2.}, {0., 100., 180.}, {0.5, std::nan(""), 0.95}); },
                 "a non-finite growth rate raises");

    // The growth rate need not be monotonic.
    const DistanceTable wiggly({0., 1., 2.}, {0., 100., 180.}, {0.8, 0.5, 0.9});
    check_close(wiggly.growthRateAt(1.), 0.5, 1.e-12, "a non-monotonic f is accepted");
  }

  group("toCartesian round-trips through redshiftAt, and its errors name the object");
  {
    const DistanceTable t(0.3, 0.7, -0.9, 0.1, 0., 2., 2000);
    const double pi = 3.14159265358979323846;

    // Objects all over the sky and the whole redshift range, including
    // both poles' neighbourhoods and both ends of the table.
    std::vector<double> sky;
    for (int i = 0; i < 400; ++i) {
      sky.push_back(2. * pi * (i % 37) / 37.);
      sky.push_back(-1.5 + 3. * (i % 23) / 22.);
      sky.push_back(2. * (i % 101) / 100.);
    }

    const std::vector<double> xyz = toCartesian(sky, t);
    check(xyz.size() == sky.size(), "three Cartesian coordinates per object");

    double zErr = 0., raErr = 0., decErr = 0.;
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      const double x = xyz[3*i], y = xyz[3*i+1], z = xyz[3*i+2];
      const double r = std::sqrt(x*x + y*y + z*z);
      zErr = std::max(zErr, std::fabs(t.redshiftAt(r) - sky[3*i+2]));
      if (r > 0.) {
        double ra = std::atan2(y, x);
        if (ra < 0.) ra += 2. * pi;
        raErr = std::max(raErr, std::fabs(std::remainder(ra - sky[3*i], 2. * pi)));
        decErr = std::max(decErr, std::fabs(std::asin(z / r) - sky[3*i+1]));
      }
    }
    check_close(zErr, 0., 1.e-12, "redshiftAt of the Cartesian radius recovers every redshift");
    check_close(raErr, 0., 1.e-12, "the right ascension is recovered");
    check_close(decErr, 0., 1.e-12, "and so is the declination");

    check(toCartesian({}, t).empty(), "an empty sky array gives an empty result");
    check_throws([&] { toCartesian({0.1, 0.2}, t); }, "a size that is not a multiple of three raises");
    check_throws([&] { toCartesian({0.1, 0.2, std::nan("")}, t); }, "a non-finite entry raises");

    bool named = false;
    try {
      toCartesian({0.1, 0.2, 0.5, 0.1, 0.2, 2.5}, t);
    }
    catch (const Error& e) {
      named = std::string(e.what()).find("object 1") != std::string::npos;
    }
    check(named, "a redshift outside the table raises, naming the object");
  }

  group("a value within a few rounding errors of either end is taken to be at it");
  {
    const double eps = std::numeric_limits<double>::epsilon();
    const DistanceTable t(0.3, 0.7, -1., 0., 0.1, 2., 500);
    const double dLo = t.distanceAt(0.1), dHi = t.distanceAt(2.);

    check(t.distanceAt(2. * (1. + 2. * eps)) == dHi, "just above the largest redshift");
    check(t.distanceAt(0.1 * (1. - 2. * eps)) == dLo, "just below the smallest redshift");
    check(t.redshiftAt(dHi * (1. + 4. * eps)) == 2., "just above the largest distance");
    check(t.redshiftAt(dLo * (1. - 4. * eps)) == 0.1, "just below the smallest distance");
    check(t.growthRateAt(2. * (1. + 2. * eps)) == t.growthRateAt(2.), "the growth rate too");

    check_throws([&] { t.distanceAt(2. + 1.e-12); }, "further above the range still raises");
    check_throws([&] { t.distanceAt(0.1 - 1.e-12); }, "and further below it");
    check_throws([&] { t.redshiftAt(dHi * (1. + 1.e-12)); }, "for distances as well");
    check_throws([&] { t.growthRateAt(2. + 1.e-12); }, "and for the growth rate");

    const DistanceTable supplied({0., 1., 2.}, {0., 100., 180.});
    check(supplied.redshiftAt(180. * (1. + 4. * eps)) == 2., "a supplied table has the same tolerance");
    check(supplied.distanceAt(-1.e-16) == 0., "at its lower end too");
    check_throws([&] { supplied.distanceAt(-1.e-12); }, "and raises beyond it");
  }

  group("the default table: relative error below 1e-6 for z >= 0.01, absolute below 2e-5 Mpc/h");
  {
    // Direct comoving distance, by composite Simpson on 1/E(z).
    auto direct = [] (const double z, const double OmegaM, const double w0, const double wa) {
      auto invE = [&] (const double x) {
        const double opz = 1. + x;
        const double de = std::exp(3. * (1. + w0 + wa) * std::log(opz) - 3. * wa * x / opz);
        return 1. / std::sqrt(OmegaM * opz * opz * opz + (1. - OmegaM) * de);
      };
      const int m = 2000;
      const double h = z / m;
      double sum = invE(0.) + invE(z);
      for (int k = 1; k < m; ++k) sum += (k % 2 ? 4. : 2.) * invE(k * h);
      return 2997.92458 * sum * h / 3.;
    };

    // The interpolation error peaks at the midpoints between samples: this
    // gives the midpoint of the default grid's interval holding z.
    const double step = 10. / (50000. - 1.);
    auto midpoint = [step] (const double z) { return (std::floor(z / step) + 0.5) * step; };

    for (const double wa : {0., 0.5}) {
      const double w0 = (wa == 0.) ? -1. : -0.8;
      const DistanceTable t(0.3, 0.7, w0, wa);
      check(t.minRedshift() == 0. && t.maxRedshift() == 10., "the default range is [0, 10]");

      for (const double z0 : {0.01, 0.37, 3.3, 9.99}) {
        const double z = midpoint(z0), d = direct(z, 0.3, w0, wa);
        check_close(t.distanceAt(z) / d, 1., 1.e-6, "for z >= 0.01, distanceAt is within 1e-6 relative");
        check_close(t.redshiftAt(d) / z, 1., 1.e-6, "and redshiftAt too");
      }

      // Over the whole range, including the first interval, where the
      // relative error is largest, and the points of largest absolute error.
      for (const double z0 : {0., 0.0093, 0.2313, 1., 9.99}) {
        const double z = midpoint(z0), d = direct(z, 0.3, w0, wa);
        check_close(t.distanceAt(z), d, 2.e-5, "distanceAt is within 2e-5 Mpc/h");
        check_close(direct(t.redshiftAt(d), 0.3, w0, wa), d, 2.e-5,
                    "and redshiftAt, as a distance along the line of sight, too");
      }
    }

    const DistanceTable lcdm(0.3, 0.7, -1., 0.);
    for (const double z0 : {0., 0.5, 2., 9.}) {
      const double z = midpoint(z0), reference = f_lcdm(1. / (1. + z), 0.3);
      check_close(lcdm.growthRateAt(z) / reference, 1., 1.e-6,
                  "growthRateAt is within 1e-6 relative over the whole range");
    }
  }

  return report("test_distance_table");
}
