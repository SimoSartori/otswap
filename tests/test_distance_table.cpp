/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file tests/test_distance_table.cpp
 *
 *  @brief The sampled redshift to comoving distance relation.
 */

#include <cmath>
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

  return report("test_distance_table");
}
