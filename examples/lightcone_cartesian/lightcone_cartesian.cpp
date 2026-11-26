// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ============================================================================
// Example code: how to run the OT reconstruction on a lightcone in Cartesian
// coordinates, and correct it for redshift-space distortions
// ============================================================================

/*
  This example explains how to reconstruct a survey catalogue that is
  given in comoving Cartesian coordinates, with the observer at the
  origin, instead of right ascension, declination and redshift. Its
  tracers and randoms are those of the lightcone and
  full_recons_lightcone examples, converted to Cartesian coordinates
  with the same cosmology, and it does what both do: it measures the
  mean particle separation alone, to choose the bins, then runs the
  reconstruction, here with a cut on the distance, and the correction.
  It first finds the distance and the direction of each tracer, with
  pointToSky, which is where the cut and the mask act.

  Without redshifts the mean particle separation is measured as a
  function of the comoving distance from the origin, in uniform distance
  bins, and the mask and the optional distance cut act on the angles and
  the distance of each object. The distance table is optional for the
  reconstruction, and only fills the redshifts of the result; it is
  needed here for the correction, which takes the growth rate f and the
  bias b(z) at the tracers' redshifts. The example runs the
  reconstruction a second time without it, to show that nothing else
  changes.

  Usage: lightcone_cartesian [data_dir]. The catalogues are read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing.
*/

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "otswap/OT.h"


// ==========================================================================


int main (int argc, char** argv)
{
  try {

    // -------------------------------------------------------------
    // ------------ Locate the input and output folders ------------
    // -------------------------------------------------------------

    // OTSWAP_EXAMPLES_DATA and OTSWAP_EXAMPLE_OUTPUT are set by this
    // example's CMakeLists.txt to the data folder and to output/ next to
    // this file. A first command-line argument replaces the data folder
    const std::string data = (argc > 1) ? argv[1] : OTSWAP_EXAMPLES_DATA;
    const std::string output = OTSWAP_EXAMPLE_OUTPUT;
    std::filesystem::create_directories(output);


    // ---------------------------------------------------------------
    // ------------ Read the tracer and random catalogues ------------
    // ---------------------------------------------------------------

    // x, y, z for each object, in Mpc/h, the observer at the origin
    const std::vector<double> tracers = otswap::io::read(data + "/lightcone_tracers_cartesian.dat", {"0", "1", "2"}).values;
    const std::vector<double> randoms = otswap::io::read(data + "/lightcone_randoms_cartesian.dat", {"0", "1", "2"}).values;

    const std::size_t nTracers = tracers.size() / 3;
    const std::size_t nRandoms = randoms.size() / 3;

    std::cout << "Tracers read: " << nTracers << std::endl;
    std::cout << "Randoms read: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer" << std::endl;


    // ---------------------------------------
    // ------------ Read the mask ------------
    // ---------------------------------------

    const otswap::Mask mask(data + "/lightcone_mask.fits");

    std::cout << "Mask: NSIDE " << mask.nside() << ", " << mask.skyAreaDeg2() << " deg^2" << std::endl;


    // ---------------------------------------------------------------
    // ------------ Find the distances and the directions ------------
    // ---------------------------------------------------------------

    // pointToSky gives, for one point, its distance from the origin, in
    // coordinate units, and its right ascension and declination, in
    // radians: the direction the mask is applied to, and the distance the
    // distance cut and the mps bins are taken on. It needs no distance
    // table. The range of the distances is where to choose the cut.
    // otswap::skyToDegrees converts the angles of an array ordered right
    // ascension, declination, then a third value left unchanged, here the
    // distance: the degrees of Python's point_to_sky, bit for bit
    std::vector<double> directions(3*nTracers);
    std::vector<double> distance(nTracers);
    for (std::size_t i = 0; i < nTracers; ++i) {
      otswap::pointToSky(tracers[3*i], tracers[3*i+1], tracers[3*i+2],
                         distance[i], directions[3*i], directions[3*i+1]);
      directions[3*i+2] = distance[i];
    }
    directions = otswap::skyToDegrees(std::move(directions));

    std::vector<double> ra(nTracers), dec(nTracers);
    for (std::size_t i = 0; i < nTracers; ++i) {
      ra[i] = directions[3*i];
      dec[i] = directions[3*i+1];
    }

    std::cout << "Tracer distances: " << *std::min_element(distance.begin(), distance.end()) << " to "
              << *std::max_element(distance.begin(), distance.end()) << " Mpc/h" << std::endl;

    // One row per tracer, in input order: distance, RA, Dec
    otswap::io::write(output + "/tracers_sky_lightcone_cartesian.dat",
                      {{"distance", 'D', "distance from the origin", "", distance, {}},
                       {"ra", 'D', "right ascension", "deg", ra, {}},
                       {"dec", 'D', "declination", "deg", dec, {}}});

    std::cout << "Written: output/tracers_sky_lightcone_cartesian.dat" << std::endl;


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // The table gives the linear growth rate f(z) of the correction, and
    // the redshifts of the result. Its distances are in Mpc/h, so the
    // coordinates must be in Mpc/h too: the unit of the coordinates is
    // that of any distance table used with them. It must cover the
    // corrected distances too: the default range, z = 0 to 10, leaves room
    // for that
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // --------------------------------------------
    // ------------ Set the parameters ------------
    // --------------------------------------------

    // A fixed seed, and the displacements crossing an unobserved pixel of
    // the mask rejected. 8 realizations, which the example randoms are
    // enough for, and a convergence threshold of 1e-2 keep the example
    // fast; for an analysis, start from the defaults, 30 and 1e-3 (the
    // guide's section "Choosing the parameters"). The report goes to
    // std::cout, in order with this program's lines
    otswap::Config config;
    config.nRealizations = 8;
    config.convergence = 1.e-2;
    config.seed = 12345;
    config.rejectCrossings = true;
    config.maxUnobservedPixelsCrossed = 0;
    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;

    // The mean particle separation is measured in 30 distance bins of
    // equal width, each with enough tracers, over the distances of the
    // tracers kept
    const unsigned nBins = 30;

    // Tracers and randoms whose comoving distance from the origin lies
    // outside this closed range, in coordinate units, are left out of the
    // reconstruction: the randoms are dropped, and the tracers keep their
    // row, flagged in Result::outsideDistanceCut. The catalogue spans 2126
    // to 2500 Mpc/h. The default, {}, cuts nothing
    otswap::DistanceCut cut;
    cut.min = 2200.;
    cut.max = 2450.;


    // --------------------------------------------------------------
    // ------------ Measure the mean particle separation ------------
    // --------------------------------------------------------------

    // The Cartesian measureMpsProfile measures mps(d) from the tracers
    // alone, as the reconstruction does, from the same mask, bins and cut,
    // and in a fraction of its time: the profile is Result::mpsProfile of
    // the reconstruction below, bit for bit. Use it to choose nBins before
    // a reconstruction. It needs no distance table
    const otswap::MpsProfile profile =
      otswap::measureMpsProfile(tracers, mask, nBins, cut, otswap::Verbosity::Normal, &std::cout);

    std::cout << "mps(d): " << profile.mps.size() << " bins, mps from "
              << *std::min_element(profile.mps.begin(), profile.mps.end()) << " to "
              << *std::max_element(profile.mps.begin(), profile.mps.end()) << " Mpc/h" << std::endl;

    // Without a mask, the overloads that take the sky area, in square
    // degrees, leave out no object on the sky and filter no crossing: the
    // catalogues must then lie within the footprint already. These do, so
    // with the mask's area the profile is the same
    const otswap::MpsProfile fromArea =
      otswap::measureMpsProfile(tracers, mask.skyAreaDeg2(), nBins, cut, otswap::Verbosity::Silent);

    std::cout << "The same profile from the sky area: " << (fromArea.mps == profile.mps ? "yes" : "no") << std::endl;

    // One row per bin: its centre, in distance, its mps and its number of
    // tracers, with the range and the representative mps as keywords
    otswap::io::writeMpsProfile(output + "/mps_lightcone_cartesian.dat", profile);

    std::cout << "Written: output/mps_lightcone_cartesian.dat" << std::endl;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // The table, given by pointer since it is optional, fills the
    // redshifts of the result
    const otswap::Result result = otswap::reconstructLightcone(tracers, randoms, mask, nBins, config, &distances, cut);

    // The profile is in distance, in coordinate units (Mpc/h here): the
    // one measured above
    std::cout << "The profile of the result is the one measured: "
              << (result.mpsProfile.mps == profile.mps ? "yes" : "no") << std::endl;


    // ----------------------------------------------------------------
    // ------------ Reconstruct without the distance table ------------
    // ----------------------------------------------------------------

    // Without a table, a null pointer, the reconstruction does the same
    // work, from the same seed, and gives the same displacements, valid
    // flags and profile, bit for bit; only the redshifts of the result
    // are NaN, and its files hold no distance table
    const otswap::Result withoutTable = otswap::reconstructLightcone(tracers, randoms, mask, nBins, config,
                                                                     nullptr, cut);

    // Equal values, NaN counting as equal to NaN: the tracers left out by
    // the cut, and those without a valid realization, have NaN means
    auto sameValues = [] (const std::vector<double>& a, const std::vector<double>& b) {
      if (a.size() != b.size()) return false;
      for (std::size_t k = 0; k < a.size(); ++k)
        if (!(a[k] == b[k] || (std::isnan(a[k]) && std::isnan(b[k])))) return false;
      return true;
    };
    const bool same = sameValues(withoutTable.displacement, result.displacement) &&
                      withoutTable.valid == result.valid &&
                      withoutTable.mpsProfile.mps == result.mpsProfile.mps;
    std::cout << "Without the table: " << (same ? "the same" : "NOT the same")
              << " displacements, valid flags and profile" << std::endl;


    // ---------------------------------------------
    // ------------ Read the bias table ------------
    // ---------------------------------------------

    const otswap::BiasTable bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");


    // ----------------------------------------------------------------
    // ------------ Correct the redshift-space distortions ------------
    // ----------------------------------------------------------------

    // As in full_recons_lightcone: realSpaceLightcone takes the result of
    // either overload of reconstructLightcone. The distance table is
    // required here, even for a result reconstructed without one: the
    // growth rate f and the bias b(z) are taken at the tracers' redshifts,
    // which the table gives from their distances.
    //
    // sigma, the width of the Gaussian average in coordinate units (Mpc/h
    // here): 10 Mpc/h is a starting value, and the value must be
    // calibrated on the sample (the guide's section "Choosing the
    // parameters")
    const double sigma = 10.;
    otswap::CorrectionConfig correction;
    correction.weightByRealizations = false;
    correction.verbosity = otswap::Verbosity::Normal;
    correction.log = &std::cout;
    const otswap::RealSpaceCatalog catalogue = otswap::realSpaceLightcone(result, distances, bias, sigma, correction);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    double sum = 0., squares = 0., length = 0.;
    std::size_t corrected = 0, withMean = 0;
    for (std::size_t i = 0; i < nTracers; ++i) {
      const double* d = &result.meanDisplacement[3*i];
      if (!std::isnan(d[0])) {
        ++withMean;
        length += std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      }
      if (std::isnan(catalogue.shift[i])) continue;
      ++corrected;
      sum += catalogue.shift[i];
      squares += catalogue.shift[i] * catalogue.shift[i];
    }
    std::cout << "Mean length of the mean displacement over " << withMean << " tracers: "
              << length / (double)withMean << " Mpc/h" << std::endl;
    std::cout << "Shift along the line of sight over the " << corrected << " corrected tracers: mean "
              << sum / (double)corrected << ", rms " << std::sqrt(squares / (double)corrected) << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // The mean displacement, with the mps(d) profile and the distance
    // table in twin files, and the corrected catalogue; the columns are
    // those of the lightcone and full_recons_lightcone examples, but the
    // mean file adds the column outsideDistanceCut, and the profile is
    // binned in distance: a column distance and the keywords DISTMIN and
    // DISTMAX, in coordinate units, instead of the redshift ones. The
    // file of the result without the table has NaN redshift columns
    // (tracRed, lagrRed) and no distance table twin
    otswap::io::writeDisplacementMean(output + "/displacement_lightcone_cartesian.dat", result);
    otswap::io::writeDisplacementMean(output + "/displacement_lightcone_cartesian_no_table.dat", withoutTable);
    otswap::io::writeRealSpaceCatalog(output + "/reconstructed_catalogue_lightcone_cartesian.dat", catalogue);

    std::cout << "Written: output/displacement_lightcone_cartesian.dat" << std::endl;
    std::cout << "Written: output/displacement_lightcone_cartesian_no_table.dat" << std::endl;
    std::cout << "Written: output/reconstructed_catalogue_lightcone_cartesian.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
