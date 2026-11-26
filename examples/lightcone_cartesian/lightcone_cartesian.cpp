// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ============================================================================
// Example code: how to run the OT reconstruction on a lightcone in Cartesian
// coordinates, and correct it for redshift-space distortions
// ============================================================================

/*
  This example explains how to reconstruct a survey catalogue that is
  given in comoving Cartesian coordinates, with the observer at the
  origin, instead of right ascension, declination and redshift. It
  is the lightcone and full_recons_lightcone examples on the same
  tracers and randoms, converted to Cartesian coordinates with the same
  cosmology.

  Without redshifts the mean particle separation is measured as a
  function of the comoving distance from the origin, in uniform distance
  bins, and the mask and the optional distance cut act on the angles and
  the distance of each object. The distance table is optional for the
  reconstruction, and only fills the redshifts of the result; it is
  needed here for the growth rate of the correction.

  Usage: lightcone_cartesian [data_dir]. The catalogues are read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing.
*/

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
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


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // The table gives the linear growth rate f(z) of the correction, and
    // the redshifts of the result. It must cover the corrected distances
    // too: the default range, z = 0 to 10, leaves room for that
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // 8 realizations, a convergence threshold of 1e-2, a fixed seed, and
    // the displacements crossing an unobserved pixel of the mask rejected.
    // The report goes to std::cout, in order with this program's lines
    otswap::Config config;
    config.nRealizations = 8;
    config.convergence = 1.e-2;
    config.seed = 12345;
    config.rejectCrossings = true;
    config.maxUnobservedPixelsCrossed = 0;
    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;

    // The mean particle separation is measured in 30 distance bins of
    // equal width, each with enough tracers. A cut on the comoving
    // distance, otswap::DistanceCut {dMin, dMax}, would leave out the
    // objects outside it; the default, {}, cuts nothing
    const unsigned nBins = 30;

    const otswap::Result result = otswap::reconstructLightcone(tracers, randoms, mask, nBins, config, &distances);

    // The profile is in distance, in Mpc/h
    std::cout << "mps(d): " << result.mpsProfile.mps.size() << " nodes, from " << result.mpsProfile.mps.front()
              << " to " << result.mpsProfile.mps.back() << " Mpc/h" << std::endl;


    // ---------------------------------------------
    // ------------ Read the bias table ------------
    // ---------------------------------------------

    const otswap::BiasTable bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");


    // ----------------------------------------------------------------
    // ------------ Correct the redshift-space distortions ------------
    // ----------------------------------------------------------------

    // As in full_recons_lightcone: realSpaceLightcone takes the result of
    // either overload of reconstructLightcone
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
    // those of the lightcone and full_recons_lightcone examples
    otswap::io::writeDisplacementMean(output + "/displacement_lightcone_cartesian.dat", result);
    otswap::io::writeRealSpaceCatalog(output + "/reconstructed_catalogue_lightcone_cartesian.dat", catalogue);

    std::cout << "Written: output/displacement_lightcone_cartesian.dat" << std::endl;
    std::cout << "Written: output/reconstructed_catalogue_lightcone_cartesian.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
