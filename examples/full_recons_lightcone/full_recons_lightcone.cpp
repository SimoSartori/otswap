// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ===========================================================================
// Example code: how to move a lightcone catalogue from redshift to real space
// ===========================================================================

/*
  This example explains how to correct a survey catalogue for
  redshift-space distortions with the OT reconstruction, in lightcone
  geometry. It runs the reconstruction of the lightcone example, then
  moves each tracer along its line of sight to the position it would
  have in real space.

  The correction of tracer i is a shift along its line of sight,

    s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5) * <Psi . r_hat>_i ,

  with f the linear growth rate, b the tracers' linear bias, and
  <Psi . r_hat>_i the line-of-sight component of the reconstructed mean
  displacement, averaged with a gaussian of width sigma over the
  tracers around i. The tracer keeps its right ascension and
  declination, and moves to the redshift of comoving distance
  d(z_i) + s_i.

  Usage: full_recons_lightcone [data_dir]. The catalogues are read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing.
*/

#include <algorithm>
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

    // Right ascension, declination and redshift for each object, one row
    // after the other. The angles are read in degrees and converted to
    // radians by skyToRadians, with the factor numpy.deg2rad and otswap's
    // Python interface use
    const std::vector<double> tracersSky =
      otswap::skyToRadians(otswap::io::read(data + "/lightcone_tracers.dat", {"0", "1", "2"}).values);
    const std::vector<double> randomsSky =
      otswap::skyToRadians(otswap::io::read(data + "/lightcone_randoms.dat", {"0", "1", "2"}).values);

    const std::size_t nTracers = tracersSky.size() / 3;
    const std::size_t nRandoms = randomsSky.size() / 3;

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

    // The table gives the comoving distances and the linear growth rate
    // f(z) of the correction, computed from Omega_m, h, w0 and wa. It must
    // cover the corrected distances too, which can fall slightly outside
    // the catalogue's redshift range: the default range, z = 0 to 10,
    // leaves room for that
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // As in the lightcone example, without a redshift cut: 8
    // realizations, a convergence threshold of 1e-2, a fixed seed, 30
    // redshift bins for the mean particle separation, and the
    // displacements crossing an unobserved pixel of the mask rejected.
    // The report goes to std::cout, in order with this program's lines
    otswap::Config config;
    config.nRealizations = 8;
    config.convergence = 1.e-2;
    config.seed = 12345;
    config.rejectCrossings = true;
    config.maxUnobservedPixelsCrossed = 0;
    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;

    const otswap::Result result = otswap::reconstructLightcone(tracersSky, randomsSky, mask, 30, distances, config);


    // ---------------------------------------------
    // ------------ Read the bias table ------------
    // ---------------------------------------------

    // The tracers' linear bias b(z), as nodes in redshift: b is
    // interpolated linearly between them, and extrapolated linearly beyond
    // the first and last, which the correction reports and counts in
    // RealSpaceCatalog::nExtrapolated. This table covers the whole
    // catalogue
    const otswap::BiasTable bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");


    // ----------------------------------------------------------------
    // ------------ Correct the redshift-space distortions ------------
    // ----------------------------------------------------------------

    // sigma, the width of the gaussian average in Mpc/h: 10 Mpc/h is a
    // reasonable starting value, but the best value depends on the sample
    // and should be checked in each analysis. With
    // weightByRealizations, each neighbour would be weighted by its number
    // of valid realizations. Detailed reports, besides the line of the
    // call, how many tracers moved with their neighbours' average and why
    // the others were left uncorrected.
    //
    // The correction reads the tracers from the result. A tracer without a
    // valid realization still receives the average of its neighbours; one
    // with no valid neighbour within 3 sigma, or left out of the
    // reconstruction by the mask, is left uncorrected, with NaN
    // coordinates; RealSpaceCatalog::status says which
    const double sigma = 10.;
    otswap::CorrectionConfig correction;
    correction.weightByRealizations = false;
    correction.verbosity = otswap::Verbosity::Detailed;
    correction.log = &std::cout;
    const otswap::RealSpaceCatalog catalogue = otswap::realSpaceLightcone(result, distances, bias, sigma, correction);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // The factor f/(b + 3f/5) at the redshifts of the tracers that took
    // part in the reconstruction, and the shift applied along the line of
    // sight to each corrected tracer, both kept in the catalogue
    double factorMin = INFINITY, factorMax = -INFINITY, sum = 0., squares = 0.;
    std::size_t corrected = 0;
    for (std::size_t i = 0; i < nTracers; ++i) {
      if (!std::isnan(catalogue.factor[i])) {
        factorMin = std::min(factorMin, catalogue.factor[i]);
        factorMax = std::max(factorMax, catalogue.factor[i]);
      }
      if (std::isnan(catalogue.shift[i])) continue;
      ++corrected;
      sum += catalogue.shift[i];
      squares += catalogue.shift[i] * catalogue.shift[i];
    }
    std::cout << "RSD factor f/(b + 3f/5): from " << factorMin << " to " << factorMax << std::endl;
    std::cout << "Shift along the line of sight over the " << corrected << " corrected tracers: mean "
              << sum / (double)corrected << ", rms " << std::sqrt(squares / (double)corrected) << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One row per tracer, in the order of the input: the right ascension
    // and declination, unchanged, and the corrected redshift; the corrected
    // Cartesian position; the number of valid realizations, the two
    // diagnostics of the average and the status of the tracer; NaN
    // positions for an uncorrected tracer
    otswap::io::writeRealSpaceCatalog(output + "/reconstructed_catalogue_lightcone.dat", catalogue);

    std::cout << "Written: output/reconstructed_catalogue_lightcone.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
