// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =======================================================================
// Example code: how to run the OT reconstruction with lightcone geometry
// =======================================================================

/*
  This example explains how to reconstruct the displacement field of a
  survey catalogue with lightcone geometry. The lightcone geometry takes
  observed coordinates (right ascension, declination, redshift), a mean
  particle separation that changes with redshift, which otswap measures
  from the tracers, and the survey's footprint, here a HEALPix mask.

  The reconstruction matches the tracers to randoms by optimal
  transport: in each realization every tracer is paired with one
  random, and its displacement points from the tracer (the Eulerian,
  observed position) to that random (the Lagrangian position).
  Independent realizations are then averaged.

  The example first measures the mean particle separation alone, the
  quick way to choose its redshift bins, and at the end reads the file
  it wrote back.

  The randoms carry the survey's geometry and selection function, which
  otswap does not model: use the survey's own random catalogue.

  Usage: lightcone [data_dir]. The catalogues are read from data_dir, by
  default the examples' data folder; the output is written to the
  folder output/ next to this file, created if missing.
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

    // io::read reads the listed columns of a table: for ASCII files
    // 0-based indices, or the names of the "###" line of a file otswap
    // wrote; for FITS files (.fits, .fit, .fits.gz) the column names.
    // Table::values holds the rows one after the other: right ascension,
    // declination and redshift for each object, the flat layout every
    // otswap function takes. The catalogues give the angles in degrees,
    // and otswap takes radians: skyToRadians converts them with the factor
    // that otswap's Python interface uses
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

    // A HEALPix map in RING or NESTED ordering, read from a FITS file: a
    // pixel is observed when its value is greater than 0. A map already in
    // memory can be passed to the constructor taking an array instead.
    // When a mask is available, give it to the reconstruction rather than
    // a sky area: otswap then takes the area from it, and leaves out the
    // objects on unobserved pixels
    const otswap::Mask mask(data + "/lightcone_mask.fits");

    std::cout << "Mask: NSIDE " << mask.nside() << ", " << mask.skyAreaDeg2() << " deg^2" << std::endl;


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // The comoving distances, from Omega_m, h, w0 and wa of a flat
    // cosmology, tabulated by default from z = 0 to 10. The table must
    // cover the redshifts of tracers and randoms
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ------------------------------------------------------------------
    // ------------ Set the parameters of the reconstruction ------------
    // ------------------------------------------------------------------

    otswap::Config config;

    // Independent realizations. Each one consumes as many randoms as there
    // are tracers, disjoint from those of the others, so the randoms kept
    // must number at least nRealizations times the tracers kept. The
    // redshift cut below keeps 202,184 of the 242,548 randoms and 24,708 of
    // the 29,579 tracers (the mask keeps them all): enough for 8
    // realizations
    config.nRealizations = 8;

    // The sweeps stop once the fraction of successful swaps in a sweep,
    // swaps per tracer visited, no longer exceeds this threshold. Larger
    // values are faster and less accurate; 1e-2 keeps the example fast.
    // For an analysis, start from the defaults, 30 realizations and 1e-3:
    // the guide's section "Choosing the parameters" advises on both
    config.convergence = 1.e-2;

    // A fixed seed makes the run reproducible; 0 draws a new seed at each
    // run, and Result::config.seed records the one drawn. The number of
    // threads follows OMP_NUM_THREADS, and with a fixed seed the result is
    // the same whatever the number of threads
    config.seed = 12345;

    // The grid cell, in units of the mean particle separation, affects the
    // speed only, never the result; the default is 3
    config.cellSize = 3.;

    // After the reconstruction, a displacement whose path on the sky (the
    // great-circle arc from the tracer to its random) crosses more than
    // maxUnobservedPixelsCrossed unobserved pixels of the mask in a row is
    // marked invalid. The count is of pixels, so the same threshold means a
    // different angle at a different NSIDE. These are the defaults
    config.rejectCrossings = true;
    config.maxUnobservedPixelsCrossed = 0;

    // What the reconstruction reports: with Normal, the default, how many
    // tracers and randoms the redshift cut and the mask left out, how many
    // displacements the crossing filter rejected, and then the time of the
    // call and the tracers left without a valid displacement. Detailed adds
    // the mps(z) profile, Silent prints nothing; the counts stay in
    // Result::selection either way. The report goes to std::clog by
    // default; here to std::cout, in order with this program's own lines
    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;

    // The mean particle separation is measured in nBins redshift bins of
    // equal width, over the redshifts of the tracers kept. Each bin needs
    // at least 278 tracers: if one is refused, the message gives an upper
    // bound on the number of bins the catalogue supports
    const unsigned nBins = 30;

    // Tracers and randoms outside this closed redshift range are left out
    // of the reconstruction; the catalogue spans 0.885 to 1.10, and this
    // range reconstructs a slice of it. The default, {}, cuts nothing
    const otswap::RedshiftCut cut {0.90, 1.08};


    // --------------------------------------------------------------
    // ------------ Measure the mean particle separation ------------
    // --------------------------------------------------------------

    // measureMpsProfile measures mps(z) from the tracers alone, as the
    // reconstruction does, from the same mask, bins, table and cut, and in
    // a fraction of its time: the profile is Result::mpsProfile of the
    // reconstruction below, bit for bit. Use it to choose nBins before a
    // reconstruction. It reports the tracers the cut and the mask left out
    const otswap::MpsProfile profile =
      otswap::measureMpsProfile(tracersSky, mask, nBins, distances, cut, otswap::Verbosity::Normal, &std::cout);

    std::cout << "mps(z): " << profile.mps.size() << " bins from z = " << profile.redshiftMin << " to "
              << profile.redshiftMax << ", mps from " << *std::min_element(profile.mps.begin(), profile.mps.end())
              << " to " << *std::max_element(profile.mps.begin(), profile.mps.end()) << " Mpc/h" << std::endl;

    // Without a mask, the overloads that take the sky area, in square
    // degrees, leave out no object on the sky and filter no crossing: the
    // catalogues must then lie within the footprint already. These do, so
    // with the mask's area the profile is the same
    const otswap::MpsProfile fromArea =
      otswap::measureMpsProfile(tracersSky, mask.skyAreaDeg2(), nBins, distances, cut, otswap::Verbosity::Silent);

    std::cout << "The same profile from the sky area: " << (fromArea.mps == profile.mps ? "yes" : "no") << std::endl;

    // One row per bin: its centre, its mps and its number of tracers, with
    // the range and the representative mps as keywords
    otswap::io::writeMpsProfile(output + "/mps_lightcone.dat", profile);

    std::cout << "Written: output/mps_lightcone.dat" << std::endl;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // Tracers left out by the cut or the mask keep their row in the
    // result, flagged in outsideRedshiftCut or outsideMask, with NaN
    // displacements. The result also holds the tracers' Cartesian
    // positions (tracers), their mean Lagrangian positions (lagrangian, and
    // lagrangianSky on the sky) and the measured mps(z) (mpsProfile)
    const otswap::Result result = otswap::reconstructLightcone(tracersSky, randomsSky, mask, nBins,
                                                               distances, config, cut);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // valid holds, for each realization and tracer, whether the
    // displacement survived the crossing filter; validRealizations counts
    // them per tracer, and meanDisplacement averages the valid ones, NaN
    // where there are none
    for (unsigned rec = 0; rec < result.nRealizations; ++rec) {
      std::size_t n = 0;
      for (std::size_t i = 0; i < nTracers; ++i) n += result.valid[rec * nTracers + i];
      std::cout << "Realization " << rec << ": " << n << " valid displacements of " << nTracers << std::endl;
    }

    std::vector<std::size_t> byCount(result.nRealizations + 1, 0);
    std::size_t leftOut = 0, noMean = 0;
    double length = 0.;
    for (std::size_t i = 0; i < nTracers; ++i) {
      ++byCount[result.validRealizations[i]];
      leftOut += (result.outsideRedshiftCut[i] || result.outsideMask[i]) ? 1 : 0;
      const double* d = &result.meanDisplacement[3*i];
      if (std::isnan(d[0])) ++noMean;
      else length += std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    std::cout << "Tracers by number of valid realizations:";
    for (std::size_t n = 0; n < byCount.size(); ++n)
      std::cout << (n ? ", " : " ") << n << ": " << byCount[n];
    std::cout << std::endl;
    std::cout << "Tracers left out (redshift cut or mask): " << leftOut << std::endl;
    std::cout << "Tracers with a NaN mean displacement: " << noMean << std::endl;
    std::cout << "Mean length of the mean displacement: " << length / (double)(nTracers - noMean) << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One row per tracer, in the order of the input: its sky coordinates
    // and those of its mean Lagrangian position, in degrees, its Cartesian
    // position, its mean Lagrangian position and its mean displacement, in
    // Mpc/h, its number of valid realizations and the two selection flags,
    // with a header describing the columns (angles in degrees, lengths in
    // coordinate units) and the parameters of the run, and, in two twin
    // files, the mps(z) profile and the distance table, so that
    // io::readDisplacementMean reads the file back whole. A tracer left
    // out, or without a valid realization, has NaN Lagrangian coordinates
    // and displacement. io::writeDisplacementComplete writes every
    // realization instead, losslessly (see the box example)
    otswap::io::writeDisplacementMean(output + "/displacement_lightcone.dat", result);

    std::cout << "Written: output/displacement_lightcone.dat, "
                 "output/displacement_lightcone.mps.dat and "
                 "output/displacement_lightcone.distances.dat" << std::endl;


    // --------------------------------------------
    // ------------ Read the file back ------------
    // --------------------------------------------

    // readDisplacementMean reads the file and its twins into a Result whose
    // only realization is the mean, with nRealizationsMean the number of
    // realizations averaged, the profile and the distance table: the values
    // as written, with 9 significant digits in ASCII, exactly in FITS
    const otswap::Result readBack = otswap::io::readDisplacementMean(output + "/displacement_lightcone.dat");

    double largest = 0.;
    for (std::size_t k = 0; k < 3 * nTracers; ++k)
      if (!std::isnan(result.meanDisplacement[k]))
        largest = std::max(largest, std::fabs(readBack.meanDisplacement[k] - result.meanDisplacement[k]));
    std::cout << "Read back: the mean of " << readBack.nRealizationsMean << " realizations, a profile of "
              << readBack.mpsProfile.mps.size() << " bins, a distance table from z = "
              << readBack.distances->minRedshift() << " to " << readBack.distances->maxRedshift()
              << "; largest difference from the mean written: " << largest << " Mpc/h" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
