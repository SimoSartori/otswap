# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =======================================================================
# Example code: how to run the OT reconstruction with lightcone geometry
# =======================================================================

"""This example explains how to reconstruct the displacement field of a
survey catalogue with lightcone geometry. The lightcone geometry takes
observed coordinates (right ascension, declination, redshift), a mean
particle separation that changes with redshift, which otswap measures from
the tracers, and the survey's footprint, here a HEALPix mask.

The reconstruction matches the tracers to randoms by optimal transport: in
each realization every tracer is paired with one random, and its
displacement points from the tracer (the Eulerian, observed position) to
that random (the Lagrangian position). Independent realizations are then
averaged.

The example first measures the mean particle separation alone, the quick
way to choose its redshift bins, and at the end reads the file it wrote
back.

The randoms carry the survey's geometry and selection function, which
otswap does not model: use the survey's own random catalogue.

Usage: python lightcone.py [data_dir]. The catalogues are read from
data_dir, by default the examples' data folder; the output is written to the
folder output/ next to this file, created if missing.
"""

import sys
from pathlib import Path

import numpy as np

import otswap


def main():

    # -------------------------------------------------------------
    # ------------ Locate the input and output folders ------------
    # -------------------------------------------------------------

    # The data folder is the examples' one, next to this example's folder,
    # unless a first command-line argument replaces it
    here = Path(__file__).resolve().parent
    data = Path(sys.argv[1]) if len(sys.argv) > 1 else here.parent / "data"
    output = here / "output"
    output.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------
    # ------------ Read the tracer and random catalogues ------------
    # ---------------------------------------------------------------

    # otswap.io.read reads the listed columns of a table: for ASCII files
    # 0-based indices, or the names of the "###" line of a file otswap
    # wrote; for FITS files the column names. Table.values has shape (N, 3):
    # right ascension and declination in degrees, and redshift. They are
    # passed as they are, with angle_unit="deg"
    tracers_sky = otswap.io.read(data / "lightcone_tracers.dat", [0, 1, 2]).values
    randoms_sky = otswap.io.read(data / "lightcone_randoms.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers_sky)}")
    print(f"Randoms read: {len(randoms_sky)}, {len(randoms_sky) / len(tracers_sky):g} per tracer")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    # A HEALPix map in RING or NESTED ordering, read from a FITS file: a
    # pixel is observed when its value is greater than 0. A map already in
    # memory can be passed to otswap.Mask.from_array instead. When a mask is
    # available, give it to the reconstruction rather than a sky area:
    # otswap then takes the area from it, and leaves out the objects on
    # unobserved pixels
    mask = otswap.Mask(str(data / "lightcone_mask.fits"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # The comoving distances, from Omega_m, h, w0 and wa of a flat
    # cosmology, tabulated by default from z = 0 to 10. The table must cover
    # the redshifts of tracers and randoms
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # ------------------------------------------------------------------
    # ------------ Set the parameters of the reconstruction ------------
    # ------------------------------------------------------------------

    # The first seven are the fields of an otswap.Config, the struct
    # otswap::Config of C++, passed to the call as config; a field left out
    # keeps its default. n_bins and redshift_cut are arguments of the call.
    #
    # n_realizations: independent realizations. Each one consumes as many
    # randoms as there are tracers, disjoint from those of the others, so
    # the randoms kept must number at least n_realizations times the
    # tracers kept. The redshift cut below keeps 202,184 of the 242,548
    # randoms and 24,708 of the 29,579 tracers (the mask keeps them all):
    # enough for 8 realizations.
    #
    # convergence: the sweeps stop once the fraction of successful swaps in a
    # sweep, swaps per tracer visited, no longer exceeds this threshold.
    # Larger values are faster and less accurate; 1e-2 keeps the example
    # fast. For an analysis, start from the defaults, 30 realizations and
    # 1e-3: the guide's section "Choosing the parameters" advises on both.
    #
    # seed: a fixed seed makes the run reproducible; 0 draws a new seed at
    # each run, and result.config.seed records the one drawn. The number
    # of threads follows OMP_NUM_THREADS, and with a fixed seed the result
    # is the same whatever the number of threads.
    #
    # cell_size: the grid cell, in units of the mean particle separation,
    # affects the speed only, never the result; the default is 3.
    #
    # reject_crossings, max_unobserved_pixels_crossed: after the
    # reconstruction, a displacement whose path on the sky (the great-circle
    # arc from the tracer to its random) crosses more than this many
    # unobserved pixels of the mask in a row is marked invalid. The count is
    # of pixels, so the same threshold means a different angle at a
    # different NSIDE. These are the defaults.
    #
    # verbosity: what the reconstruction prints. With "normal", the
    # default, how many tracers and randoms the redshift cut and the mask
    # left out, how many displacements the crossing filter rejected, and
    # then the time of the call and the tracers left without a valid
    # displacement. "detailed" adds the mps(z) profile, "silent" prints
    # nothing; the counts stay in result.selection either way.
    #
    # n_bins: the mean particle separation is measured in n_bins redshift
    # bins of equal width, over the redshifts of the tracers kept. Each bin
    # needs at least 278 tracers: if one is refused, the message gives an
    # upper bound on the number of bins the catalogue supports.
    #
    # redshift_cut: tracers and randoms outside this closed redshift range
    # are left out of the reconstruction; the catalogue spans 0.885 to 1.10,
    # and this range reconstructs a slice of it. The default, None, cuts
    # nothing; a (min, max) pair is accepted too
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, cell_size=3.0,
                           reject_crossings=True, max_unobserved_pixels_crossed=0, verbosity="normal")
    n_bins = 30
    redshift_cut = otswap.RedshiftCut(0.90, 1.08)

    # --------------------------------------------------------------
    # ------------ Measure the mean particle separation ------------
    # --------------------------------------------------------------

    # measure_mps_profile measures mps(z) from the tracers alone, as the
    # reconstruction does, from the same mask, bins, table and cut, and in a
    # fraction of its time: the profile is Result.mps_profile of the
    # reconstruction below, bit for bit. Use it to choose n_bins before a
    # reconstruction. It reports the tracers the cut and the mask left out
    profile = otswap.measure_mps_profile(tracers_sky, mask=mask, n_bins=n_bins, distances=distances,
                                         redshift_cut=redshift_cut, angle_unit="deg")

    print(f"mps(z): {len(profile.mps)} bins from z = {profile.redshift_min:g} to {profile.redshift_max:g}, "
          f"mps from {profile.mps.min():g} to {profile.mps.max():g} Mpc/h")

    # Without a mask, the keyword sky_area_deg2 gives the sky area, in square
    # degrees: no object is then left out on the sky and no crossing is
    # filtered, so the catalogues must lie within the footprint already.
    # These do, so with the mask's area the profile is the same
    from_area = otswap.measure_mps_profile(tracers_sky, sky_area_deg2=mask.sky_area_deg2, n_bins=n_bins,
                                           distances=distances, redshift_cut=redshift_cut, angle_unit="deg",
                                           verbosity="silent")

    print(f"The same profile from the sky area: {'yes' if np.array_equal(from_area.mps, profile.mps) else 'no'}")

    # One row per bin: its centre, its mps and its number of tracers, with
    # the range and the representative mps as keywords
    otswap.io.write_mps_profile(output / "mps_lightcone.dat", profile)

    print("Written: output/mps_lightcone.dat")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # Tracers left out by the cut or the mask keep their row in the result,
    # flagged in outside_redshift_cut or outside_mask, with NaN
    # displacements. The result also holds the tracers' Cartesian positions
    # (tracers), their mean Lagrangian positions (lagrangian, and
    # lagrangian_sky on the sky, in the angle unit of the call) and the
    # measured mps(z) (mps_profile)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=n_bins,
                                          distances=distances, angle_unit="deg",
                                          redshift_cut=redshift_cut, config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # valid holds, for each realization and tracer, whether the displacement
    # survived the crossing filter; valid_realizations counts them per
    # tracer, and mean_displacement averages the valid ones, NaN where there
    # are none
    n = len(tracers_sky)
    for rec in range(result.n_realizations):
        print(f"Realization {rec}: {np.count_nonzero(result.valid[rec])} valid displacements of {n}")

    counts = np.bincount(result.valid_realizations, minlength=result.n_realizations + 1)
    print("Tracers by number of valid realizations: " + ", ".join(f"{k}: {c}" for k, c in enumerate(counts)))
    left_out = np.count_nonzero(result.outside_redshift_cut | result.outside_mask)
    print(f"Tracers left out (redshift cut or mask): {left_out}")
    mean = result.mean_displacement
    no_mean = np.isnan(mean[:, 0])
    print(f"Tracers with a NaN mean displacement: {np.count_nonzero(no_mean)}")
    print(f"Mean length of the mean displacement: {np.linalg.norm(mean[~no_mean], axis=1).mean():g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: its sky coordinates
    # and those of its mean Lagrangian position, in degrees, its Cartesian
    # position, its mean Lagrangian position and its mean displacement, in
    # Mpc/h, its number of valid realizations and the two selection flags,
    # with a header describing the columns (angles in degrees, lengths in
    # coordinate units) and the parameters of the run, and, in two twin
    # files, the mps(z) profile and the distance table, so that
    # otswap.io.read_displacement_mean reads the file back whole. A tracer
    # left out, or without a valid realization, has NaN Lagrangian
    # coordinates and displacement. The writer is the C++ library's, so the
    # files are the C++ example's, byte for byte.
    # otswap.io.write_displacement_complete writes every realization instead,
    # losslessly (see the box example)
    otswap.io.write_displacement_mean(output / "displacement_lightcone.dat", result)

    print("Written: output/displacement_lightcone.dat, "
          "output/displacement_lightcone.mps.dat and "
          "output/displacement_lightcone.distances.dat")

    # --------------------------------------------
    # ------------ Read the file back ------------
    # --------------------------------------------

    # read_displacement_mean reads the file and its twins into a Result whose
    # only realization is the mean, with n_realizations_mean the number of
    # realizations averaged, the profile and the distance table: the values
    # as written, with 9 significant digits in ASCII, exactly in FITS. The
    # angles of the Result read are in angle_unit, which is required
    read_back = otswap.io.read_displacement_mean(output / "displacement_lightcone.dat", angle_unit="deg")

    largest = np.nanmax(np.abs(read_back.mean_displacement - result.mean_displacement))
    table = read_back.distances
    print(f"Read back: the mean of {read_back.n_realizations_mean} realizations, a profile of "
          f"{len(read_back.mps_profile.mps)} bins, a distance table from z = {table.min_redshift:g} to "
          f"{table.max_redshift:g}; largest difference from the mean written: {largest:g} Mpc/h")


if __name__ == "__main__":
    main()
