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

    # n_realizations: independent realizations. Each one consumes as many
    # randoms as there are tracers, disjoint from those of the others, so
    # the randoms kept must number at least n_realizations times the
    # tracers kept.
    #
    # convergence: the swap loop stops when a sweep changes fewer than this
    # fraction of the pairs. Larger values are faster and less accurate; the
    # default is 1e-3.
    #
    # seed: a fixed seed makes the run reproducible; 0 draws a new seed at
    # each run, and result.seed records the one drawn. The number of threads
    # follows OMP_NUM_THREADS, and with a fixed seed the result is the same
    # whatever the number of threads.
    #
    # cell_size: the grid cell, in units of the mean particle separation,
    # affects the speed only, never the result; the default is 4.
    #
    # reject_crossings, max_unobserved_pixels_crossed: after the
    # reconstruction, a displacement whose path on the sky (the great-circle
    # arc from the tracer to its random) crosses more than this many
    # unobserved pixels of the mask is marked invalid. The count is of
    # pixels, so the same threshold means a different angle at a different
    # NSIDE. These are the defaults.
    #
    # verbosity: what the reconstruction prints. With "normal", the
    # default, how many tracers and randoms the redshift cut and the mask
    # left out, how many displacements the crossing filter rejected, and
    # then the time of the call and the tracers left without a valid
    # displacement. "detailed" adds the mps(z) profile, "silent" prints
    # nothing; the counts stay in result.selection either way.
    #
    # n_bins: the mean particle separation is measured in n_bins redshift
    # bins of equal width. Each bin needs at least 278 tracers: if one is
    # refused, the message says how many bins the catalogue supports.
    #
    # redshift_cut: tracers and randoms outside this redshift range are left
    # out of the reconstruction. The range here spans the whole catalogue:
    # narrow it to reconstruct a slice. The default, None, cuts nothing
    options = dict(n_realizations=8, convergence=1e-2, seed=12345, cell_size=4.0,
                   reject_crossings=True, max_unobserved_pixels_crossed=0, verbosity="normal",
                   n_bins=30, redshift_cut=(0.885, 1.10))

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # Tracers left out by the cut or the mask keep their row in the result,
    # flagged in outside_redshift_cut or outside_mask, with NaN
    # displacements. The result also holds the tracers' Cartesian positions
    # (tracers), their mean Lagrangian positions (lagrangian, and
    # lagrangian_sky on the sky, in the angle unit of the call) and the
    # measured mps(z) (mps_profile)
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, distances=distances,
                                          angle_unit="deg", **options)

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
    # with a header giving the units and the parameters of the run. A tracer
    # left out, or without a valid realization, has NaN Lagrangian
    # coordinates and displacement. The writer is the C++ library's, so the
    # file is the C++ example's, byte for byte.
    # otswap.io.write_displacement_field writes every realization instead,
    # losslessly, and otswap.io.write_mps_profile the mps(z)
    otswap.io.write_displacements(output / "displacement_lightcone.dat", result)

    print("Written: output/displacement_lightcone.dat")


if __name__ == "__main__":
    main()
