# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ===========================================================================
# Example code: how to move a lightcone catalogue from redshift to real space
# ===========================================================================

"""This example explains how to correct a survey catalogue for
redshift-space distortions with the OT reconstruction, in lightcone
geometry. It runs the reconstruction of the lightcone example, then moves
each tracer along its line of sight to the position it would have in real
space.

The correction of tracer i is a shift along its line of sight,

    s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5) * <Psi . r_hat>_i ,

with f the linear growth rate, b the tracers' linear bias, and
<Psi . r_hat>_i the line-of-sight component of the reconstructed mean
displacement, averaged with a gaussian of width sigma over the tracers
around i. The tracer keeps its right ascension and declination, and moves to
the redshift of comoving distance d(z_i) + s_i.

Usage: python full_recons_lightcone.py [data_dir]. The catalogues are read
from data_dir, by default the examples' data folder; the output is written to
the folder output/ next to this file, created if missing.
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

    # Right ascension and declination in degrees, and redshift; passed as
    # they are, with angle_unit="deg"
    tracers_sky = otswap.io.read(data / "lightcone_tracers.dat", [0, 1, 2]).values
    randoms_sky = otswap.io.read(data / "lightcone_randoms.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers_sky)}")
    print(f"Randoms read: {len(randoms_sky)}, {len(randoms_sky) / len(tracers_sky):g} per tracer")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    mask = otswap.Mask(str(data / "lightcone_mask.fits"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # The table gives the comoving distances and the linear growth rate f(z)
    # of the correction, computed from Omega_m, h, w0 and wa. It must cover
    # the corrected distances too, which can fall slightly outside the
    # catalogue's redshift range: the default range, z = 0 to 10, leaves
    # room for that
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # As in the lightcone example, without a redshift cut: 8 realizations,
    # a convergence threshold of 1e-2, a fixed seed, 30 redshift bins for
    # the mean particle separation, and the displacements crossing an
    # unobserved pixel of the mask rejected
    result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=30,
                                          distances=distances, angle_unit="deg", n_realizations=8,
                                          convergence=1e-2, seed=12345, reject_crossings=True,
                                          max_unobserved_pixels_crossed=0, verbosity="normal")

    # ---------------------------------------------
    # ------------ Read the bias table ------------
    # ---------------------------------------------

    # The tracers' linear bias b(z), as nodes in redshift: b is interpolated
    # linearly between them, and extrapolated linearly beyond the first and
    # last, with an otswap.ExtrapolationWarning saying for how many tracers.
    # This table covers the whole catalogue
    bias = otswap.io.read_bias_table(data / "lightcone_bias.dat")

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # sigma, the width of the gaussian average in Mpc/h: 10 Mpc/h is a
    # reasonable starting value, but the best value depends on the sample
    # and should be checked in each analysis. With
    # weight_by_realizations=True, each neighbour would be weighted by its
    # number of valid realizations. "detailed" prints, besides the line of
    # the call, how many tracers moved with their neighbours' average and
    # why the others were left uncorrected.
    #
    # The correction reads the tracers from the result. A tracer without a
    # valid realization still receives the average of its neighbours; one
    # with no valid neighbour within 3 sigma, or left out of the
    # reconstruction by the mask, is left uncorrected, with NaN coordinates;
    # catalogue.status says which
    catalogue = otswap.real_space_lightcone(result, distances=distances, bias=bias, sigma=10.0,
                                            weight_by_realizations=False, verbosity="detailed")

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # The factor f/(b + 3f/5) at the redshifts of the tracers that took part
    # in the reconstruction, and the shift applied along the line of sight
    # to each corrected tracer, both kept in the catalogue
    factor = catalogue.factor[~np.isnan(catalogue.factor)]
    print(f"RSD factor f/(b + 3f/5): from {factor.min():g} to {factor.max():g}")
    shift = catalogue.shift[~np.isnan(catalogue.shift)]
    print(f"Shift along the line of sight over the {len(shift)} corrected tracers: mean {shift.mean():g}, "
          f"rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: the right ascension and
    # declination, unchanged, and the corrected redshift; the corrected
    # Cartesian position; the number of valid realizations, the two
    # diagnostics of the average and the status of the tracer; NaN positions
    # for an uncorrected tracer. The writer is the C++ library's, so the file
    # is the C++ example's, byte for byte
    otswap.io.write_real_space_catalog(output / "reconstructed_catalogue_lightcone.dat", catalogue)

    print("Written: output/reconstructed_catalogue_lightcone.dat")


if __name__ == "__main__":
    main()
