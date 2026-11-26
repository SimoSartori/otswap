# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ============================================================================
# Example code: how to run the OT reconstruction on a lightcone in Cartesian
# coordinates, and correct it for redshift-space distortions
# ============================================================================

"""This example explains how to reconstruct a survey catalogue that is given
in comoving Cartesian coordinates, with the observer at the origin, instead
of right ascension, declination and redshift. It is the lightcone and
full_recons_lightcone examples on the same tracers and randoms, converted to
Cartesian coordinates with the same cosmology.

Without redshifts the mean particle separation is measured as a function of
the comoving distance from the origin, in uniform distance bins, and the mask
and the optional distance cut act on the angles and the distance of each
object. The distance table is optional for the reconstruction, and only fills
the redshifts of the result; it is needed here for the growth rate of the
correction.

Usage: python lightcone_cartesian.py [data_dir]. The catalogues are read from
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

    # x, y, z for each object, in Mpc/h, the observer at the origin
    tracers = otswap.io.read(data / "lightcone_tracers_cartesian.dat", [0, 1, 2]).values
    randoms = otswap.io.read(data / "lightcone_randoms_cartesian.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers)}")
    print(f"Randoms read: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    mask = otswap.Mask(str(data / "lightcone_mask.fits"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # The table gives the linear growth rate f(z) of the correction, and the
    # redshifts of the result. It must cover the corrected distances too:
    # the default range, z = 0 to 10, leaves room for that
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # 8 realizations, a convergence threshold of 1e-2, a fixed seed, and the
    # displacements crossing an unobserved pixel of the mask rejected.
    #
    # The mean particle separation is measured in 30 distance bins of equal
    # width, each with enough tracers. A cut on the comoving distance, the
    # keyword distance_cut, a (d_min, d_max) pair, would leave out the
    # objects outside it; the default cuts nothing. angle_unit is required:
    # it is the unit of the angles in the result's sky coordinates
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, reject_crossings=True,
                           max_unobserved_pixels_crossed=0, verbosity="normal")
    result = otswap.reconstruct_lightcone_cartesian(tracers, randoms, mask=mask, n_bins=30,
                                                    distances=distances, angle_unit="deg", config=config)

    # The profile is in distance, in Mpc/h
    profile = result.mps_profile
    print(f"mps(d): {len(profile.mps)} nodes, from {profile.mps[0]:g} to {profile.mps[-1]:g} Mpc/h")

    # ---------------------------------------------
    # ------------ Read the bias table ------------
    # ---------------------------------------------

    bias = otswap.io.read_bias_table(data / "lightcone_bias.dat")

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # As in full_recons_lightcone: real_space_lightcone takes the result of
    # either function
    correction = otswap.CorrectionConfig(weight_by_realizations=False, verbosity="normal")
    catalogue = otswap.real_space_lightcone(result, distances=distances, bias=bias, sigma=10.0,
                                            correction=correction)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    mean = result.mean_displacement
    length = np.linalg.norm(mean, axis=1)
    length = length[~np.isnan(length)]
    print(f"Mean length of the mean displacement over {len(length)} tracers: {length.mean():g} Mpc/h")
    shift = catalogue.shift[~np.isnan(catalogue.shift)]
    print(f"Shift along the line of sight over the {len(shift)} corrected tracers: mean {shift.mean():g}, "
          f"rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # The mean displacement, with the mps(d) profile and the distance table
    # in twin files, and the corrected catalogue; the columns are those of
    # the lightcone and full_recons_lightcone examples. The writers are the
    # C++ library's, so the files are the C++ example's, byte for byte
    otswap.io.write_displacement_mean(output / "displacement_lightcone_cartesian.dat", result)
    otswap.io.write_real_space_catalog(output / "reconstructed_catalogue_lightcone_cartesian.dat", catalogue)

    print("Written: output/displacement_lightcone_cartesian.dat")
    print("Written: output/reconstructed_catalogue_lightcone_cartesian.dat")


if __name__ == "__main__":
    main()
