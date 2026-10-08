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


def write_table(path, columns, rows):
    """Write rows in the ASCII layout of otswap's C++ io::write: a comment
    line per described column, the column names, then one line per row with
    9 significant digits, and integer columns as integers. columns is a list
    of (name, is_integer, description)."""
    described = [c for c in columns if c[2]]
    width = max((len(c[0]) for c in described), default=0)
    with open(path, "w") as f:
        for name, _, description in described:
            f.write(f"# {name}{' ' * (width - len(name) + 3)}{description}\n")
        if described:
            f.write("#\n")
        f.write("###" + "".join(f"   {c[0]}" for c in columns) + "\n")
        for row in rows:
            f.write(" ".join(str(int(v)) if c[1] else "%.9g" % v for c, v in zip(columns, row)) + "\n")


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
    tracers_sky = np.loadtxt(data / "lightcone_tracers.dat")
    randoms_sky = np.loadtxt(data / "lightcone_randoms.dat")

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
                                          max_unobserved_pixels_crossed=0, verbose=True)

    # ---------------------------------------------
    # ------------ Read the bias table ------------
    # ---------------------------------------------

    # The tracers' linear bias b(z), as nodes in redshift: b is interpolated
    # linearly between them, and extrapolated linearly beyond the first and
    # last, with an otswap.ExtrapolationWarning saying for how many tracers.
    # This table covers the whole catalogue
    bias_table = np.loadtxt(data / "lightcone_bias.dat")
    bias_redshift, bias = bias_table[:, 0], bias_table[:, 1]

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # sigma, the width of the gaussian average in Mpc/h: 10 Mpc/h is a
    # reasonable starting value, but the best value depends on the sample
    # and should be checked in each analysis. With
    # weight_by_realizations=True, each neighbour would be weighted by its
    # number of valid realizations.
    #
    # Pass the sky coordinates the reconstruction was run on, in the same
    # order. A tracer without a valid realization still receives the average
    # of its neighbours; one with no valid neighbour within 3 sigma, or left
    # out of the reconstruction by the mask, is left uncorrected, with NaN
    # coordinates, and listed in uncorrected
    catalogue = otswap.real_space_lightcone(result, tracers_sky, distances=distances,
                                            bias_redshift=bias_redshift, bias=bias, sigma=10.0,
                                            angle_unit="deg", weight_by_realizations=False)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # The factor f/(b + 3f/5) at the redshifts of the tracers that took part
    # in the reconstruction
    factor = otswap.rsd_factor(tracers_sky[~result.outside_mask, 2], distances,
                               bias_redshift=bias_redshift, bias=bias)
    print(f"RSD factor f/(b + 3f/5): from {factor.min():g} to {factor.max():g}")

    positions = catalogue.positions
    corrected = ~np.isnan(positions[:, 2])
    masked = int(np.count_nonzero(result.outside_mask))
    n_uncorrected = len(catalogue.uncorrected)
    print(f"Tracers corrected: {np.count_nonzero(corrected)}")
    print("  of which moved with the average of their neighbours (no valid realization): "
          f"{np.count_nonzero(corrected & (result.valid_realizations == 0))}")
    print(f"Tracers left uncorrected: {n_uncorrected} (left out by the mask: {masked}, "
          f"no valid tracer within 3 sigma: {n_uncorrected - masked})")
    shift = distances.distance_at(positions[corrected, 2]) - distances.distance_at(tracers_sky[corrected, 2])
    print(f"Shift along the line of sight, d(z') - d(z): mean {shift.mean():g}, "
          f"rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: the right ascension and
    # declination, unchanged, and the corrected redshift; NaN for an
    # uncorrected tracer. Then the number of valid realizations and the two
    # diagnostics of the average
    columns = [("RA", False, "corrected position, in degrees"),
               ("Dec", False, "corrected position, in degrees"),
               ("z", False, "corrected redshift"),
               ("nValidRec", True, "number of valid OT realizations of the tracer"),
               ("nNeighbours", True, "number of tracers with a valid OT realization averaged within "
                                     "3 sigma, the tracer included if valid"),
               ("nRealizationsAveraged", True, "sum of nValidRec over those tracers")]
    rows = np.column_stack([positions, result.valid_realizations, catalogue.n_neighbours,
                            catalogue.n_realizations_averaged])
    write_table(output / "reconstructed_catalogue_lightcone.dat", columns, rows)

    print("Written: output/reconstructed_catalogue_lightcone.dat")


if __name__ == "__main__":
    main()
