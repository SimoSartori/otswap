# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =====================================================================
# Example code: how to move a box catalogue from redshift to real space
# =====================================================================

"""This example explains how to correct a catalogue in box geometry for
redshift-space distortions with the OT reconstruction. The tracers are given
in redshift space, distorted along one Cartesian axis, the line of sight.
The example runs the reconstruction of the box example on them, then moves
each tracer along that axis to the position it would have in real space.

The correction of tracer i is a shift along the axis,

    s_i = f(z) / (b + 3 f(z) / 5) * <Psi_axis>_i ,

with f the linear growth rate at the box's redshift z, b the tracers' linear
bias, and <Psi_axis>_i the component along the axis of the reconstructed
mean displacement, averaged with a gaussian of width sigma over the tracers
around i.

The box is not periodic, here or in the reconstruction: nothing flows
through its faces, so modes on the scale of the box itself are not
reconstructed, and the corrected positions are not wrapped. A box cut from a
periodic simulation is treated as a plain sub-volume.

Usage: python full_recons_box.py [data_dir]. The catalogues are read from
data_dir, by default the examples' data folder; the output is written to the
folder output/ next to this file, created if missing.
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

    # Comoving Cartesian coordinates in Mpc/h, shape (N, 3). The tracers are
    # in redshift space, with the line of sight along y; the randoms are
    # uniform in the volume
    tracers = np.loadtxt(data / "box_halos_redshift_space.dat")
    randoms = np.loadtxt(data / "box_randoms.dat")

    print(f"Tracers: {len(tracers)}")
    print(f"Randoms: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # --------------------------------------------------------------
    # ------------ Compute the mean particle separation ------------
    # --------------------------------------------------------------

    # (V/N)^(1/3), with V the volume of the tracers' bounding box, as in the
    # box example
    extent = (tracers.max(axis=0) - tracers.min(axis=0)).tolist()
    mps = (extent[0] * extent[1] * extent[2] / len(tracers)) ** (1.0 / 3.0)

    print(f"Mean particle separation: {mps:g} Mpc/h")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # As in the box example: 8 realizations, a convergence threshold of 1e-2
    # and a fixed seed
    result = otswap.reconstruct_box(tracers, randoms, mps=mps, n_realizations=8, convergence=1e-2,
                                    seed=12345)

    print(f"Realizations: {result.n_realizations}")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # Only the linear growth rate f at the box's redshift is needed. The
    # table computes it from Omega_m, h, w0 and wa of a flat cosmology
    distances = otswap.DistanceTable.flat(0.3186, 0.67, w0=-1.0, wa=0.0)

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # The box's redshift, the tracers' linear bias, and the line of sight:
    # axis 0, 1 or 2 for x, y or z
    redshift, bias, axis = 1.0, 1.91255, 1

    # sigma, the width of the gaussian average in Mpc/h: 10 Mpc/h is a
    # reasonable starting value, but the best value depends on the sample
    # and should be checked in each analysis. With
    # weight_by_realizations=True, each neighbour would be weighted by its
    # number of valid realizations.
    #
    # Pass the positions the reconstruction was run on, in the same order. A
    # tracer with no valid neighbour within 3 sigma, which in a box happens
    # only for an isolated one, is left uncorrected, with NaN coordinates,
    # and listed in uncorrected
    catalogue = otswap.real_space_box(result, tracers, axis=axis, redshift=redshift,
                                      distances=distances, bias=bias, sigma=10.0,
                                      weight_by_realizations=False)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    print(f"RSD factor f/(b + 3f/5): {otswap.rsd_factor_box(redshift, distances, bias=bias):g}")

    positions = catalogue.positions
    corrected = ~np.isnan(positions[:, axis])
    shift = positions[corrected, axis] - tracers[corrected, axis]
    print(f"Tracers corrected: {np.count_nonzero(corrected)}")
    print(f"Tracers left uncorrected: {len(catalogue.uncorrected)}")
    print(f"Shift along the line of sight: mean {shift.mean():g}, rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: the corrected position,
    # NaN for an uncorrected tracer, then the number of valid realizations
    # and the two diagnostics of the average
    columns = [("X", False, "corrected position, in Mpc/h"),
               ("Y", False, "corrected position, in Mpc/h"),
               ("Z", False, "corrected position, in Mpc/h"),
               ("nValidRec", True, "number of valid OT realizations of the tracer"),
               ("nNeighbours", True, "number of tracers with a valid OT realization averaged within "
                                     "3 sigma, the tracer included if valid"),
               ("nRealizationsAveraged", True, "sum of nValidRec over those tracers")]
    rows = np.column_stack([positions, result.valid_realizations, catalogue.n_neighbours,
                            catalogue.n_realizations_averaged])
    write_table(output / "reconstructed_catalogue_box.dat", columns, rows)

    print("Written: output/reconstructed_catalogue_box.dat")


if __name__ == "__main__":
    main()
