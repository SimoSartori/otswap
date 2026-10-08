# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =================================================================
# Example code: how to run the OT reconstruction with box geometry
# =================================================================

"""This example explains how to reconstruct the displacement field of a
tracer catalogue with box geometry. The box geometry can be used with any
catalogue in comoving Cartesian coordinates, including non-cubic ones,
provided that the mean particle separation is approximately constant
throughout the catalogue.

The reconstruction matches the tracers to randoms by optimal transport: in
each realization every tracer is paired with one random, and its
displacement points from the tracer (the Eulerian, observed position) to
that random (the Lagrangian position). Independent realizations are then
averaged.

The box is not periodic: nothing flows through its faces, so modes on the
scale of the box itself are not reconstructed, and the displacements are
least reliable near the faces. A box cut from a periodic simulation is
treated as a plain sub-volume.

Usage: python box.py [data_dir]. The catalogues are read from data_dir, by
default the examples' data folder; the output is written to the folder
output/ next to this file, created if missing.
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
    # wrote; for FITS files the column names. Table.values has shape (N, 3),
    # x, y, z for each object; coordinates are comoving, in Mpc/h
    tracers = otswap.io.read(data / "box_halos_real_space.dat", [0, 1, 2]).values
    randoms = otswap.io.read(data / "box_randoms.dat", [0, 1, 2]).values

    print(f"Tracers: {len(tracers)}")
    print(f"Randoms: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # ------------------------------------------------------------------
    # ------------ Set the parameters of the reconstruction ------------
    # ------------------------------------------------------------------

    # n_realizations: independent realizations. Each one consumes N
    # randoms, disjoint from those of the others, so the random catalogue
    # must hold at least n_realizations x N objects: the catalogue read here
    # has exactly 8 per tracer. Omitting the randoms lets otswap draw them
    # itself, uniformly in the tracers' bounding box.
    #
    # convergence: the sweeps stop once the fraction of successful swaps in a
    # sweep, swaps per tracer visited, no longer exceeds this threshold.
    # Larger values are faster and less accurate; the default is 1e-3.
    #
    # seed: a fixed seed makes the run reproducible; 0 draws a new seed at
    # each run, and Result.seed records the one drawn. The number of threads
    # follows OMP_NUM_THREADS, and with a fixed seed the result is the same
    # whatever the number of threads.
    #
    # cell_size: the grid cell, in units of mps, affects the speed only,
    # never the result; the default is 4.
    #
    # verbosity: what the reconstruction prints; "detailed" adds, before the
    # line of the call, a line with the mean particle separation and its
    # source
    options = dict(n_realizations=8, convergence=1e-2, seed=12345, cell_size=4.0,
                   verbosity="detailed")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # In box geometry the mean particle separation (mps) is a single number.
    # Without one, as here, otswap takes (V/N)^(1/3), with V the volume of
    # the tracers' bounding box and N their number, and records it in
    # Result.mps; the keyword mps gives it instead. It sets the scale of the
    # search for swap partners, so it must describe the whole catalogue:
    # with a strongly varying density use the lightcone geometry, which
    # measures it as a function of redshift
    result = otswap.reconstruct_box(tracers, randoms, **options)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # In a box result every realization is valid for every tracer, and
    # mean_displacement is the average over all of them. Its average over
    # the tracers is the offset between the centres of mass of the randoms
    # used and of the tracers
    mean = result.mean_displacement
    all_valid = int(np.count_nonzero(result.valid_realizations == result.n_realizations))
    print(f"Realizations: {result.n_realizations}, all valid for {all_valid} tracers")
    print(f"Mean length of the mean displacement: {np.linalg.norm(mean, axis=1).mean():g} Mpc/h")
    print("Average of the mean displacement over the tracers: ({:g}, {:g}, {:g}) Mpc/h".format(*mean.mean(axis=0)))
    print("Displacement of the first tracer in the first realization: ({:g}, {:g}, {:g}) Mpc/h"
          .format(*result.displacement[0, 0]))

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: its position, its mean
    # Lagrangian position (the position plus the mean displacement), the
    # mean displacement and the number of valid realizations, with a header
    # giving the units and the parameters of the run. The writer is the C++
    # library's, so the file is the C++ example's, byte for byte;
    # otswap.io.write_displacement_field writes every realization instead,
    # losslessly
    otswap.io.write_displacements(output / "displacement_box.dat", result)

    print("Written: output/displacement_box.dat")


if __name__ == "__main__":
    main()
