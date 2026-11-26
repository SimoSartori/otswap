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

The example then keeps every realization in a file, with the catalogue's
IDs, reads it back, and uses the realizations to check whether their number
is enough.

Without the periodic option, the box is not periodic: nothing flows
through its faces, so modes on the scale of the box itself are not
reconstructed, and the displacements are least reliable near the faces. A
box cut from a periodic simulation is treated as a plain sub-volume; a
whole simulation box with periodic boundary conditions is reconstructed
with the periodic option: see the periodic_box example.

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

    # They are the fields of an otswap.Config, the struct otswap::Config of
    # C++, passed to the call as config; a field left out keeps its default.
    #
    # n_realizations: independent realizations. Each one consumes N
    # randoms, disjoint from those of the others, so the random catalogue
    # must hold at least n_realizations x N objects: the catalogue read here
    # has exactly 8 per tracer (404,376 randoms for 50,547 halos), hence 8
    # realizations. Omitting the randoms lets otswap draw them itself,
    # uniformly in the box limits, or without them in the bounding box of
    # the tracers.
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
    # cell_size: the grid cell, in units of mps, affects the speed only,
    # never the result; the default is 3.
    #
    # verbosity: what the reconstruction prints: with box limits, how many
    # tracers and randoms lie outside them, then the line of the call;
    # "detailed" adds, before that line, the mean particle separation and
    # its source
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, cell_size=3.0,
                           verbosity="detailed")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # The catalogues fill a known volume, the cube [269, 669] Mpc/h, so it is
    # given as the box limits, the keyword box, closed on every face.
    # Tracers outside them would keep their row, flagged in
    # Result.outside_box with NaN displacements, and randoms outside them
    # would be dropped; here none is. Without limits the box is the bounding
    # box of the tracers and the randoms, as in the full_recons_box example.
    #
    # In box geometry the mean particle separation (mps) is a single number.
    # Without one, as here, otswap takes (V/N)^(1/3), with V the volume of
    # the limits and N the number of tracers inside them, and records it in
    # Result.mps; the keyword mps gives it instead. It sets the radius of the
    # first guess, 8 mps, within which each tracer is paired with a free
    # random, and the cell of the neighbour grids, so the density must be
    # roughly uniform across the box
    box = ([269.0, 269.0, 269.0], [669.0, 669.0, 669.0])

    result = otswap.reconstruct_box(tracers, randoms, box=box, config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # Every realization is valid for every tracer inside the limits, and
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
    # mean displacement, the number of valid realizations and the flag
    # outsideBox, with a header describing the columns (lengths in
    # coordinate units), the box and the parameters of the run. The writer
    # is the C++ library's, so the file is the C++ example's, byte for byte;
    # otswap.io.read_displacement_mean reads it back (see the lightcone
    # example)
    otswap.io.write_displacement_mean(output / "displacement_box.dat", result)

    print("Written: output/displacement_box.dat")

    # --------------------------------------------------------------
    # ------------ Keep every realization, with the IDs ------------
    # --------------------------------------------------------------

    # A catalogue usually names its objects by an ID column, read exactly by
    # otswap.io.read with integer_columns, and set as Result.ids. These
    # files have none, so here the line of each halo in its file, counted
    # from 1, stands for it. The reconstruction sets no IDs and no
    # computation reads them: the writers put them in a first column, id,
    # and real_space_box and real_space_lightcone copy them into the
    # corrected catalogue. None, or exactly 0..N-1, means none
    result.ids = np.arange(1, len(tracers) + 1)

    # write_displacement_complete writes every realization of every tracer,
    # with its valid flag, the whole configuration and the box, exactly:
    # with 17 significant digits in ASCII, in binary in FITS, as here, which
    # is smaller and faster to read. read_displacement_complete rebuilds the
    # Result, every double identical, the means recomputed from the
    # realizations; angle_unit is required, and applies to lightcone files
    otswap.io.write_displacement_complete(output / "displacement_complete_box.fits", result)

    print("Written: output/displacement_complete_box.fits")

    read_back = otswap.io.read_displacement_complete(output / "displacement_complete_box.fits", angle_unit="deg")

    identical = (np.array_equal(read_back.displacement, result.displacement)
                 and np.array_equal(read_back.valid, result.valid)
                 and np.array_equal(read_back.mean_displacement, result.mean_displacement)
                 and np.array_equal(read_back.ids, result.ids))
    ids = read_back.ids
    print(f"Read back: {read_back.n_realizations} realizations of {read_back.n_objects} tracers, "
          f"IDs {ids[0]} to {ids[-1]}, {'identical' if identical else 'NOT identical'} to the result written")

    # ----------------------------------------------------------
    # ------------ Check the number of realizations ------------
    # ----------------------------------------------------------

    # otswap.recompute_means recomputes valid_realizations, mean_displacement
    # and the Lagrangian positions from displacement and valid, for instance
    # after a filter of one's own has cleared some entries of valid. The
    # arrays of a Result are read-only in Python: Result.from_arrays builds
    # a new one from edited arrays, its means computed as recompute_means
    # computes them. Here it keeps the first half of the realizations. The
    # rms difference between the mean over half of them and the mean over
    # all of them shows how much the mean still depends on their number: the
    # larger it is against the mean length printed above, the more
    # realizations are worth adding
    n_half = read_back.n_realizations // 2
    valid = read_back.valid.copy()
    valid[n_half:] = False
    half = otswap.Result.from_arrays(read_back.displacement, valid, read_back.tracers)

    difference = half.mean_displacement - read_back.mean_displacement
    print(f"rms difference between the mean displacements of {n_half} and {read_back.n_realizations} "
          f"realizations: {np.sqrt((difference ** 2).sum() / len(tracers)):g} Mpc/h")


if __name__ == "__main__":
    main()
