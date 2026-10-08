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

    # Arrays of shape (N, 3), x, y, z for each object; coordinates are
    # comoving, in Mpc/h
    tracers = np.loadtxt(data / "box_halos_real_space.dat")
    randoms = np.loadtxt(data / "box_randoms.dat")

    print(f"Tracers: {len(tracers)}")
    print(f"Randoms: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # --------------------------------------------------------------
    # ------------ Compute the mean particle separation ------------
    # --------------------------------------------------------------

    # In box geometry the mean particle separation (mps) is a single number
    # given by the caller: here (V/N)^(1/3), with V the volume of the
    # tracers' bounding box and N their number. It sets the scale of the
    # search for swap partners, so it must describe the whole catalogue:
    # with a strongly varying density use the lightcone geometry, which
    # measures it as a function of redshift
    extent = (tracers.max(axis=0) - tracers.min(axis=0)).tolist()
    mps = (extent[0] * extent[1] * extent[2] / len(tracers)) ** (1.0 / 3.0)

    print(f"Mean particle separation: {mps:g} Mpc/h")

    # ------------------------------------------------------------------
    # ------------ Set the parameters of the reconstruction ------------
    # ------------------------------------------------------------------

    # n_realizations: independent realizations. Each one consumes N
    # randoms, disjoint from those of the others, so the random catalogue
    # must hold at least n_realizations x N objects: the catalogue read here
    # has exactly 8 per tracer. Omitting the randoms lets otswap draw them
    # itself, uniformly in the tracers' bounding box.
    #
    # convergence: the swap loop stops when a sweep changes fewer than this
    # fraction of the pairs. Larger values are faster and less accurate; the
    # default is 1e-3.
    #
    # seed: a fixed seed makes the run reproducible; 0 draws a new seed at
    # each run. The number of threads follows OMP_NUM_THREADS, and with a
    # fixed seed the result is the same whatever the number of threads.
    #
    # cell_size: the grid cell, in units of mps, affects the speed only,
    # never the result; the default is 4
    options = dict(n_realizations=8, convergence=1e-2, seed=12345, cell_size=4.0)

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    result = otswap.reconstruct_box(tracers, randoms, mps=mps, **options)

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
    # Lagrangian position (the position plus the mean displacement), and the
    # mean displacement
    columns = [("tracX", False, "tracer position, in Mpc/h"), ("tracY", False, ""), ("tracZ", False, ""),
               ("lagrX", False, "mean Lagrangian position, tracer + mean displacement, in Mpc/h"),
               ("lagrY", False, ""), ("lagrZ", False, ""),
               ("displX", False, "mean displacement over the valid realizations, in Mpc/h"),
               ("displY", False, ""), ("displZ", False, "")]
    write_table(output / "displacement_box.dat", columns, np.hstack([tracers, tracers + mean, mean]))

    print("Written: output/displacement_box.dat")


if __name__ == "__main__":
    main()
