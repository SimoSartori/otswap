# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ======================================================================
# Example code: how to run the OT reconstruction in a periodic box
# ======================================================================

"""This example explains how to reconstruct the displacement field of halos
in a cube with periodic boundary conditions, such as a simulation box. The
box is given by its limits; every distance is then the minimum image, so
nothing is cut at the faces and the displacements are as reliable at the
faces as inside, and modes on the scale of the box are reconstructed.

Each side of a periodic box must be an integer multiple of the shortest one,
so a cube, or a cuboid such as 500 x 500 x 1000, is accepted. Tracers
outside the limits are wrapped into the box, by a whole number of periods
per axis, and reported; the displacements point along the shortest way,
within half a period on each axis.

Usage: python periodic_box.py [data_dir]. The catalogue is read from
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

    # ----------------------------------------------------
    # ------------ Read the tracer catalogue ------------
    # ----------------------------------------------------

    # Halos of a simulation in a cube of side 1000 Mpc/h, [0, 1000) on each
    # axis, comoving coordinates in Mpc/h; Table.values has shape (N, 3)
    tracers = otswap.io.read(data / "periodic_halos.dat", [0, 1, 2]).values

    print(f"Tracers: {len(tracers)}")

    # ------------------------------------------------------------------
    # ------------ Set the parameters of the reconstruction ------------
    # ------------------------------------------------------------------

    # n_realizations: independent realizations, each on its own set of
    # randoms. With no randoms given, as here, otswap draws them itself,
    # uniformly in the box, a different set for each realization.
    #
    # convergence: the sweeps stop once the fraction of successful swaps in
    # a sweep no longer exceeds this threshold; the default is 1e-3.
    #
    # seed: a fixed seed makes the run reproducible, whatever the number of
    # threads.
    #
    # cell_size: the grid cell, in units of mps, affects the speed only,
    # never the result. In a periodic box the cells divide the shortest side
    # a whole number of times, so the cell used can be smaller:
    # result.cell_side
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, cell_size=3.0,
                           verbosity="normal")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # The box is periodic: the limits, the keyword box, are required, and
    # the mean particle separation, (V/N)^(1/3) with V the volume of the
    # limits, comes from them. A tracer outside the limits would be wrapped
    # into the box and counted in result.selection; here none is
    result = otswap.reconstruct_box(tracers, box=((0.0, 0.0, 0.0), (1000.0, 1000.0, 1000.0)),
                                    periodic=True, config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    print(f"Periodic: {'yes' if result.periodic else 'no'}, tracers wrapped: {result.selection.tracers_wrapped}")
    print(f"Mean particle separation: {result.mps:g} Mpc/h, grid cell side: {result.cell_side:g} Mpc/h")

    # Each displacement component lies within half a period, [-500, 500)
    # Mpc/h here, however far apart the tracer and its random lie in the
    # plain coordinates
    mean = result.mean_displacement
    print(f"Mean length of the mean displacement: {np.linalg.norm(mean, axis=1).mean():g} Mpc/h")
    print(f"Largest component of a mean displacement: {np.abs(mean).max():g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: its position, its mean
    # Lagrangian position, wrapped into the box, the mean displacement and
    # the number of valid realizations, with a header giving the units, the
    # box, the periodicity and the parameters of the run. The writer is the
    # C++ library's, so the file is the C++ example's, byte for byte
    otswap.io.write_displacement_mean(output / "displacement_periodic_box.dat", result)

    print("Written: output/displacement_periodic_box.dat")


if __name__ == "__main__":
    main()
