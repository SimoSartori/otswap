# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.
"""

from typing import Literal

from ._otswap import (
    DistanceTable,
    Error,
    ExtrapolationWarning,
    Mask,
    RealSpaceCatalog,
    Result,
    SelectionCounts,
    line_of_sight_projection,
    neighbour_average,
    real_space_box,
    real_space_lightcone,
    reconstruct_box,
    reconstruct_lightcone,
    reject_mask_crossings,
    rsd_factor,
    rsd_factor_box,
    shift_along_line_of_sight,
    to_cartesian,
    to_sky,
)

AngleUnit = Literal["deg", "rad"]

__all__ = [
    "AngleUnit",
    "DistanceTable",
    "Error",
    "ExtrapolationWarning",
    "Mask",
    "RealSpaceCatalog",
    "Result",
    "SelectionCounts",
    "line_of_sight_projection",
    "neighbour_average",
    "real_space_box",
    "real_space_lightcone",
    "reconstruct_box",
    "reconstruct_lightcone",
    "reject_mask_crossings",
    "rsd_factor",
    "rsd_factor_box",
    "shift_along_line_of_sight",
    "to_cartesian",
    "to_sky",
]
