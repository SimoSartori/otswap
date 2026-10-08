Python reference
================

Generated from the type stubs ``otswap/__init__.pyi`` and ``otswap/io.pyi``,
whose docstrings are the ones an IDE shows.

otswap
------

Errors
~~~~~~

.. autoapiexception:: otswap.Error

.. autoapiexception:: otswap.ExtrapolationWarning

Distances
~~~~~~~~~

.. autoapiclass:: otswap.DistanceTable
   :members:

.. autoapifunction:: otswap.to_cartesian

.. autoapifunction:: otswap.to_sky

Mask
~~~~

.. autoapiclass:: otswap.Mask
   :members:

Reconstruction
~~~~~~~~~~~~~~

.. autoapifunction:: otswap.reconstruct_box

.. autoapifunction:: otswap.reconstruct_lightcone

.. autoapifunction:: otswap.reject_mask_crossings

.. autoapifunction:: otswap.recompute_means

Results
~~~~~~~

.. autoapiclass:: otswap.Result
   :members:

.. autoapiclass:: otswap.MpsProfile
   :members:

.. autoapiclass:: otswap.SelectionCounts
   :members:

Redshift-space correction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoapifunction:: otswap.real_space_lightcone

.. autoapifunction:: otswap.real_space_box

.. autoapiclass:: otswap.RealSpaceCatalog
   :members:

.. autoapiclass:: otswap.CorrectionStatus
   :members:

.. autoapifunction:: otswap.radial_projection

.. autoapifunction:: otswap.axis_projection

.. autoapifunction:: otswap.neighbour_average

.. autoapiclass:: otswap.NeighbourAverage
   :members:

.. autoapifunction:: otswap.rsd_factor

.. autoapifunction:: otswap.rsd_factor_box

.. autoapiclass:: otswap.BiasTable
   :members:

.. autoapifunction:: otswap.shift_radially

.. autoapifunction:: otswap.shift_along_axis

Aliases
~~~~~~~

.. autoapidata:: otswap.AngleUnit

.. autoapidata:: otswap.Verbosity

otswap.io
---------

Tables
~~~~~~

.. autoapiclass:: otswap.io.Column
   :members:

.. autoapiclass:: otswap.io.Table
   :members:

.. autoapifunction:: otswap.io.read

.. autoapifunction:: otswap.io.write

.. autoapifunction:: otswap.io.read_bias_table

Writers
~~~~~~~

.. autoapifunction:: otswap.io.write_displacements

.. autoapidata:: otswap.io.DisplacementGroup

.. autoapifunction:: otswap.io.write_displacement_field

.. autoapifunction:: otswap.io.write_real_space_catalog

.. autoapidata:: otswap.io.CatalogGroup

.. autoapifunction:: otswap.io.write_mps_profile
