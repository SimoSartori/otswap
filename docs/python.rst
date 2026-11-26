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

Configuration
~~~~~~~~~~~~~

.. autoapiclass:: otswap.Config
   :members:

.. autoapiclass:: otswap.CorrectionConfig
   :members:

.. autoapiclass:: otswap.RedshiftCut
   :members:

.. autoapiclass:: otswap.DistanceCut
   :members:

Distances
~~~~~~~~~

.. autoapiclass:: otswap.DistanceTable
   :members:

.. autoapifunction:: otswap.to_cartesian

.. autoapifunction:: otswap.to_sky

.. autoapifunction:: otswap.point_to_sky

Mask
~~~~

.. autoapiclass:: otswap.Mask
   :members:

Reconstruction
~~~~~~~~~~~~~~

.. autoapiclass:: otswap.BoxLimits
   :members:

.. autoapifunction:: otswap.reconstruct_box

.. autoapifunction:: otswap.reconstruct_lightcone

.. autoapifunction:: otswap.reconstruct_lightcone_cartesian

.. autoapifunction:: otswap.measure_mps_profile

.. autoapifunction:: otswap.measure_mps_profile_cartesian

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

OpenMP runtime
~~~~~~~~~~~~~~

.. autoapifunction:: otswap.openmp_root

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

.. autoapifunction:: otswap.io.read_keywords

.. autoapifunction:: otswap.io.column_names

.. autoapiclass:: otswap.io.TablePart
   :members:

.. autoapifunction:: otswap.io.write_parts

.. autoapifunction:: otswap.io.read_part

.. autoapifunction:: otswap.io.has_part

.. autoapifunction:: otswap.io.twin_file

.. autoapifunction:: otswap.io.read_bias_table

Writers
~~~~~~~

.. autoapifunction:: otswap.io.write_displacement_mean

.. autoapidata:: otswap.io.DisplacementGroup

.. autoapifunction:: otswap.io.read_displacement_mean

.. autoapifunction:: otswap.io.write_displacement_complete

.. autoapifunction:: otswap.io.read_displacement_complete

.. autoapifunction:: otswap.io.write_real_space_catalog

.. autoapidata:: otswap.io.CatalogGroup

.. autoapifunction:: otswap.io.write_mps_profile

.. _same-openmp-runtime:

Linking the same OpenMP runtime
-------------------------------

A process on macOS cannot hold two copies of the LLVM OpenMP runtime: the
second to start aborts it with ``OMP: Error #15``. A Python package that
compiles OpenMP code of its own, imported in the same process as otswap,
must therefore load the runtime otswap loads.

The macOS wheels of otswap ship it in the ``openmp`` directory of the
package, the directory ``otswap.openmp_root()`` returns:

- ``lib/libomp.dylib``, with install name ``@rpath/libomp.dylib``; the
  module ``_otswap`` finds it through its rpath ``@loader_path/openmp/lib``;
- ``include/omp.h`` and ``include/ompx.h``;
- ``LICENSE.TXT``, LLVM's licence: the macOS wheels include LLVM's libomp.

The runtime is the LLVM libomp otswap was built with; its version may change
between minor releases of otswap. On Linux, and in a build from source,
``otswap.openmp_root()`` returns None and otswap uses the runtime it was
built against.

A package built with CMake links the same file as follows.

1. It requires otswap at build time as well as at run time, in the same
   minor version, so that the directory is there when CMake runs:
   ``otswap >= 0.3, < 0.4`` in ``build-system.requires`` and in
   ``dependencies``.
2. Before ``find_package(OpenMP)``, it takes ``OpenMP_ROOT`` from otswap
   whenever otswap gives one, in a wheel build and in a build from source
   alike:

   .. code-block:: cmake

      execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import otswap; print(otswap.openmp_root() or '')"
        OUTPUT_VARIABLE otswap_openmp_root OUTPUT_STRIP_TRAILING_WHITESPACE)
      if (otswap_openmp_root)
        set(OpenMP_ROOT "${otswap_openmp_root}")
      endif()
      find_package(OpenMP REQUIRED COMPONENTS CXX)

   Its extension then records ``@rpath/libomp.dylib``.
3. It gives its extension the rpath of otswap's directory seen from its own,
   for two packages installed side by side in one ``site-packages``:

   .. code-block:: cmake

      set_target_properties(my_extension PROPERTIES
        INSTALL_RPATH "@loader_path/../otswap/openmp/lib")

   With otswap installed elsewhere (another ``--target``, an editable
   install), that rpath finds nothing.
4. Its wheel repair leaves the file where it is: with delocate,
   ``delocate-wheel --exclude libomp``. A copy in the package's own
   ``.dylibs/`` is a second runtime.
5. It imports otswap before its own extension.

This covers packages that follow it. Another package that ships its own
libomp, and is imported in the same process, can still clash with otswap and
with them.

``KMP_DUPLICATE_LIB_OK=TRUE`` silences the error but keeps two runtimes in
one process, which can crash or give wrong results: it is not a fix.
