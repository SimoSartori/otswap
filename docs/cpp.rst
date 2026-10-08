C++ reference
=============

Generated from the public headers ``otswap/OT.h``, ``otswap/RSD.h`` and
``otswap/io.h`` by Doxygen, with only what carries a comment included.

otswap/OT.h
-----------

Errors
~~~~~~

.. doxygenclass:: otswap::Error
   :members:

Configuration
~~~~~~~~~~~~~

.. doxygenstruct:: otswap::Config
   :members:

.. doxygenenum:: otswap::Verbosity

.. doxygenenum:: otswap::Geometry

.. doxygenstruct:: otswap::RedshiftCut
   :members:

Distances
~~~~~~~~~

.. doxygenclass:: otswap::DistanceTable
   :members:

.. doxygenfunction:: otswap::toCartesian

.. doxygenfunction:: otswap::toSky

.. doxygenfunction:: otswap::skyToRadians

.. doxygenfunction:: otswap::skyToDegrees

Mask
~~~~

.. doxygenclass:: otswap::Mask
   :members:

.. doxygenenum:: otswap::PixelOrdering

Reconstruction
~~~~~~~~~~~~~~

.. doxygenfunction:: otswap::reconstructBox(const std::vector<double>&, const std::vector<double>&, double, const Config&)

.. doxygenfunction:: otswap::reconstructBox(const std::vector<double>&, const std::vector<double>&, const Config&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, double, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, const Mask&, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, double, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const Mask&, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::rejectMaskCrossings

.. doxygenfunction:: otswap::recomputeMeans

Results
~~~~~~~

.. doxygenstruct:: otswap::Result
   :members:

.. doxygenstruct:: otswap::MpsProfile
   :members:

.. doxygenstruct:: otswap::SelectionCounts
   :members:

otswap/RSD.h
------------

The chain
~~~~~~~~~

.. doxygenfunction:: otswap::realSpaceLightcone

.. doxygenfunction:: otswap::realSpaceBox

.. doxygenstruct:: otswap::CorrectionConfig
   :members:

.. doxygenstruct:: otswap::RealSpaceCatalog
   :members:

.. doxygenenum:: otswap::CorrectionStatus

The four steps
~~~~~~~~~~~~~~

.. doxygenfunction:: otswap::radialProjection

.. doxygenfunction:: otswap::axisProjection

.. doxygenfunction:: otswap::neighbourAverage

.. doxygenstruct:: otswap::NeighbourAverage
   :members:

.. doxygenfunction:: otswap::rsdFactor

.. doxygenfunction:: otswap::rsdFactorBox

.. doxygenstruct:: otswap::BiasTable
   :members:

.. doxygenfunction:: otswap::shiftRadially

.. doxygenfunction:: otswap::shiftAlongAxis

otswap/io.h
-----------

Tables
~~~~~~

.. doxygenstruct:: otswap::io::Column
   :members:

.. doxygenstruct:: otswap::io::Table
   :members:

.. doxygenstruct:: otswap::io::Keyword
   :members:

.. doxygenstruct:: otswap::io::WriteOptions
   :members:

.. doxygenfunction:: otswap::io::read(const std::string&, const std::vector<std::string>&, char, char)

.. doxygenfunction:: otswap::io::read(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, char, char)

.. doxygenfunction:: otswap::io::write(const std::string&, const std::vector<Column>&, const WriteOptions&)

.. doxygenfunction:: otswap::io::write(const std::string&, const std::vector<Column>&, std::size_t, const std::function<void(std::size_t, std::vector<double>&)>&, const WriteOptions&)

.. doxygenfunction:: otswap::io::write(const std::string&, const std::vector<Column>&, std::size_t, const std::function<void(std::size_t, std::vector<double>&, std::vector<std::int64_t>&)>&, const WriteOptions&)

.. doxygenfunction:: otswap::io::readBiasTable

Writers
~~~~~~~

``io::writeDisplacements``, ``io::writeDisplacementField``,
``io::writeRealSpaceCatalog``, ``io::writeMpsProfile`` and their group enums
are documented on the :doc:`formats` page, with the layouts they write.
