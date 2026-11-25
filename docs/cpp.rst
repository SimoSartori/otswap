C++ reference
=============

Generated from the public headers ``otswap/OT.h``, ``otswap/RSD.h``,
``otswap/io.h``, ``otswap/ParameterFile.h``, ``otswap/Maths.h`` and
``otswap/Random.h`` by Doxygen,
with only what carries a comment included.

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

.. doxygenstruct:: otswap::DistanceCut
   :members:

.. doxygenstruct:: otswap::BoxLimits
   :members:

Distances
~~~~~~~~~

.. doxygenclass:: otswap::DistanceTable
   :members:

.. doxygenfunction:: otswap::toCartesian

.. doxygenfunction:: otswap::toSky

.. doxygenfunction:: otswap::pointToSky

.. doxygenfunction:: otswap::skyToRadians

.. doxygenfunction:: otswap::skyToDegrees

Mask
~~~~

.. doxygenclass:: otswap::Mask
   :members:

.. doxygenenum:: otswap::PixelOrdering

Reconstruction
~~~~~~~~~~~~~~

.. doxygenfunction:: otswap::reconstructBox(const std::vector<double>&, const std::vector<double>&, double, const Config&, const std::optional<BoxLimits>&, bool)

.. doxygenfunction:: otswap::reconstructBox(const std::vector<double>&, const std::vector<double>&, const Config&, const std::optional<BoxLimits>&, bool)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, double, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, const Mask&, unsigned, const DistanceTable&, const Config&, const RedshiftCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, double, unsigned, const Config&, const DistanceTable*, const DistanceCut&)

.. doxygenfunction:: otswap::reconstructLightcone(const std::vector<double>&, const std::vector<double>&, const Mask&, unsigned, const Config&, const DistanceTable*, const DistanceCut&)

.. doxygenfunction:: otswap::measureMpsProfile(const std::vector<double>&, double, unsigned, const DistanceTable&, const RedshiftCut&, Verbosity, std::ostream*)

.. doxygenfunction:: otswap::measureMpsProfile(const std::vector<double>&, const Mask&, unsigned, const DistanceTable&, const RedshiftCut&, Verbosity, std::ostream*)

.. doxygenfunction:: otswap::measureMpsProfile(const std::vector<double>&, double, unsigned, const DistanceCut&, Verbosity, std::ostream*)

.. doxygenfunction:: otswap::measureMpsProfile(const std::vector<double>&, const Mask&, unsigned, const DistanceCut&, Verbosity, std::ostream*)

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

.. doxygenfunction:: otswap::neighbourAverage(const std::vector<double>&, const std::vector<double>&, const std::vector<unsigned>&, double, bool)

.. doxygenfunction:: otswap::neighbourAverage(const std::vector<double>&, const std::vector<double>&, const std::vector<unsigned>&, double, bool, const BoxLimits&, bool)

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

.. doxygenfunction:: otswap::io::readKeywords

.. doxygenfunction:: otswap::io::columnNames

.. doxygenstruct:: otswap::io::TablePart
   :members:

.. doxygenfunction:: otswap::io::writeParts

.. doxygenfunction:: otswap::io::readPart

.. doxygenfunction:: otswap::io::hasPart

.. doxygenfunction:: otswap::io::twinFile

.. doxygenfunction:: otswap::io::readBiasTable

Writers
~~~~~~~

``io::writeDisplacementMean``, ``io::writeDisplacementComplete``,
``io::writeRealSpaceCatalog``, ``io::writeMpsProfile`` and their group enums
are documented on the :doc:`formats` page, with the layouts they write, and
so are ``io::readDisplacementMean`` and ``io::readDisplacementComplete``,
which read a mean and a complete displacement file back.

otswap/ParameterFile.h
----------------------

Not included by ``otswap/OT.h``; the syntax of the file is on the
:doc:`formats` page.

.. doxygenclass:: otswap::io::ParameterFile
   :members:

otswap/Maths.h
--------------

.. doxygenfile:: Maths.h

otswap/Random.h
---------------

Not included by ``otswap/OT.h``. Random draws over ``std::mt19937`` that give
the same values on every supported platform and compiler, unlike the standard
distributions and ``std::shuffle``. Their algorithms are fixed: no later
version changes them, so a program that draws with them keeps its values when
otswap is updated.

.. doxygenfunction:: otswap::uniformInt

.. doxygenfunction:: otswap::uniformReal

.. doxygenfunction:: otswap::shuffle
