# otswap

<!-- docs:about-start -->
otswap reconstructs the Lagrangian displacement field of a catalogue of
large-scale structure tracers, such as galaxies, quasars or dark matter
haloes, by discrete optimal transport. Each tracer is paired one to one
with a point drawn from a random catalogue that samples the same volume
uniformly, standing for the homogeneous initial conditions, so that the
sum of squared comoving distances between paired points is minimal. The
pairing starts from a local first guess, which matches each tracer to a
nearby free random, and is then refined by sweeps of local swaps: for
four neighbouring tracers, all 24 permutations of their partners are
tried and the cheapest is kept, until the fraction of successful swaps
in a sweep no longer exceeds a chosen threshold.

The reconstruction is repeated on independent sets of N points extracted
from the random catalogue. A tracer's displacement is the mean, over
these realizations, of the vector from the tracer to its partner. As the
number of realizations grows, the mean Lagrangian position tends to the
centroid of the equal-volume Laguerre cell that optimal transport
assigns to the tracer: the Lagrangian volume from which the observed
object has evolved.

otswap is built for survey data as much as for simulations. It handles
boxes and lightcones; in a lightcone the footprint enters through a
HEALPix mask and through a random catalogue reproducing the selection
function, with an optional redshift cut, and displacements crossing
unobserved regions are rejected. From the reconstructed displacements it
corrects the catalogue for redshift-space distortions, given a fiducial
cosmology and a bias b(z), and it writes its results as ASCII or FITS
tables. Results are identical bit for bit on every supported platform
and compiler, whatever the number of threads (see Reproducibility).
C++17, with a Python package on the same code.
<!-- docs:about-end -->

It is the middle of three packages: it depends on
[meshsearch](https://github.com/SimoSartori/meshsearch) for its neighbour
searches, and the Back-in-time Void Finder will be built on it. Extracted from
CosmoBolognaLib, it depends on nothing from it.

Documentation, for both APIs: **<https://simosartori.github.io/otswap/>**

## Install

```
pip install otswap
```

Wheels are built for CPython 3.9 to 3.14, on Linux (x86_64, aarch64) and macOS
(x86_64, arm64). From a clone, `pip install .` builds the package from source,
with the requirements below.

## Build from source

Needs CMake 3.16 or later, a C++17 compiler and cfitsio; OpenMP when found.
meshsearch 0.2 is used when installed, and fetched from GitHub otherwise.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix
```

Set `CFITSIO_ROOT` (and `ZLIB_ROOT` for a static cfitsio) where cfitsio is not
in a system path. The options:

- `OTSWAP_BUILD_TESTS`: the C++ test suite, on when otswap is the top-level
  project;
- `OTSWAP_BUILD_PYTHON`: the Python extension, off (`pip install .` sets it);
- `OTSWAP_INSTALL`: the install rules;
- `OTSWAP_MESHSEARCH_REPOSITORY`, `OTSWAP_MESHSEARCH_TAG`: where meshsearch is
  fetched from, when it is not installed; the repository may be a local path.

### Using it from CMake

```cmake
find_package(otswap 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE otswap::otswap)
```

Configure your project with `-DCMAKE_PREFIX_PATH=/path/to/prefix`. otswap can
also be built as part of your project, with `add_subdirectory` or
FetchContent, and the same target.

## Reproducibility

<!-- docs:reproducibility-start -->
With a fixed seed, a reconstruction gives the same result whatever the number
of threads, bit for bit. The seed 0 draws one, and the seed used is recorded in
the result, so any run can be repeated. Every translation unit is compiled with
`-ffp-contract=off`, and the math functions the results depend on come from a
vendored copy of FreeBSD's msun rather than the platform's library, so the same
inputs give the same bits on every supported platform and compiler.
<!-- docs:reproducibility-end -->

## Tests

```
ctest --test-dir build
OTSWAP_CPP_REFERENCE=build/cpp_reference pytest tests/
```

The Python tests run against the installed package; those comparing it with
the C++ library need the `cpp_reference` program of the C++ build, and are
skipped without it. `tests/test_determinism.py` checks the hashes of the
results recorded in `tests/determinism_hashes.json`, with one thread and with
four.

## Example

```python
import otswap

tracers = otswap.io.read("tracers.dat", [0, 1, 2]).values   # (N, 3), Mpc/h

result = otswap.reconstruct_box(tracers, n_realizations=8, seed=12345)
print(result.mean_displacement[:3])

otswap.io.write_displacements("displacements.dat", result)
```

`examples/` holds four examples, each in C++ and in Python, with the data to
run them: box and lightcone reconstructions, and the redshift-space correction
in each geometry. The documentation walks through them.

## Python

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift, with
the angular unit given explicitly by `angle_unit`. Inputs are converted to
float64; the reconstruction works on its own copy of them. The number of
threads follows `OMP_NUM_THREADS`, and the reconstruction functions release the
GIL while they run. Tables are read and written by `otswap.io`.

## Citing

Please cite Sartori et al. 2026,
[arXiv:2601.15378](https://arxiv.org/abs/2601.15378), and CosmoBolognaLib,
Marulli, Veropalumbo & Moresco 2016, Astronomy and Computing, 14, 35.
`CITATION.cff` has the details.

## Licence

GPL-2.0-or-later; see `LICENSE`. Vendored code keeps its own terms:

- `extern/fdlibm`: FreeBSD msun, under the Sun Microsystems notice in each file,
  a permissive licence (see its README);
- `extern/healpix`: HEALPix 3.82, GPL-2.0-or-later, built only with the tests.

The example data in `examples/data` are not covered by this licence: they may
be redistributed unchanged with the examples and used only to run them (see
`examples/data/README.md`).
