# otswap

```{include} ../README.md
:start-after: <!-- docs:about-start -->
:end-before: <!-- docs:about-end -->
```

The two APIs are one library: the C++ library and the Python package run the
same code, and write the same files from the same values; angles given in
degrees in Python are converted at the boundary, within rounding. [Names side
by side](naming.md) maps one to the other, and the [guide](guide.md) covers
what applies to both.

## Install

```{include} ../README.md
:start-after: <!-- docs:install-start -->
:end-before: <!-- docs:install-end -->
```

## In a few lines

A box reconstruction of tracers given as an array, its mean displacement, and
a file with one row per tracer:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```python
import otswap

config = otswap.Config(n_realizations=8, seed=12345)

result = otswap.reconstruct_box(tracers, config=config)
print(result.mean_displacement[0])   # tracer 0
otswap.io.write_displacement_mean("displacements.dat", result)
```
:::

:::{tab-item} C++
:sync: cpp

```cpp
#include <iostream>
#include "otswap/OT.h"

otswap::Config config;
config.nRealizations = 8;
config.seed = 12345;

otswap::Result result = otswap::reconstructBox(tracers, {}, config);
std::cout << result.meanDisplacement[0] << ' ' << result.meanDisplacement[1]
          << ' ' << result.meanDisplacement[2] << '\n';   // tracer 0
otswap::io::writeDisplacementMean("displacements.dat", result);
```
:::

::::

`tracers` is an array of shape (N, 3) in Python, a `std::vector<double>` of
`3 * N` entries in C++, in coordinate units. In a box, without randoms, otswap
draws them in the bounding box of the tracers. Eight realizations keep the
snippet fast; the default is 30. The [file formats](formats.md) page gives the
columns of the file.

This box is not periodic. For the whole volume of a simulation with periodic
boundary conditions, give the box limits and `periodic`, as in the [periodic
box example](examples.md#periodic-box); the [guide](guide.md) describes each
geometry step by step, and the [examples](examples.md) read real catalogues
and show every step.

## The pages

- [Guide](guide.md): the concepts and their use, in the order a new user
  needs them, with a path of steps for each geometry and starting values for
  the parameters.
- [Method](method.md): the problem otswap solves, how it solves it, and the
  redshift-space correction.
- [Examples](examples.md): six programs, each in C++ and in Python, on the
  data of the repository.
- [The executable](executable.md): a reconstruction and its correction from
  a parameter file, with no code to write.
- [File formats](formats.md): the tables otswap reads and writes.
- [C++ reference](cpp.rst) and [Python reference](python.rst): every member,
  its arguments and the errors it raises.
- [Names side by side](naming.md): the C++ and Python spelling of each member.
- [Using otswap from another package](packages.md): linking the library,
  sharing its Python types, and loading the same OpenMP runtime.

```{toctree}
:maxdepth: 2
:hidden:

guide
method
examples
executable
formats
cpp
python
naming
packages
changelog
citing
```
