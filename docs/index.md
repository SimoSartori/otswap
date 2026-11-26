# otswap

```{include} ../README.md
:start-after: <!-- docs:about-start -->
:end-before: <!-- docs:about-end -->
```

The two APIs are one library: the C++ library and the Python package run the
same code and give the same results, bit for bit. [Names side by
side](naming.md) maps one to the other, and the [guide](guide.md) covers what
applies to both.

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
`3 * N` entries in C++, in coordinate units. Without randoms, otswap draws them in the
bounding box of the tracers. The [examples](examples.md) read real catalogues and
show every step.

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
changelog
citing
```
