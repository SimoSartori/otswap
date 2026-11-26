# Using otswap from another package

This page is for the authors of a package built on otswap, in three
sections: a C++ program that links the library; a Python extension that
takes and returns otswap's types; and, on macOS, a Python extension with
OpenMP code of its own, imported beside otswap. A user of the Python package
alone needs none of it.

## From C++ with CMake

An installed otswap is found with `find_package`, and one built as part of
the project, with `add_subdirectory` or FetchContent, gives the same target:

```cmake
find_package(otswap 0.4 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE otswap::otswap)
```

`otswap::otswap` is a static library. Its CMake package finds what it needs
in turn: meshsearch 0.3 (0.3.0 or a later 0.3.x), cfitsio (honouring
`CFITSIO_ROOT`, and `ZLIB_ROOT` for a static cfitsio) and, when otswap was
built with it, OpenMP. **The program must use the OpenMP runtime family otswap
was built with**: libgomp for a library built with GCC, LLVM's libomp for
one built with Clang. With Apple's Clang, unless `OpenMP_ROOT` is set, as a
CMake or an environment variable, the package takes Homebrew's libomp
(`brew install libomp`) from `/opt/homebrew/opt/libomp` or
`/usr/local/opt/libomp`; without either, `find_package(otswap)` fails on
OpenMP.

The version given to `find_package` accepts the patch releases of that minor
version only: 0.4 accepts 0.4.x. The public structs (`Config`, `Result`,
`RealSpaceCatalog` and the others) keep their layout across them; code
compiled against the headers of another minor version must be rebuilt.

## From another nanobind module

The Python types of otswap are the C++ structs themselves, so another
nanobind module can take and return them, naming the C++ types in its own
signatures: `Config`, `CorrectionConfig`, `RedshiftCut`, `DistanceCut`,
`BoxLimits`, `Mask`, `DistanceTable`, `MpsProfile`, `BiasTable`, `Result` and
`RealSpaceCatalog`. For example, `double f (const otswap::Result&)` takes an
`otswap.Result`, and a function returning an `otswap::Config` returns an
`otswap.Config`. The other module binds none of these types itself.

The module includes otswap's headers and links the library, which come from
an installed otswap ([From C++ with CMake](#from-c-with-cmake)) or a build of
it from source: the wheel holds no otswap header, library or CMake file. With
CMake, nanobind found as otswap's own bindings find it:

```cmake
find_package(otswap 0.4 CONFIG REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(my_module NB_STATIC my_module.cpp)
target_link_libraries(my_module PRIVATE otswap::otswap)
```

At run time the module uses the otswap wheel of the same minor version.

- **The same nanobind.** Both modules are built with nanobind 3.1.0 and the
  same nanobind ABI (the same compiler family and C++ standard library), and
  neither defines `NB_DOMAIN`.
- **Import otswap first.** `import otswap` registers the types, so it comes
  before the other module's functions are used. The pairs `(min, max)` and
  `((x0, y0, z0), (x1, y1, z1))` are then converted in the other module's
  functions too.
- **The same minor version.** The headers the module is built against and
  the otswap it runs with share the layouts of the structs, which hold
  across the patch releases of a minor version only (see the C++ section): a
  module built against the headers of 0.3 must be rebuilt for 0.4.
- **The log stream.** A `Config` or `CorrectionConfig` received from Python
  holds `log = &std::clog`; the other module sets it before using it.
- **Angles.** A `Result` or `RealSpaceCatalog` made by the other module is
  in radians: its `angle_unit` is `"rad"`. The unit of one made by otswap
  from angles in degrees is kept by otswap, not in the object.
- **Results changed in place.** Take a `Result` as `const otswap::Result&`.
  One changed in place behind otswap's back keeps the degree arrays already
  read from it; to filter in place, call otswap's own functions.
- **Errors.** An `otswap::Error` thrown by otswap's Python functions reaches
  Python as `otswap.Error`. One thrown where the other module calls otswap in
  C++ is that module's to translate: otswap's translator does not recognise
  it with every standard library (with libc++, two modules that each link
  otswap do not share the type information of `otswap::Error`), and it then
  arrives as a plain `RuntimeError`. Catch it where the module calls otswap
  and rethrow the module's own error, its message prefixed by `"otswap: "`,
  as otswap does with the errors of meshsearch.
- **Malformed objects.** For a malformed `Result` or catalogue, which only
  another module can make, the library's consistency messages name the C++
  members (`meanDisplacement`), and the array properties raise
  `otswap.Error` naming the Python ones (`mean_displacement`).

## Linking the same OpenMP runtime

A process on macOS cannot hold two copies of the LLVM OpenMP runtime: the
second to start aborts it with `OMP: Error #15`. A Python package that
compiles OpenMP code of its own, imported in the same process as otswap,
must therefore load the runtime otswap loads.

The macOS wheels of otswap ship it in the `openmp` directory of the package,
the directory {py:func}`otswap.openmp_root` returns:

- `lib/libomp.dylib`, with install name `@rpath/libomp.dylib`; the module
  `_otswap` finds it through its rpath `@loader_path/openmp/lib`;
- `include/omp.h` and `include/ompx.h`;
- `LICENSE.TXT`, LLVM's licence: the macOS wheels include LLVM's libomp.

The runtime is the LLVM libomp otswap was built with; its version may change
between minor releases of otswap. On Linux, and in a build from source,
`otswap.openmp_root()` returns None and otswap uses the runtime it was built
against.

A package built with CMake links the same file as follows.

1. It requires otswap at build time as well as at run time, in the same
   minor version, so that the directory is there when CMake runs:
   `otswap >= 0.4, < 0.5` in `build-system.requires` and in `dependencies`.
2. Before `find_package(OpenMP)`, it takes `OpenMP_ROOT` from otswap
   whenever otswap gives one, in a wheel build and in a build from source
   alike. `Python_EXECUTABLE` is set by scikit-build-core, and otherwise by
   `find_package(Python COMPONENTS Interpreter)`:

   ```cmake
   execute_process(
     COMMAND "${Python_EXECUTABLE}" -c "import otswap; print(otswap.openmp_root() or '')"
     OUTPUT_VARIABLE otswap_openmp_root OUTPUT_STRIP_TRAILING_WHITESPACE)
   if (otswap_openmp_root)
     set(OpenMP_ROOT "${otswap_openmp_root}")
   endif()
   find_package(OpenMP REQUIRED COMPONENTS CXX)
   ```

   Its extension then records `@rpath/libomp.dylib`.
3. It gives its extension the rpath of otswap's directory seen from its own,
   for two packages installed side by side in one `site-packages`:

   ```cmake
   if (APPLE)
     set_target_properties(my_extension PROPERTIES
       INSTALL_RPATH "@loader_path/../otswap/openmp/lib")
   endif()
   ```

   The path holds for an extension placed directly in its package
   directory, one level below `site-packages`; each level deeper adds one
   `../`. With otswap installed elsewhere (another `--target`, an editable
   install), that rpath finds nothing.
4. Its wheel repair leaves the file where it is: with delocate,
   `delocate-wheel --exclude libomp`. A copy in the package's own `.dylibs/`
   is a second runtime.
5. It imports otswap before its own extension.

On Linux none of this is needed: `otswap.openmp_root()` returns None, so
step 2 leaves `OpenMP_ROOT` alone, and step 3 applies to macOS only.

This covers packages that follow it. Another package that ships its own
libomp, and is imported in the same process, can still clash with otswap and
with them.

`KMP_DUPLICATE_LIB_OK=TRUE` silences the error but keeps two runtimes in one
process, which can crash or give wrong results: it is not a fix.
