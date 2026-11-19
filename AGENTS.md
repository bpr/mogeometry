# AGENTS.md

Project context for AI coding assistants working on mogeometry. Tool-specific
workflow instructions live in each tool's own file, e.g.
CLAUDE.md.

## What this is

mogeometry is a pure Mojo port of geometry-central
(https://github.com/nmwsharp/geometry-central, MIT), the C++ geometry
processing library used in Keenan Crane's discrete differential geometry
course: halfedge surface meshes, their intrinsic and embedded geometry,
and algorithms on them (heat method distance, vector heat method,
direction fields, intrinsic triangulations, parameterization, remeshing,
point clouds). The port starts from geometry-central commit `019669d`
(2026-06-13; full hash in `NOTICE`).

It runs on Mojo 1.3 and builds on MojoDEC (a sibling checkout,
`../MojoDEC`), a discrete exterior calculus library by the same author,
for sparse linear algebra. Much of the course's material is DEC, so the
two overlap (the cotan Laplacian, dual areas, `d0`/`d1` and the Hodge
stars); mogeometry works on halfedge meshes, MojoDEC on simplicial
complexes of any dimension.

Status: the start of layer 1 is ported (`Vector2`, `Vector3`, elementary
geometry, the scalar helpers of `utilities.h`), with the reference-data
pipeline that tests it against geometry-central. `docs/plan.md` gives
the order and the design questions to settle at layer 2.

## Commands

All tooling runs through pixi (linux-64 only):

```sh
pixi install                 # Mojo 1.3 nightly (MojoDEC's version) + MAX, Python
pixi run deps                # precompile NuMojo and MojoDEC -> build/*.mojoc
pixi run test                # every tests/test_*.mojo (runs deps first)
pixi run test-release        # tests with -D ASSERT=none
pixi run package             # build/mogeometry.mojoc
pixi run format              # mojo format mogeometry tests
pixi run doc-check           # doc string warnings vs docs/docstring-baseline.txt
pixi run -e ref ref-data     # rebuild tests/data/*/reference.txt from geometry-central
```

The `ref` environment has a C++ compiler and Eigen 3.3.9 (the version
geometry-central's CMake downloads), and no Mojo. It compiles each
`tests/data/<area>/reference.cpp` against geometry-central at
`$GEOMETRY_CENTRAL_PATH` (default `../../c++/geometry-central`, checked
out at the commit in `NOTICE`); see `tests/data/README.md`.

Run a single test file from the repo root (after `pixi run deps`):

```sh
pixi run mojo run -I . -I build tests/test_build.mojo
```

Property tests (`tests/test_properties*.mojo`, none yet) need proptest at
`../proptest`; set `PROPTEST_PATH` to use another checkout.

## Project layout

Packages follow geometry-central's directories, so a C++ file maps to one
Mojo module of the same name (`src/surface/heat_method_distance.cpp` ->
`mogeometry/surface/heat_method_distance.mojo`):

- `mogeometry/utilities`: vectors, quaternions, elementary geometry,
  disjoint sets, nearest neighbors.
- `mogeometry/surface`: meshes, mesh data, geometry, algorithms.
- `mogeometry/pointcloud`: point clouds.
- `mogeometry/__init__.mojo` re-exports the public API (add new public
  names there).
- geometry-central's `numerical/` has no counterpart: use MojoDEC.
- `tests/`: one `test_*.mojo` per area; reference data from
  geometry-central in `tests/data/<area>/` with the C++ driver that
  wrote it; meshes in `tests/meshes/`.
- `tools/build_deps.sh`: the `deps` task.
- `docs/`: the plan, the doc string guide, the doc-check baseline.

## Tech stack and toolchain

- **MojoDEC dependency.** `pixi run deps` precompiles MojoDEC's vendored
  NuMojo and MojoDEC itself from `$MOJODEC_PATH` (default `../MojoDEC`)
  into `build/numojo.mojoc` and `build/mojodec.mojoc`, with this
  workspace's Mojo (a `.mojoc` loads only in the Mojo that built it), so
  code compiles with `-I . -I build`. Import from `mojodec` directly
  (`from mojodec import SparseCholesky`, `from mojodec.core.sparse import
  CSRMatrix`). Read `../MojoDEC/AGENTS.md` before using its solvers:
  `SparseLDLT` (refinement steps, supernodal layout), `SparseCholesky`
  (SPD, no refinement), `cg`/`minres` (keyword options, `nullspace=`),
  `eigsh`, and NuMojo's slow elementwise arithmetic are all described
  there. A fix MojoDEC needs goes into MojoDEC.
- **Pinned Mojo.** `pixi.toml` pins `mojo` and `max-core` to the exact
  nightly MojoDEC's `pixi.lock` records (1.3.0.dev2026110305), because
  MojoDEC and NuMojo are compiled from source here and newer nightlies
  break them (2026111805 removed `std.memory.UnsafePointer`). When
  MojoDEC moves to a new nightly, move the pins with it.
- `-D ASSERT=none` removes `debug_assert` bounds checks. `mojo run`/`build`
  accept it (the `*-release` tasks pass it via `MOJO_FLAGS`); `mojo
  precompile` does not, so it applies only to code compiled from source.

### Known toolchain bugs

These are from MojoDEC's experience with the same toolchain; its
AGENTS.md has the details and more.

- Write floats as text with `write_float` (MojoDEC's
  `routines/io/scan.mojo`), not
  `String(x)`: Mojo prints some doubles with a digit too few, so they
  don't read back. Parse integer text with `atol`, not `atof`.
- `assert_almost_equal` has a default `rtol=1e-5` on top of `atol`, so
  every call passes `rtol` explicitly (normally `rtol=0`);
  `tests/run_all.sh` fails on a call without one.
- Write `sqrt(x)` (`std.math`), never `x ** 0.5`: Mojo's `**` with a
  fractional exponent is off by up to ~5e-12 relative. For other powers
  use `libm_pow` (`mogeometry.utilities.utilities`), the C library's.
- Mojo's `std.math.exp` on `Float64` is off by up to ~35,000 ulps
  (~8e-12 relative; measured against libm on 200,000 points in
  [-10, 10]); use `libm_exp`. `log`, `sin`, `cos`, `atan2`, `acos`, and
  `hypot` are within an ulp of libm.
- Mojo's float operations allow contraction, so LLVM may fuse `a * b + c`
  into one rounding, or not, depending on the code around it. Where two
  code paths must agree bit for bit, write explicit `x.fma(y, acc)`.
  It is why the port does not match geometry-central bit for bit: even
  `Vector3.norm2()` differs in the last bit (the reference is built with
  `-ffp-contract=off`). Tests compare in roundoffs instead (Testing).
- Mojo does not merge overloads of a name imported from two modules
  (`from .a import f` and `from .b import f` is an "ambiguous import"
  error), so a function geometry-central overloads across headers (a
  free `dot` for `Vector2` and for `Vector3`) becomes a method.
- Mojo's `>>` by 64 or more is undefined (it gave garbage): bound shift
  loops.
- NuMojo's elementwise array arithmetic is slow (50 times a sparse
  product for 4,000 entries); in hot loops, write the loop over
  `unsafe_get`/`unsafe_set`.

## Architecture

`docs/plan.md` lists the layers in porting order and the design
questions to settle before the mesh (element handles, containers that
follow mesh mutation, cached geometry quantities, implicit twins).
Record each decision here, with its reason, when it is made.

- **Small vectors are three (or two) `Float64` fields.** `Vector2` and
  `Vector3` are `TrivialRegisterPassable` structs with fields `x, y(, z)`,
  as geometry-central's are; this keeps its order of operations. It is
  provisional: a `SIMD[DType.float64, 4]` layout may be much faster, and
  `docs/plan.md` (Performance experiments) records the comparison to
  run. A struct's `__getattr__(self, name: StringLiteral)` lets a SIMD
  layout keep the `v.x` syntax, so the swap can stay inside the module.
- **Vector operations are methods, not free functions:** `u.dot(v)`,
  `u.cross(v)`, `v.norm()`, `v.unit()` (Mojo cannot import two free
  `dot`s; see Known toolchain bugs). Geometry code ported later follows
  the same rule.

## Porting

- geometry-central is MIT: port its code directly. Each ported file's
  header names the geometry-central file(s) it follows and the license
  (`NOTICE` has the full text and commit), as MojoDEC's ports of faer
  do. Its bundled dependencies are permissive too (happly MIT, nanoflann
  BSD, nanort MIT, Eigen MPL2); credit any that a port follows.
- Keep geometry-central's algorithms and numerics (the same formulas,
  the same order of operations where results are compared to the
  reference data), but write Mojo, not transliterated C++: snake_case
  names (`compute_distance`, `vertex_dual_areas`), structs instead of
  class hierarchies (traits where geometry-central has virtual
  interfaces), raising errors instead of exceptions or silent garbage.
  Record each renaming that is not mechanical in the module's doc
  string.
- Parallelize only after the serial port matches the reference, and
  then with MojoDEC's `parallel_ranges` (`mojodec.core.parallel`),
  keeping results bit for bit serial (see MojoDEC's AGENTS.md, Parallel
  loops). geometry-central itself is serial.

## Coding standards

- Use American spelling everywhere (code, comments, doc strings, docs,
  commit messages): "-ize", not "-ise" (`normalize`, `initialize`).
- Return a small struct with named fields rather than a tuple, and an
  enum-like struct with `comptime` members rather than a Bool or String
  flag.
- Public symbols get Google-style doc strings (`Args:`, `Parameters:`,
  `Returns:`, `Raises:`; see `docs/writing-docs.md`). The doc-check
  warning count must not get worse than `docs/docstring-baseline.txt`.

## Testing

- Each test file ends with
  `TestSuite.discover_tests[__functions_in_module()]().run()`, so any
  `def test_*() raises` in the file is picked up automatically.
- Tests compare against geometry-central's own results (its GoogleTest
  cases, and reference data its C++ writes), and against MojoDEC where
  both compute the same quantity.
- Floating-point checks use `assert_almost_equal(..., atol=..., rtol=0)`
  with a tolerance from the problem's conditioning, not a loose default.
- Reference data stores each double as its bit pattern, so the expected
  values are exactly geometry-central's. The port differs from them in
  the last bits (FMA contraction), so `tests/test_utilities.mojo`
  compares in roundoffs, `eps * max(1, |value|)`, and prints the worst
  difference per family. The drivers keep random inputs well-conditioned
  (no nearly parallel lines, no flat tetrahedra) so the tolerance can be
  tight (16 roundoffs); degenerate input gets hand-written tests.
