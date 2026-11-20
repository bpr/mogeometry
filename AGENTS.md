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
pixi run format              # mojo format mogeometry tests benchmarks
pixi run doc-check           # doc string warnings vs docs/docstring-baseline.txt
pixi run -e ref ref-data     # rebuild tests/data/*/reference.txt from geometry-central
pixi run bench-handles       # the element-handle experiment (release build)
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
- `benchmarks/`: timing experiments behind design decisions (each
  file's doc string says what it measures; AGENTS.md, Architecture,
  records the results).
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
- **Element handles are bare typed indices; the mesh navigates.**
  `Vertex`, `Halfedge`, `Edge`, `Face` will be one-`Int` structs
  (distinct types, so a vertex cannot be passed as a face), and
  navigation goes through the mesh: `mesh.next(h)`, `mesh.twin(h)`,
  `mesh.vertex(h)`, `for h in mesh.outgoing_halfedges(v)`.
  geometry-central's handles carry a mesh pointer (`h.twin().vertex()`);
  `benchmarks/handles.mojo` (`pixi run bench-handles`) measured both
  styles in Mojo, with safe `Pointer`s, on a torus of 1M faces
  (i7-10875H; 2026-11-18):
  - Bare indices ran as fast as raw array code in every kernel (one-ring
    sums, face areas, cotan weights, navigation only). Pointer handles
    were 7 to 10% slower on the one-ring position sum (in grid order,
    repeatably) and equal elsewhere. With shuffled element order, all
    variants were memory-bound and within run-to-run noise (+-25%).
  - Iterators (`for h in mesh.outgoing(v)`) cost nothing over a manual
    `do ... while` loop, in either style.
  - A pointer handle is 16 bytes, an index 8: stored handles (paths,
    surface points, lists of elements) would double.
  - A pointer handle's type carries the mesh's origin
    (`Halfedge[o: ImmOrigin]`), so every struct and function that holds
    or passes handles would take that parameter, and a struct that owns
    a mesh could not hold handles into it (it would refer to itself).
  - Mojo does not stop mutation while pointer handles are alive (its
    exclusivity checks are per call, not across statements), so they
    would not make mutation safer; it does reject using one after the
    mesh is moved.
  The cost is wordier navigation (`mesh.vertex(mesh.twin(h))`), eased
  by iterators and by the mesh's compound queries (geometry-central's
  `tipVertex`, `degree`, `isBoundary` become `mesh.tip_vertex(h)`,
  `mesh.degree(v)`, `mesh.is_boundary(v)`). An index alone cannot tell
  which mesh it came from; the next two items are how that is checked.
- **Compile-time branding was tried and rejected.** A handle can carry
  its mesh's origin as a zero-size parameter (`Halfedge[o: ImmOrigin]`,
  still 8 bytes; mesh methods `def next[o: ImmOrigin](ref [o] self, h:
  Halfedge[o])`), and then a handle from one mesh passed to another is
  a compile error. But Mojo treats the brand as a reference into the
  mesh, so a mutating method cannot take a handle: `mesh.flip(e)` fails
  with "aliasing values passed mutably to 'self' and immutably to 'e'".
  The brand also spreads the origin parameter, as pointer handles do.
  (Checked on 1.3.0.dev2026110305.)
- **Debug builds tag handles with their mesh.** Each mesh gets a tag
  when it is built (copies keep it: they have the same elements). In
  debug builds, the handles a mesh returns carry the tag in the high
  bits of their `Int` (indices below 2^44 leave 20 bits), and every
  mesh access checks the tag, the index range, and that the element is
  not deleted. `h.index` masks the tag off. A handle made from a plain
  integer (`Vertex(5)`) has tag 0 and is not checked. Part of the tag
  is a generation that `compress()` bumps, so handles from before a
  renumbering fail the check. With `-D ASSERT=none` the tag is compiled
  out: handles are plain indices, at no cost.
- **Mesh data containers sync explicitly.** Mutations (flip, split,
  collapse) never renumber: new elements go at the end, deleted ones
  are marked. Only `compress()` renumbers, and the mesh logs each
  renumbering with the generation it starts. Containers hold no
  reference to the mesh (geometry-central's `MeshData` holds a pointer
  and registers callbacks); they record the generation they match, and
  `data.sync(mesh)` grows them and applies the pending renumberings. In
  debug builds `data[v]` checks `v`'s generation against the
  container's, so a container not synced after a `compress()` raises
  instead of returning wrong values. This keeps containers free of an
  origin parameter, lets a struct own a mesh together with its data,
  and costs nothing per access in release builds; the price is a
  `sync` call after a phase of mutation.
- **The manifold mesh only, for now.** The mesh is geometry-central's
  `ManifoldSurfaceMesh`: implicit twins (`twin(h) = h ^ 1`,
  `edge(h) = h / 2`), no sibling arrays, no runtime twin flag (its
  general `SurfaceMesh` branches on `useImplicitTwinFlag` in every
  `twin`). Most of the course's algorithms need no more. Algorithms
  navigate only through mesh methods, so a general mesh can be added
  later as a second type with the same methods.
- **Connectivity is stored as `Int32`; handles stay `Int`.** The
  storage type is one module-level constant (`comptime MESH_INDEX =
  DType.int32`), private behind the mesh methods, which widen every
  read to `Int`; switching all meshes to 64-bit is a one-line change.
  `benchmarks/handles.mojo` also runs the raw and index kernels on a
  copy of the mesh whose arrays hold `Int32`. Three runs on the
  1M-face torus (i7-10875H; 2026-11-18), grid order: one-ring sums 20
  to 25% faster, navigation only 25 to 40%, face areas 5 to 20%, cotan
  weights 5 to 8%; on a 20,000-face torus, which fits in cache, 10 to
  15% on the one-ring kernels. Shuffled, `Int32` was usually about 10%
  faster but within run-to-run noise (+-25%). The connectivity arrays
  drop from 144 to 72 MB per 1M faces; the limit is 2^31 halfedges
  (about 700M faces). A per-mesh width parameter
  (`SurfaceMesh[I: DType = DType.int32]`) was considered: functions
  taking a bare `SurfaceMesh` are generic over it for free, but every
  struct that holds a mesh, and every return type, would have to bind
  or carry the width, and the unused width would need its own tests.
  Since the width is private, adding the parameter later touches only
  those structs.
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
