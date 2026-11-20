# Plan

The port follows geometry-central's layers, bottom up, and each layer is
tested before the next one starts. Line counts are geometry-central's
(`src/` plus its headers' inline code), to size the work; the whole
library is about 40,000 lines of C++.

## Order

1. **Utilities** (`utilities/`): `Vector2` (also complex arithmetic,
   as geometry-central uses it for tangent vectors), `Vector3`,
   `Quaternion`, `elementary_geometry` (triangle and segment helpers),
   `DisjointSets`, k-nearest neighbors (geometry-central wraps nanoflann;
   a kd-tree in Mojo).
   *Done:* `Vector2`, `Vector3`, `elementary_geometry` except
   `triTriIntersection` (744 lines, `tri_tri_intersect.cpp`; no users
   before layer 5), and `utilities.h`'s `PI`, `clamp`,
   `regularizeAngle`. *Deferred until a user needs them:* `Quaternion`,
   `UnitVector3`, `DisjointSets`, the random-number helpers, and k-nearest
   neighbors (point clouds, layer 7; MojoDEC's `KDTree`,
   `mojodec.core.kd_tree`, may serve, after a check of its speed against
   nanoflann's).
2. **Surface mesh** (`surface/surface_mesh`, `manifold_surface_mesh`,
   `halfedge_element_types`, `surface_mesh_factories`, `meshio`):
   the halfedge arrays, element handles and their navigation,
   iteration (`vertices()`, `outgoingHalfedges()`, ...), the data
   containers (`MeshData<E, T>`: `VertexData`, `FaceData`, ...),
   construction from a polygon soup, `SimplePolygonMesh`, OBJ/PLY/OFF/STL
   input and output, then the mutations (`flip`, `splitEdge`,
   `collapseEdge`, `insertVertex`, ...) and `MutationManager`.
   *Done:* `ManifoldSurfaceMesh` from face lists (exactly
   geometry-central's element numbering), the handles with debug mesh
   tags, navigation, element ranges and navigation iterators, the
   whole-mesh queries (counts, Euler characteristic, genus, components,
   index maps, `validateConnectivity`), and `MeshData` with growth-only
   `sync`. *Next:* mesh input and output (OBJ, PLY, OFF, STL; tests then
   use files), `SimplePolygonMesh`, the mutations and `compress()` with
   its renumbering log for `sync`, and the construction from a twin
   list. *Deferred:* the general `SurfaceMesh` (decided: manifold only,
   for now), `MeshData` arithmetic operators, `RichSurfaceMeshData`.
3. **Geometry** (`base_geometry_interface`, `intrinsic_`,
   `extrinsic_`, `embedded_geometry_interface`,
   `edge_length_geometry`, `vertex_position_geometry`): edge lengths,
   angles, areas, cotan weights, curvatures, tangent bases and
   transport, the Laplacians and mass matrices (cotan,
   connection, Crouzeix-Raviart, polygon), and the DEC operators
   (`d0`, `d1`, `hodge0..2`).
4. **Algorithms on fixed meshes**: heat method distance, vector heat
   method, signed heat method, direction fields, fast marching, exact
   geodesics (MMP), trace geodesic and surface points, parameterization
   (BFF, `parameterize`), stripe patterns, tufted Laplacian,
   uniformization, marching triangles, mesh graph algorithms, ray
   tracing.
5. **Mesh-changing algorithms**: subdivision, remeshing, quadric error
   simplification, surgery, intrinsic mollification.
6. **Intrinsic triangulations**: simple IDT, signpost, integer
   coordinates (normal coordinates), common subdivision, flip
   geodesics, geodesic centroidal Voronoi tessellation, embed convex.
7. **Point clouds** (`pointcloud/`): neighborhoods, local
   triangulations, the point cloud heat solver, sampling.

`numerical/` (Eigen and SuiteSparse wrappers) is not ported: its users
call MojoDEC's `CSRMatrix`, `SparseLDLT`/`SparseCholesky`, `cg`,
`minres`, and `eigsh`. Complex-valued systems (the connection
Laplacians, direction fields, the vector heat method) need complex
sparse matrices or their real 2n x 2n form; geometry-central uses
`std::complex<double>`. Decide when layer 3 needs it: the real form
works with MojoDEC as it is.

## Performance experiments

Comparisons to run before a choice hardens. Each records what to
measure and where; the result and the decision go in AGENTS.md
(Architecture).

- **`Vector3` as `SIMD[DType.float64, 4]`.** `Vector3` is three
  `Float64` fields (AGENTS.md, Architecture), which keeps
  geometry-central's operations in its order. A 4-lane SIMD layout (the
  fourth lane zero) could make arithmetic, `dot`, `cross`, and `norm`
  one or two vector instructions each, which may be a big win in the
  geometry layer's loops. To measure:
  - Write `Vector3Simd` with the same API in a scratch module.
    `__getattr__(self, name: StringLiteral)` keeps `v.x` working
    (checked on this toolchain), and `__setattr__` would have to keep
    `v.x = ...` working; check in the generated assembly that both cost
    nothing.
  - Microbenchmark the kernels that dominate layer 3, on a mesh with
    about 1M faces: face normals and areas, cotan weights, assembly of
    the cotan Laplacian, and a gather and scatter over vertex positions.
    Use real data access patterns: positions are read through face
    indices, not sequentially.
  - Compare storage too. `List[Vector3]` is 24 bytes per entry with
    fields, 32 with SIMD (a third more memory traffic); also try a
    structure of arrays (separate x, y, z lists) for the hot loops.
  - Check accuracy against `tests/data/utilities/reference.txt`. SIMD
    reductions (`reduce_add` in `dot` and `norm2`) sum in a different
    order, so results move by a few roundoffs; the test's 16-roundoff
    tolerance should still hold, and the worst-case numbers it prints
    show by how much.
  - When: the microbenchmark can run now. The decision is best made
    once the first layer 3 kernels exist, so it is measured on code that
    will ship. If the swap stays inside `vector3.mojo` (the `__getattr__`
    route), deciding later costs little; if it needs call sites to
    change, decide before layer 3 grows.

- **Element order.** The handle experiment found memory order matters
  far more than handle style: on a 1M-face torus, the same kernels run
  15 to 20 times slower when vertices and faces are randomly relabeled
  than in grid order (`pixi run bench-handles`). geometry-central keeps
  the order of the input file. Worth measuring, once meshes load from
  files: how ordered real meshes are (the test meshes, Thingi10k-style
  scans), and what an optional reordering at load costs and saves
  (vertices along a space-filling curve or by reverse Cuthill-McKee,
  faces sorted by their lowest vertex, halfedges by face). It must
  stay optional and explicit: indices are user-visible, and results
  compared with geometry-central's need its order.

## Design questions to settle at layer 2

- **Element handles.** *Decided: bare typed indices, navigation on the
  mesh* (AGENTS.md, Architecture, has the reasons and the numbers from
  `benchmarks/handles.mojo`). geometry-central's `Vertex`, `Halfedge`,
  ... carry a pointer to their mesh, so `v.halfedge().twin().vertex()`
  works. Measured in Mojo, that style is no faster, doubles the size of
  a stored handle, and puts an origin parameter on every type that
  holds a handle.
- **Checking a handle's mesh.** *Decided: a mesh tag in debug builds*
  (AGENTS.md, Architecture). Compile-time branding was tried and fails
  on mutating calls.
- **Containers that follow mutation.** *Decided: explicit sync*
  (AGENTS.md, Architecture). geometry-central's `MeshData` registers
  callbacks with its mesh so it grows, shrinks, and permutes when the
  mesh changes; here containers hold no reference to the mesh and
  `data.sync(mesh)` catches up, with a generation check in debug builds.
- **Index width in storage.** *Decided: `Int32`, set by one constant*
  (AGENTS.md, Architecture). Handles stay `Int`; navigation is 20 to
  40% faster on a 1M-face mesh in grid order, 5 to 20% on face and
  edge kernels.
- **Cached geometry.** geometry-central computes quantities on demand
  (`requireEdgeLengths()`, then `geometry.edgeLengths`) and keeps them
  until `unrequire`. MojoDEC caches nothing; here the cache is the
  library's API, so port it, with the dependency graph between
  quantities explicit.
- **Implicit twins.** *Decided: the manifold mesh only, for now*
  (AGENTS.md, Architecture). `ManifoldSurfaceMesh` stores no twin array
  (the twin of halfedge `h` is `h ^ 1`); the general `SurfaceMesh`
  stores siblings. The navigation API stays one a general mesh could
  implement later.

## Testing

- geometry-central's GoogleTest suite (`test/src/`, 16 files) is the
  starting point: port each test with the code it covers.
- Reference values come from geometry-central itself: a small C++
  driver per area, built in the pixi `ref` environment (a C++ compiler
  and Eigen 3.3.9), writes the expected results to `tests/data/<area>/`
  (as MojoDEC does for faer's AMD and amgcl). Commit the data and the
  driver, not the build. `tests/data/README.md` describes the format;
  `tests/data/utilities/` is the first. Layers that need geometry-central
  compiled (not just its headers) list the sources to build in
  `tests/data/build_reference.sh`; once that list gets long, build the
  library with CMake instead.
- Check the license of each of geometry-central's test meshes
  (`test/assets/`) before copying it in; Keenan Crane's models (spot,
  bob) are CC0.
- Where a quantity also exists in MojoDEC (the DEC operators, cotan
  Laplacian, dual areas), compare the two on the same mesh, minding
  orientation and sign conventions.
