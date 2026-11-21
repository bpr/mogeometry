# mogeometry

Geometry processing on surface meshes and point clouds, in Mojo: a port
of [geometry-central](https://github.com/nmwsharp/geometry-central), the
C++ library behind Keenan Crane's discrete differential geometry course.
Halfedge meshes, their intrinsic and embedded geometry, and the
algorithms built on them: geodesic distance by the heat method, vector
transport, direction fields, intrinsic triangulations,
parameterization, remeshing.

Sparse linear algebra (matrices, sparse Cholesky and LDLᵀ, Krylov and
eigenvalue solvers) comes from [molina](https://github.com/bpr/molina),
the linear algebra of MojoDEC, a sibling project on discrete exterior
calculus.

**Status:** started. Nothing is ported yet; `docs/plan.md` gives the
order.

## Setup

```sh
pixi install      # also builds molina from git
pixi run test
```

Linux x86-64 only, as MojoDEC. Needs pixi 0.84 or newer (`pixi
self-update`) to build molina. The Mojo version is pinned to the one
MojoDEC's `pixi.lock` records.

## License

MIT (`LICENSE`). geometry-central is MIT too; its notice is in `NOTICE`.
