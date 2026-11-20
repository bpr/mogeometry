# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# The connectivity arrays follow geometry-central's ManifoldSurfaceMesh
# (include/geometrycentral/surface/surface_mesh.h, commit 019669d, MIT
# License, Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
# contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""The element-handle experiment (docs/plan.md, Design questions): what
geometry-central's chained navigation style costs in Mojo.

A minimal manifold triangle mesh with implicit twins (`twin(h) = h ^ 1`,
`edge(h) = h / 2`), and the same four kernels written five ways:

- `raw`: indices into the arrays, no handle types (the lower bound).
- `index`: handles are bare typed indices (`Halfedge(index)`), and the
  mesh navigates: `mesh.vertex(mesh.twin(h))`.
- `index-iter`: the same, iterating with `for h in mesh.outgoing(v)`.
- `pointer`: handles carry a safe `Pointer` to the mesh, as
  geometry-central's carry a `SurfaceMesh*`: `h.twin().vertex()`.
- `pointer-iter`: the same, iterating with `for h in v.outgoing()`.

Each variant without pointer handles also runs on a copy of the mesh
whose connectivity arrays hold `Int32` instead of `Int` (`raw32`,
`index32`, `index-iter32`): handles stay `Int`, only the storage
narrows, which halves the memory traffic of navigation.

Kernels: `ring` sums neighbor positions over each vertex's one-ring (a
Laplacian smoothing step); `connect` sums neighbor and face indices over
each one-ring (navigation only, no geometry); `faces` computes each face's
area from its halfedge loop; `cotan` computes each edge's cotan weight
from the two corners opposite it.

The mesh is a torus grid, in grid order (good locality) and with its
vertices and faces randomly relabeled (the order of a mesh read from a
file that was not optimized). Every variant's checksum is printed, so
they can be seen to compute the same thing.

```sh
pixi run bench-handles
```
"""

from std.collections import Dict
from std.math import cos, sin, sqrt
from std.random import random_ui64, seed
from std.sys import size_of
from std.time import perf_counter_ns

from mogeometry import PI, Vector3

comptime RUNS = 9
"""Each timing is the best of this many runs."""


# ===----------------------------------------------------------------------=== #
# The mesh
# ===----------------------------------------------------------------------=== #


struct Mesh[I: DType](Movable):
    """The manifold connectivity arrays of geometry-central, implicit twins.

    Parameters:
        I: The integer type the connectivity arrays store.
    """

    var he_next: List[Scalar[Self.I]]
    var he_vertex: List[Scalar[Self.I]]  # the halfedge's tail
    var he_face: List[Scalar[Self.I]]
    var v_halfedge: List[Scalar[Self.I]]  # an outgoing halfedge
    var f_halfedge: List[Scalar[Self.I]]
    var positions: List[Vector3]
    var n_vertices: Int
    var n_faces: Int
    var n_edges: Int

    def __init__(
        out self,
        n_vertices: Int,
        triangles: List[Int],
        var positions: List[Vector3],
    ) raises:
        # pair the two halfedges of each edge: the first seen is 2e, its
        # twin 2e + 1
        self.n_vertices = n_vertices
        self.n_faces = len(triangles) // 3
        var n_he = 3 * self.n_faces
        self.he_next = List[Scalar[Self.I]](length=n_he, fill=-1)
        self.he_vertex = List[Scalar[Self.I]](length=n_he, fill=-1)
        self.he_face = List[Scalar[Self.I]](length=n_he, fill=-1)
        self.v_halfedge = List[Scalar[Self.I]](length=n_vertices, fill=-1)
        self.f_halfedge = List[Scalar[Self.I]](length=self.n_faces, fill=-1)
        self.positions = positions^
        var edge_of = Dict[Int, Int]()
        var n_edges = 0
        for f in range(self.n_faces):
            var hs = SIMD[DType.int64, 4](0)
            for k in range(3):
                var a = triangles[3 * f + k]
                var b = triangles[3 * f + (k + 1) % 3]
                var key = min(a, b) * n_vertices + max(a, b)
                var h: Int
                if key in edge_of:
                    h = 2 * edge_of[key] + 1
                else:
                    edge_of[key] = n_edges
                    h = 2 * n_edges
                    n_edges += 1
                if self.he_vertex[h] != -1:
                    raise Error("not manifold or not oriented")
                self.he_vertex[h] = Scalar[Self.I](a)
                self.he_face[h] = Scalar[Self.I](f)
                self.v_halfedge[a] = Scalar[Self.I](h)
                hs[k] = Int64(h)
            for k in range(3):
                self.he_next[Int(hs[k])] = Scalar[Self.I](hs[(k + 1) % 3])
            self.f_halfedge[f] = Scalar[Self.I](hs[0])
        self.n_edges = n_edges
        if 2 * n_edges != n_he:
            raise Error("the mesh has boundary")

    def __init__(out self, *, narrowing: Mesh[_]):
        """A copy of `other` whose arrays store `I`."""
        self.he_next = _narrow[Self.I](narrowing.he_next)
        self.he_vertex = _narrow[Self.I](narrowing.he_vertex)
        self.he_face = _narrow[Self.I](narrowing.he_face)
        self.v_halfedge = _narrow[Self.I](narrowing.v_halfedge)
        self.f_halfedge = _narrow[Self.I](narrowing.f_halfedge)
        self.positions = narrowing.positions.copy()
        self.n_vertices = narrowing.n_vertices
        self.n_faces = narrowing.n_faces
        self.n_edges = narrowing.n_edges


comptime Mesh64 = Mesh[DType.int64]
comptime Mesh32 = Mesh[DType.int32]


def _narrow[J: DType](xs: List[Scalar[_]]) -> List[Scalar[J]]:
    var ys = List[Scalar[J]](capacity=len(xs))
    for x in xs:
        ys.append(Scalar[J](x))
    return ys^


def torus(n: Int, m: Int, shuffle: Bool) raises -> Mesh64:
    """An `n x m` torus grid, each quad split in two."""
    var nv = n * m
    var perm_v = _permutation(nv, shuffle)
    var pos = List[Vector3](length=nv, fill=Vector3.zero())
    for i in range(n):
        for j in range(m):
            var u = 2 * PI * Float64(i) / Float64(n)
            var w = 2 * PI * Float64(j) / Float64(m)
            pos[perm_v[i * m + j]] = Vector3(
                (2 + cos(w)) * cos(u), (2 + cos(w)) * sin(u), sin(w)
            )
    var quads = List[Int]()
    for i in range(n):
        for j in range(m):
            var a = perm_v[i * m + j]
            var b = perm_v[((i + 1) % n) * m + j]
            var c = perm_v[((i + 1) % n) * m + (j + 1) % m]
            var d = perm_v[i * m + (j + 1) % m]
            quads.append(a)
            quads.append(b)
            quads.append(c)
            quads.append(a)
            quads.append(c)
            quads.append(d)
    var nf = len(quads) // 3
    var perm_f = _permutation(nf, shuffle)
    var tris = List[Int](length=3 * nf, fill=0)
    for f in range(nf):
        for k in range(3):
            tris[3 * perm_f[f] + k] = quads[3 * f + k]
    return Mesh64(nv, tris, pos^)


def _permutation(n: Int, shuffle: Bool) -> List[Int]:
    var p = List[Int](capacity=n)
    for i in range(n):
        p.append(i)
    if shuffle:
        for i in range(n - 1, 0, -1):
            var j = Int(random_ui64(0, UInt64(i)))
            var t = p[i]
            p[i] = p[j]
            p[j] = t
    return p^


# ===----------------------------------------------------------------------=== #
# Handles that are bare indices; the mesh navigates
# ===----------------------------------------------------------------------=== #


@fieldwise_init
struct Vertex(Equatable, TrivialRegisterPassable):
    var index: Int


@fieldwise_init
struct Halfedge(Equatable, TrivialRegisterPassable):
    var index: Int


@fieldwise_init
struct Face(Equatable, TrivialRegisterPassable):
    var index: Int


struct Outgoing[origin: ImmOrigin, I: DType](Copyable, Iterator):
    """The outgoing halfedges of a vertex, counterclockwise."""

    comptime Element = Halfedge
    var mesh: Pointer[Mesh[Self.I], Self.origin]
    var current: Int
    var start: Int
    var started: Bool

    def __init__(
        out self, mesh: Pointer[Mesh[Self.I], Self.origin], start: Int
    ):
        self.mesh = mesh
        self.current = start
        self.start = start
        self.started = False

    def __iter__(ref self) -> Self:
        return self.copy()

    def __next__(mut self) raises StopIteration -> Halfedge:
        if self.started and self.current == self.start:
            raise StopIteration()
        self.started = True
        var h = self.current
        self.current = Int(self.mesh[].he_next[h ^ 1])
        return Halfedge(h)


@always_inline
def nav_next[I: DType](mesh: Mesh[I], h: Halfedge) -> Halfedge:
    return Halfedge(Int(mesh.he_next[h.index]))


@always_inline
def nav_twin[I: DType](mesh: Mesh[I], h: Halfedge) -> Halfedge:
    return Halfedge(h.index ^ 1)


@always_inline
def nav_vertex[I: DType](mesh: Mesh[I], h: Halfedge) -> Vertex:
    return Vertex(Int(mesh.he_vertex[h.index]))


@always_inline
def nav_face[I: DType](mesh: Mesh[I], h: Halfedge) -> Face:
    return Face(Int(mesh.he_face[h.index]))


@always_inline
def nav_halfedge[I: DType](mesh: Mesh[I], v: Vertex) -> Halfedge:
    return Halfedge(Int(mesh.v_halfedge[v.index]))


@always_inline
def nav_halfedge[I: DType](mesh: Mesh[I], f: Face) -> Halfedge:
    return Halfedge(Int(mesh.f_halfedge[f.index]))


@always_inline
def outgoing[
    I: DType
](ref mesh: Mesh[I], v: Vertex) -> Outgoing[origin_of(mesh), I]:
    return Outgoing(Pointer(to=mesh), Int(mesh.v_halfedge[v.index]))


# ===----------------------------------------------------------------------=== #
# Handles that carry a pointer to the mesh (geometry-central's style)
# ===----------------------------------------------------------------------=== #


@fieldwise_init
struct PVertex[origin: ImmOrigin](TrivialRegisterPassable):
    var index: Int
    var mesh: Pointer[Mesh64, Self.origin]

    @always_inline
    def halfedge(self) -> PHalfedge[Self.origin]:
        return PHalfedge(Int(self.mesh[].v_halfedge[self.index]), self.mesh)

    @always_inline
    def outgoing(self) -> POutgoing[Self.origin]:
        return POutgoing(self.mesh, Int(self.mesh[].v_halfedge[self.index]))


@fieldwise_init
struct PHalfedge[origin: ImmOrigin](TrivialRegisterPassable):
    var index: Int
    var mesh: Pointer[Mesh64, Self.origin]

    @always_inline
    def next(self) -> Self:
        return Self(Int(self.mesh[].he_next[self.index]), self.mesh)

    @always_inline
    def twin(self) -> Self:
        return Self(self.index ^ 1, self.mesh)

    @always_inline
    def vertex(self) -> PVertex[Self.origin]:
        return PVertex(Int(self.mesh[].he_vertex[self.index]), self.mesh)

    @always_inline
    def face(self) -> PFace[Self.origin]:
        return PFace(Int(self.mesh[].he_face[self.index]), self.mesh)

    @always_inline
    def __eq__(self, other: Self) -> Bool:
        return self.index == other.index


@fieldwise_init
struct PFace[origin: ImmOrigin](TrivialRegisterPassable):
    var index: Int
    var mesh: Pointer[Mesh64, Self.origin]

    @always_inline
    def halfedge(self) -> PHalfedge[Self.origin]:
        return PHalfedge(Int(self.mesh[].f_halfedge[self.index]), self.mesh)


struct POutgoing[origin: ImmOrigin](Copyable, Iterator):
    """The outgoing halfedges of a vertex, as pointer handles."""

    comptime Element = PHalfedge[Self.origin]
    var mesh: Pointer[Mesh64, Self.origin]
    var current: Int
    var start: Int
    var started: Bool

    def __init__(out self, mesh: Pointer[Mesh64, Self.origin], start: Int):
        self.mesh = mesh
        self.current = start
        self.start = start
        self.started = False

    def __iter__(ref self) -> Self:
        return self.copy()

    def __next__(mut self) raises StopIteration -> PHalfedge[Self.origin]:
        if self.started and self.current == self.start:
            raise StopIteration()
        self.started = True
        var h = self.current
        self.current = Int(self.mesh[].he_next[h ^ 1])
        return PHalfedge(h, self.mesh)


@always_inline
def pvertex(ref mesh: Mesh64, i: Int) -> PVertex[origin_of(mesh)]:
    return PVertex(i, Pointer(to=mesh))


@always_inline
def pface(ref mesh: Mesh64, i: Int) -> PFace[origin_of(mesh)]:
    return PFace(i, Pointer(to=mesh))


@always_inline
def phalfedge(ref mesh: Mesh64, i: Int) -> PHalfedge[origin_of(mesh)]:
    return PHalfedge(i, Pointer(to=mesh))


# ===----------------------------------------------------------------------=== #
# Kernels: ring (neighbor position sums)
# ===----------------------------------------------------------------------=== #


def ring_raw[I: DType](mesh: Mesh[I], mut out: List[Vector3]):
    for v in range(mesh.n_vertices):
        var acc = Vector3.zero()
        var start = Int(mesh.v_halfedge[v])
        var h = start
        while True:
            acc += mesh.positions[Int(mesh.he_vertex[h ^ 1])]
            h = Int(mesh.he_next[h ^ 1])
            if h == start:
                break
        out[v] = acc


def ring_index[I: DType](mesh: Mesh[I], mut out: List[Vector3]):
    for i in range(mesh.n_vertices):
        var acc = Vector3.zero()
        var start = nav_halfedge(mesh, Vertex(i))
        var h = start
        while True:
            acc += mesh.positions[nav_vertex(mesh, nav_twin(mesh, h)).index]
            h = nav_next(mesh, nav_twin(mesh, h))
            if h == start:
                break
        out[i] = acc


def ring_index_iter[I: DType](mesh: Mesh[I], mut out: List[Vector3]):
    for i in range(mesh.n_vertices):
        var acc = Vector3.zero()
        for h in outgoing(mesh, Vertex(i)):
            acc += mesh.positions[nav_vertex(mesh, nav_twin(mesh, h)).index]
        out[i] = acc


def ring_pointer(mesh: Mesh64, mut out: List[Vector3]):
    for i in range(mesh.n_vertices):
        var acc = Vector3.zero()
        var start = pvertex(mesh, i).halfedge()
        var h = start
        while True:
            acc += mesh.positions[h.twin().vertex().index]
            h = h.twin().next()
            if h == start:
                break
        out[i] = acc


def ring_pointer_iter(mesh: Mesh64, mut out: List[Vector3]):
    for i in range(mesh.n_vertices):
        var acc = Vector3.zero()
        for h in pvertex(mesh, i).outgoing():
            acc += mesh.positions[h.twin().vertex().index]
        out[i] = acc


# ===----------------------------------------------------------------------=== #
# Kernels: connect (navigation only)
# ===----------------------------------------------------------------------=== #


def connect_raw[I: DType](mesh: Mesh[I]) -> Int:
    var s = 0
    for v in range(mesh.n_vertices):
        var start = Int(mesh.v_halfedge[v])
        var h = start
        while True:
            s += Int(mesh.he_vertex[h ^ 1]) + Int(mesh.he_face[h])
            h = Int(mesh.he_next[h ^ 1])
            if h == start:
                break
    return s


def connect_index[I: DType](mesh: Mesh[I]) -> Int:
    var s = 0
    for i in range(mesh.n_vertices):
        var start = nav_halfedge(mesh, Vertex(i))
        var h = start
        while True:
            s += (
                nav_vertex(mesh, nav_twin(mesh, h)).index
                + nav_face(mesh, h).index
            )
            h = nav_next(mesh, nav_twin(mesh, h))
            if h == start:
                break
    return s


def connect_index_iter[I: DType](mesh: Mesh[I]) -> Int:
    var s = 0
    for i in range(mesh.n_vertices):
        for h in outgoing(mesh, Vertex(i)):
            s += (
                nav_vertex(mesh, nav_twin(mesh, h)).index
                + nav_face(mesh, h).index
            )
    return s


def connect_pointer(mesh: Mesh64) -> Int:
    var s = 0
    for i in range(mesh.n_vertices):
        var start = pvertex(mesh, i).halfedge()
        var h = start
        while True:
            s += h.twin().vertex().index + h.face().index
            h = h.twin().next()
            if h == start:
                break
    return s


def connect_pointer_iter(mesh: Mesh64) -> Int:
    var s = 0
    for i in range(mesh.n_vertices):
        for h in pvertex(mesh, i).outgoing():
            s += h.twin().vertex().index + h.face().index
    return s


# ===----------------------------------------------------------------------=== #
# Kernels: faces (areas) and cotan (edge weights)
# ===----------------------------------------------------------------------=== #


@always_inline
def _area(a: Vector3, b: Vector3, c: Vector3) -> Float64:
    return 0.5 * (b - a).cross(c - a).norm()


@always_inline
def _cot(tail: Vector3, tip: Vector3, corner: Vector3) -> Float64:
    var u = tail - corner
    var v = tip - corner
    return u.dot(v) / u.cross(v).norm()


def faces_raw[I: DType](mesh: Mesh[I], mut out: List[Float64]):
    for f in range(mesh.n_faces):
        var h = Int(mesh.f_halfedge[f])
        var h1 = Int(mesh.he_next[h])
        var h2 = Int(mesh.he_next[h1])
        out[f] = _area(
            mesh.positions[Int(mesh.he_vertex[h])],
            mesh.positions[Int(mesh.he_vertex[h1])],
            mesh.positions[Int(mesh.he_vertex[h2])],
        )


def faces_index[I: DType](mesh: Mesh[I], mut out: List[Float64]):
    for i in range(mesh.n_faces):
        var h = nav_halfedge(mesh, Face(i))
        var h1 = nav_next(mesh, h)
        var h2 = nav_next(mesh, h1)
        out[i] = _area(
            mesh.positions[nav_vertex(mesh, h).index],
            mesh.positions[nav_vertex(mesh, h1).index],
            mesh.positions[nav_vertex(mesh, h2).index],
        )


def faces_pointer(mesh: Mesh64, mut out: List[Float64]):
    for i in range(mesh.n_faces):
        var h = pface(mesh, i).halfedge()
        out[i] = _area(
            mesh.positions[h.vertex().index],
            mesh.positions[h.next().vertex().index],
            mesh.positions[h.next().next().vertex().index],
        )


def cotan_raw[I: DType](mesh: Mesh[I], mut out: List[Float64]):
    for e in range(mesh.n_edges):
        var h = 2 * e
        var t = h + 1
        var a = mesh.positions[Int(mesh.he_vertex[h])]
        var b = mesh.positions[Int(mesh.he_vertex[t])]
        var c = mesh.positions[
            Int(mesh.he_vertex[Int(mesh.he_next[Int(mesh.he_next[h])])])
        ]
        var d = mesh.positions[
            Int(mesh.he_vertex[Int(mesh.he_next[Int(mesh.he_next[t])])])
        ]
        out[e] = 0.5 * (_cot(a, b, c) + _cot(b, a, d))


def cotan_index[I: DType](mesh: Mesh[I], mut out: List[Float64]):
    for e in range(mesh.n_edges):
        var h = Halfedge(2 * e)
        var t = nav_twin(mesh, h)
        var a = mesh.positions[nav_vertex(mesh, h).index]
        var b = mesh.positions[nav_vertex(mesh, t).index]
        var c = mesh.positions[
            nav_vertex(mesh, nav_next(mesh, nav_next(mesh, h))).index
        ]
        var d = mesh.positions[
            nav_vertex(mesh, nav_next(mesh, nav_next(mesh, t))).index
        ]
        out[e] = 0.5 * (_cot(a, b, c) + _cot(b, a, d))


def cotan_pointer(mesh: Mesh64, mut out: List[Float64]):
    for e in range(mesh.n_edges):
        var h = phalfedge(mesh, 2 * e)
        var t = h.twin()
        var a = mesh.positions[h.vertex().index]
        var b = mesh.positions[t.vertex().index]
        var c = mesh.positions[h.next().next().vertex().index]
        var d = mesh.positions[t.next().next().vertex().index]
        out[e] = 0.5 * (_cot(a, b, c) + _cot(b, a, d))


# ===----------------------------------------------------------------------=== #
# Timing
# ===----------------------------------------------------------------------=== #


def _ms(t0: Int) -> Float64:
    return Float64(Int(perf_counter_ns()) - t0) / 1e6


def _sum3(xs: List[Vector3]) -> Float64:
    var s = 0.0
    for v in xs:
        s += v.sum()
    return s


def _sum(xs: List[Float64]) -> Float64:
    var s = 0.0
    for x in xs:
        s += x
    return s


def _row(kernel: String, variant: String, ms: Float64, n: Int, check: String):
    var ns = ms * 1e6 / Float64(n)
    print(
        "  ",
        kernel,
        " " * (8 - kernel.byte_length()),
        variant,
        " " * (13 - variant.byte_length()),
        round(ms, 2),
        "ms  ",
        round(ns, 2),
        "ns/elem   check",
        check,
    )


def bench_ring[
    I: DType, //, f: def(Mesh[I], mut List[Vector3]) thin -> None
](mesh: Mesh[I], name: String):
    var out = List[Vector3](length=mesh.n_vertices, fill=Vector3.zero())
    var best = 1e30
    for _ in range(RUNS):
        var t = Int(perf_counter_ns())
        f(mesh, out)
        best = min(best, _ms(t))
    _row("ring", name, best, mesh.n_vertices, String(_sum3(out)))


def bench_connect[
    I: DType, //, f: def(Mesh[I]) thin -> Int
](mesh: Mesh[I], name: String):
    var best = 1e30
    var s = 0
    for _ in range(RUNS):
        var t = Int(perf_counter_ns())
        s = f(mesh)
        best = min(best, _ms(t))
    _row("connect", name, best, mesh.n_vertices, String(s))


def bench_faces[
    I: DType, //, f: def(Mesh[I], mut List[Float64]) thin -> None
](mesh: Mesh[I], name: String, kernel: String, n: Int):
    var out = List[Float64](length=n, fill=0.0)
    var best = 1e30
    for _ in range(RUNS):
        var t = Int(perf_counter_ns())
        f(mesh, out)
        best = min(best, _ms(t))
    _row(kernel, name, best, n, String(_sum(out)))


def run(mesh: Mesh64, m32: Mesh32):
    comptime I64 = DType.int64
    comptime I32 = DType.int32
    bench_ring[ring_raw[I64]](mesh, "raw")
    bench_ring[ring_raw[I32]](m32, "raw32")
    bench_ring[ring_index[I64]](mesh, "index")
    bench_ring[ring_index[I32]](m32, "index32")
    bench_ring[ring_index_iter[I64]](mesh, "index-iter")
    bench_ring[ring_index_iter[I32]](m32, "index-iter32")
    bench_ring[ring_pointer](mesh, "pointer")
    bench_ring[ring_pointer_iter](mesh, "pointer-iter")
    bench_connect[connect_raw[I64]](mesh, "raw")
    bench_connect[connect_raw[I32]](m32, "raw32")
    bench_connect[connect_index[I64]](mesh, "index")
    bench_connect[connect_index[I32]](m32, "index32")
    bench_connect[connect_index_iter[I64]](mesh, "index-iter")
    bench_connect[connect_index_iter[I32]](m32, "index-iter32")
    bench_connect[connect_pointer](mesh, "pointer")
    bench_connect[connect_pointer_iter](mesh, "pointer-iter")
    var nf = mesh.n_faces
    var ne = mesh.n_edges
    bench_faces[faces_raw[I64]](mesh, "raw", "faces", nf)
    bench_faces[faces_raw[I32]](m32, "raw32", "faces", nf)
    bench_faces[faces_index[I64]](mesh, "index", "faces", nf)
    bench_faces[faces_index[I32]](m32, "index32", "faces", nf)
    bench_faces[faces_pointer](mesh, "pointer", "faces", nf)
    bench_faces[cotan_raw[I64]](mesh, "raw", "cotan", ne)
    bench_faces[cotan_raw[I32]](m32, "raw32", "cotan", ne)
    bench_faces[cotan_index[I64]](mesh, "index", "cotan", ne)
    bench_faces[cotan_index[I32]](m32, "index32", "cotan", ne)
    bench_faces[cotan_pointer](mesh, "pointer", "cotan", ne)


def main() raises:
    seed(20261118)
    print(
        "handle sizes in bytes: index",
        size_of[Halfedge](),
        " pointer",
        size_of[PHalfedge[ImmutAnyOrigin]](),
    )
    print("best of", RUNS, "runs; ns per vertex (ring, connect), face, or edge")
    var sizes: List[Int] = [100, 708]
    for n in sizes:
        for shuffle in [False, True]:
            var mesh = torus(n, n, shuffle)
            var m32 = Mesh32(narrowing=mesh)
            print(
                "\ntorus",
                n,
                "x",
                n,
                ":",
                mesh.n_vertices,
                "vertices,",
                mesh.n_faces,
                "faces,",
                "shuffled" if shuffle else "grid order",
            )
            run(mesh, m32)
