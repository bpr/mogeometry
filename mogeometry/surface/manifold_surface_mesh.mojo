# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's include/geometrycentral/surface/
# surface_mesh.h and .ipp, manifold_surface_mesh.h, halfedge_element_types.ipp
# (the navigation), src/surface/surface_mesh.cpp, and
# src/surface/manifold_surface_mesh.cpp (commit 019669d, MIT License,
# Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
# contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Manifold surface meshes (`mogeometry.surface.manifold_surface_mesh`)
====================================================================
`ManifoldSurfaceMesh`, geometry-central's halfedge mesh of an oriented
manifold surface, possibly with boundary, with polygon faces: built from
a list of faces, it stores for each halfedge its next halfedge, its tail
vertex, and its face, and for each vertex and face one halfedge. Twins
are implicit (`twin(h)` is `h ^ 1`, `edge(h)` is `h / 2`), and each hole
is a boundary loop, stored as an extra face after the real ones.

Navigation is on the mesh, with handles (`Vertex`, `Halfedge`, ...) as
arguments:

```mojo
var mesh = ManifoldSurfaceMesh(faces)
for v in mesh.vertices():
    for h in mesh.outgoing_halfedges(v):
        var w = mesh.tip_vertex(h)
```

The element order (indices) is geometry-central's exactly, so results
can be compared with it element by element.

Differences from geometry-central:

- Navigation and element queries are mesh methods: `he.twin().vertex()`
  is `mesh.vertex(mesh.twin(he))`, `v.degree()` is `mesh.degree(v)`,
  `f.adjacentVertices()` is `mesh.adjacent_vertices(f)`.
- Only the manifold mesh: there is no general `SurfaceMesh` (no
  nonmanifold or unoriented meshes, no `sibling` lists), so queries
  geometry-central answers in general (`isManifold`, `isOriented`) are
  constant here.
- Connectivity is stored as `MESH_INDEX` (`Int32`), so a mesh has fewer
  than 2^31 halfedges.
- Invalid input raises (geometry-central throws from `GC_SAFETY_ASSERT`,
  or skips the check when built with `NGC_SAFETY_CHECKS`).
- Not yet ported: the mutations (`flip`, `splitEdgeTriangular`,
  `collapseEdgeTriangular`, ...), `compress()`, the construction from a
  twin list, `getVertexManifoldStatus` and the other constant status
  queries, and conversion to a general mesh.
"""

from mogeometry.surface.halfedge_element_types import (
    BoundaryLoop,
    BoundaryLoopData,
    Corner,
    CornerData,
    DiamondBoundary,
    Edge,
    EdgeAdjacentFaceNavigator,
    EdgeAdjacentHalfedgeNavigator,
    EdgeAdjacentInteriorHalfedgeNavigator,
    EdgeData,
    EdgeVertices,
    ElementRange,
    Face,
    FaceAdjacentCornerNavigator,
    FaceAdjacentEdgeNavigator,
    FaceAdjacentFaceNavigator,
    FaceAdjacentHalfedgeNavigator,
    FaceAdjacentVertexNavigator,
    FaceData,
    Halfedge,
    HalfedgeData,
    NavigationIterator,
    RangeFilter,
    Vertex,
    VertexAdjacentCornerNavigator,
    VertexAdjacentEdgeNavigator,
    VertexAdjacentFaceNavigator,
    VertexAdjacentVertexNavigator,
    VertexData,
    VertexIncomingHalfedgeNavigator,
    VertexOutgoingHalfedgeNavigator,
)
from mogeometry.utilities.disjoint_sets import DisjointSets
from mogeometry.utilities.element import (
    INVALID_IND,
    MESH_CHECKS,
    MESH_INDEX,
    ElementHandle,
    ElementKind,
    ElementMesh,
    MeshTag,
)

comptime _Index = Scalar[MESH_INDEX]
"""A stored connectivity entry."""

comptime _MAX_INDEX = Int(_Index.MAX)
"""The largest index the storage holds."""


def _fnv1a(hash: UInt64, x: Int) -> UInt64:
    """One step of 64-bit FNV-1a, over the eight bytes of `x`."""
    var h = hash
    var u = UInt64(x)
    for k in range(8):
        h = (h ^ ((u >> UInt64(8 * k)) & 0xFF)) * 0x100000001B3
    return h


def _c_div2(x: Int) -> Int:
    """`x / 2` rounded toward zero, as C++ divides `int`s."""
    return x // 2 if x >= 0 else -((-x) // 2)


struct ManifoldSurfaceMesh(Copyable, ElementMesh, Movable):
    """A halfedge mesh of an oriented manifold surface, possibly with
    boundary, with polygon faces (geometry-central's
    `ManifoldSurfaceMesh`).

    Copies (`mesh.copy()`) have the same elements, and their handles are
    interchangeable.
    """

    var _he_next_arr: List[_Index]
    """Each halfedge's next halfedge in its face (or boundary loop)."""
    var _he_vertex_arr: List[_Index]
    """Each halfedge's tail vertex."""
    var _he_face_arr: List[_Index]
    """Each halfedge's face; a boundary loop's halfedges hold its face
    index (`n_faces_capacity - 1 - loop index`)."""
    var _v_halfedge_arr: List[_Index]
    """An outgoing halfedge of each vertex: the interior one along the
    boundary, for a boundary vertex."""
    var _f_halfedge_arr: List[_Index]
    """A halfedge of each face, then of each boundary loop (from the back:
    loop 0 is the last entry)."""

    var _n_halfedges_count: Int
    """The number of live halfedges."""
    var _n_interior_halfedges_count: Int
    """The number of live halfedges in faces."""
    var _n_edges_count: Int
    """The number of live edges."""
    var _n_vertices_count: Int
    """The number of live vertices."""
    var _n_faces_count: Int
    """The number of live faces."""
    var _n_boundary_loops_count: Int
    """The number of live boundary loops."""

    var _n_vertices_capacity_count: Int
    """The size of the vertex arrays."""
    var _n_halfedges_capacity_count: Int
    """The size of the halfedge arrays."""
    var _n_edges_capacity_count: Int
    """Half the halfedge capacity."""
    var _n_faces_capacity_count: Int
    """The size of the face array (faces and boundary loops)."""

    var _n_vertices_fill_count: Int
    """One past the largest vertex index in use."""
    var _n_halfedges_fill_count: Int
    """One past the largest halfedge index in use."""
    var _n_edges_fill_count: Int
    """One past the largest edge index in use."""
    var _n_faces_fill_count: Int
    """Where the real faces end in the face array."""
    var _n_boundary_loops_fill_count: Int
    """How many boundary loops fill the face array from the back."""

    var _is_compressed: Bool
    """Whether every index space is dense."""
    var _modification_tick: Int
    """Counts changes to the mesh."""
    var _id: Int
    """The mesh's id (`MeshTag`): a hash of its input."""
    var _generation: Int
    """How many times the mesh has been compressed."""
    var _tag: Int
    """`MeshTag.make(_id, _generation)`, which handles carry in debug
    builds."""

    # ===------------------------------------------------------------------=== #
    # Construction
    # ===------------------------------------------------------------------=== #

    def __init__(out self, polygons: List[List[Int]]) raises:
        """A mesh from faces, each a list of 0-based vertex indices in
        counterclockwise order.

        Every index from 0 to the largest must appear in some face. Vertex
        and face indices are kept; edge `e` is the `e`-th distinct edge met
        reading the faces in order, and its halfedge `2e` points the way
        the first face to use it goes.

        Args:
            polygons: The faces.

        Raises:
            If a face has fewer than three vertices, an index is negative
            or unused, a face repeats a vertex consecutively, a directed
            edge appears twice (three faces at an edge, or two faces
            oriented inconsistently), a vertex is on two boundary loops or
            has disconnected neighborhoods (an hourglass), or the mesh has
            too many halfedges for `MESH_INDEX`.
        """
        # == Count, check, and hash the input
        var n_faces = len(polygons)
        var n_vertices = 0
        var n_slots = 0
        var hash: UInt64 = 0xCBF29CE484222325
        for f in range(n_faces):
            ref poly = polygons[f]
            if len(poly) < 3:
                raise Error(
                    String(
                        "ManifoldSurfaceMesh: face ",
                        f,
                        " has degree ",
                        len(poly),
                        "; faces must have degree >= 3",
                    )
                )
            hash = _fnv1a(hash, len(poly))
            for v in poly:
                if v < 0:
                    raise Error(
                        String("ManifoldSurfaceMesh: negative vertex index ", v)
                    )
                n_vertices = max(n_vertices, v)
                hash = _fnv1a(hash, v)
            n_slots += len(poly)
        n_vertices += 1  # 0-based, so the count is the largest index + 1
        if (
            n_vertices > _MAX_INDEX
            or n_slots > _MAX_INDEX
            or n_faces > _MAX_INDEX
        ):
            raise Error("ManifoldSurfaceMesh: too large for MESH_INDEX")

        # == Each face-halfedge slot: its face's first slot, tail, and tip
        var face_start = List[Int](capacity=n_faces + 1)
        var tail = List[Int](capacity=n_slots)
        var tip = List[Int](capacity=n_slots)
        for f in range(n_faces):
            ref poly = polygons[f]
            face_start.append(len(tail))
            var d = len(poly)
            for k in range(d):
                var a = poly[k]
                var b = poly[(k + 1) % d]
                if a == b:
                    raise Error(
                        String(
                            "ManifoldSurfaceMesh: self-edge in face list ",
                            a,
                            " -- ",
                            b,
                        )
                    )
                tail.append(a)
                tip.append(b)
        face_start.append(n_slots)

        # == Slots sorted by (tail, tip), by two stable counting sorts, so a
        # directed edge's slot is found by binary search among its tail's
        # (geometry-central uses a hash map; the result is the same)
        var count = List[Int](length=n_vertices + 1, fill=0)
        for s in range(n_slots):
            count[tip[s] + 1] += 1
        for v in range(n_vertices):
            count[v + 1] += count[v]
        var by_tip = List[Int](length=n_slots, fill=0)
        for s in range(n_slots):
            by_tip[count[tip[s]]] = s
            count[tip[s]] += 1
        var tail_start = List[Int](length=n_vertices + 1, fill=0)
        for s in range(n_slots):
            tail_start[tail[s] + 1] += 1
        for v in range(n_vertices):
            tail_start[v + 1] += tail_start[v]
        var fill = tail_start.copy()
        var order = List[Int](length=n_slots, fill=0)
        for k in range(n_slots):
            var s = by_tip[k]
            order[fill[tail[s]]] = s
            fill[tail[s]] += 1
        for k in range(1, n_slots):
            var s = order[k]
            var r = order[k - 1]
            if tail[s] == tail[r] and tip[s] == tip[r]:
                raise Error(
                    String(
                        "ManifoldSurfaceMesh: duplicate edge in list ",
                        tail[s],
                        " -- ",
                        tip[s],
                    )
                )

        # == Walk the faces, creating halfedges and hooking up pointers
        var he_next = List[_Index]()
        var he_vertex = List[_Index]()
        var he_face = List[_Index]()
        var v_halfedge = List[_Index](length=n_vertices, fill=-1)
        var f_halfedge = List[_Index](length=n_faces, fill=-1)
        var he_of_slot = List[Int](length=n_slots, fill=INVALID_IND)
        var vertex_used = List[Bool](length=n_vertices, fill=False)
        for f in range(n_faces):
            var prev = INVALID_IND
            var first = INVALID_IND
            for s in range(face_start[f], face_start[f + 1]):
                var a = tail[s]
                var b = tip[s]
                vertex_used[a] = True
                # the twin's slot, if that directed edge appears
                var t = INVALID_IND
                var lo = tail_start[b]
                var hi = tail_start[b + 1]
                while lo < hi:
                    var mid = (lo + hi) // 2
                    if tip[order[mid]] < a:
                        lo = mid + 1
                    else:
                        hi = mid
                if lo < tail_start[b + 1] and tip[order[lo]] == a:
                    t = order[lo]
                var he: Int
                if t != INVALID_IND and he_of_slot[t] != INVALID_IND:
                    # the twin was created first
                    he = he_of_slot[t] ^ 1
                else:
                    # a new edge: halfedge 2e here, its twin 2e + 1
                    he = len(he_next)
                    if he + 1 > _MAX_INDEX:
                        raise Error(
                            "ManifoldSurfaceMesh: too large for MESH_INDEX"
                        )
                    he_next.append(-1)
                    he_next.append(-1)
                    he_vertex.append(_Index(a))
                    he_vertex.append(_Index(b))
                    he_face.append(-1)
                    he_face.append(-1)
                he_of_slot[s] = he
                he_face[he] = _Index(f)
                v_halfedge[a] = _Index(he)
                if s == face_start[f]:
                    f_halfedge[f] = _Index(he)
                    first = he
                else:
                    he_next[prev] = _Index(he)
                prev = he
            he_next[prev] = _Index(first)

        for v in range(n_vertices):
            if not vertex_used[v]:
                raise Error(
                    String("ManifoldSurfaceMesh: unreferenced vertex ", v)
                )

        # == Each boundary neighborhood must be a disk or a half-disk
        var n_halfedges = len(he_next)
        var on_boundary = List[Bool](length=n_vertices, fill=False)
        for h in range(n_halfedges):
            if Int(he_next[h]) == INVALID_IND:
                var v = Int(he_vertex[h])
                if on_boundary[v]:
                    raise Error(
                        String(
                            "ManifoldSurfaceMesh: vertex ",
                            v,
                            " appears in more than one boundary loop",
                        )
                    )
                on_boundary[v] = True

        # == Resolve boundary loops: the halfedges with no face yet
        var n_interior = n_halfedges
        var n_loops = 0
        for i in range(n_halfedges):
            if Int(he_face[i]) != INVALID_IND:
                continue
            var loop_face = n_faces + n_loops
            f_halfedge.append(_Index(i))
            n_loops += 1
            # walk around the loop (clockwise)
            var cur = i
            var prev: Int
            var loop_count = 0
            while True:
                he_face[cur] = _Index(loop_face)
                # the twin is interior and on the boundary: make it its
                # vertex's halfedge
                var cur_t = cur ^ 1
                v_halfedge[Int(he_vertex[cur_t])] = _Index(cur_t)
                n_interior -= 1
                prev = cur
                cur = Int(he_next[cur ^ 1]) ^ 1
                var inner = 0
                while Int(he_face[cur]) != INVALID_IND:
                    if cur == i:
                        break
                    cur = Int(he_next[cur]) ^ 1
                    inner += 1
                    if inner >= n_halfedges:
                        raise Error(
                            "ManifoldSurfaceMesh: boundary infinite loop orbit"
                        )
                he_next[cur] = _Index(prev)
                loop_count += 1
                if loop_count >= n_halfedges:
                    raise Error("ManifoldSurfaceMesh: boundary infinite loop")
                if cur == i:
                    break

        # == Each vertex's outgoing halfedges must form one orbit
        var seen = List[Bool](length=n_halfedges, fill=False)
        for v in range(n_vertices):
            var first = Int(v_halfedge[v])
            var cur = first
            while True:
                if seen[cur]:
                    raise Error(
                        "ManifoldSurfaceMesh: somehow encountered outgoing"
                        " halfedge before orbiting v"
                    )
                seen[cur] = True
                cur = Int(he_next[cur ^ 1])
                if cur == first:
                    break
        for h in range(n_halfedges):
            if not seen[h]:
                raise Error(
                    String(
                        "ManifoldSurfaceMesh: mesh not manifold. Vertex ",
                        Int(he_vertex[h]),
                        " has disconnected neighborhoods incident (imagine an",
                        " hourglass)",
                    )
                )

        self._he_next_arr = he_next^
        self._he_vertex_arr = he_vertex^
        self._he_face_arr = he_face^
        self._v_halfedge_arr = v_halfedge^
        self._f_halfedge_arr = f_halfedge^
        self._n_halfedges_count = n_halfedges
        self._n_interior_halfedges_count = n_interior
        self._n_edges_count = n_halfedges // 2
        self._n_vertices_count = n_vertices
        self._n_faces_count = n_faces
        self._n_boundary_loops_count = n_loops
        self._n_vertices_capacity_count = n_vertices
        self._n_halfedges_capacity_count = n_halfedges
        self._n_edges_capacity_count = n_halfedges // 2
        self._n_faces_capacity_count = n_faces + n_loops
        self._n_vertices_fill_count = n_vertices
        self._n_halfedges_fill_count = n_halfedges
        self._n_edges_fill_count = n_halfedges // 2
        self._n_faces_fill_count = n_faces
        self._n_boundary_loops_fill_count = n_loops
        self._is_compressed = True
        self._modification_tick = 1
        self._id = MeshTag.id_from_hash(hash)
        self._generation = 0
        self._tag = MeshTag.make(self._id, self._generation)

    # ===------------------------------------------------------------------=== #
    # Element counts, capacities, and index sizes
    # ===------------------------------------------------------------------=== #

    def n_halfedges(self) -> Int:
        """The number of halfedges, interior and exterior.

        Returns:
            The count.
        """
        return self._n_halfedges_count

    def n_interior_halfedges(self) -> Int:
        """The number of halfedges in faces.

        Returns:
            The count.
        """
        return self._n_interior_halfedges_count

    def n_exterior_halfedges(self) -> Int:
        """The number of halfedges in boundary loops.

        Returns:
            The count.
        """
        return self._n_halfedges_count - self._n_interior_halfedges_count

    def n_corners(self) -> Int:
        """The number of corners (interior halfedges).

        Returns:
            The count.
        """
        return self._n_interior_halfedges_count

    def n_vertices(self) -> Int:
        """The number of vertices.

        Returns:
            The count.
        """
        return self._n_vertices_count

    def n_interior_vertices(self) -> Int:
        """The number of vertices not on the boundary (O(n)).

        Returns:
            The count.
        """
        var n = 0
        for v in self.vertices():
            if not self.is_boundary(v):
                n += 1
        return n

    def n_edges(self) -> Int:
        """The number of edges.

        Returns:
            The count.
        """
        return self._n_edges_count

    def n_faces(self) -> Int:
        """The number of faces (not counting boundary loops).

        Returns:
            The count.
        """
        return self._n_faces_count

    def n_boundary_loops(self) -> Int:
        """The number of boundary loops (holes).

        Returns:
            The count.
        """
        return self._n_boundary_loops_count

    def n_halfedges_capacity(self) -> Int:
        """The number of halfedges the arrays have room for.

        Returns:
            The capacity.
        """
        return self._n_halfedges_capacity_count

    def n_corners_capacity(self) -> Int:
        """The number of corners the arrays have room for.

        Returns:
            The capacity (the halfedges').
        """
        return self._n_halfedges_capacity_count

    def n_vertices_capacity(self) -> Int:
        """The number of vertices the arrays have room for.

        Returns:
            The capacity.
        """
        return self._n_vertices_capacity_count

    def n_edges_capacity(self) -> Int:
        """The number of edges the arrays have room for.

        Returns:
            The capacity.
        """
        return self._n_edges_capacity_count

    def n_faces_capacity(self) -> Int:
        """The number of faces the arrays have room for (boundary loops
        take the rest).

        Returns:
            The capacity.
        """
        return self._n_faces_capacity_count - self._n_boundary_loops_fill_count

    def n_boundary_loops_capacity(self) -> Int:
        """The number of boundary loops the arrays have room for.

        Returns:
            The capacity.
        """
        return self._n_faces_capacity_count - self._n_faces_fill_count

    def halfedge_index_size(self) -> Int:
        """One past the largest halfedge index in use.

        Returns:
            The index size.
        """
        return self._n_halfedges_fill_count

    def corner_index_size(self) -> Int:
        """At least one past the largest corner index in use (it may be one
        more, as in geometry-central).

        Returns:
            The index size.
        """
        return self._n_halfedges_fill_count

    def vertex_index_size(self) -> Int:
        """One past the largest vertex index in use.

        Returns:
            The index size.
        """
        return self._n_vertices_fill_count

    def edge_index_size(self) -> Int:
        """One past the largest edge index in use.

        Returns:
            The index size.
        """
        return self._n_edges_fill_count

    def face_index_size(self) -> Int:
        """One past the largest face index in use (boundary loops' face
        indices come after).

        Returns:
            The index size.
        """
        return self._n_faces_fill_count

    def boundary_loop_index_size(self) -> Int:
        """One past the largest boundary loop index in use.

        Returns:
            The index size.
        """
        return self._n_boundary_loops_fill_count

    def is_compressed(self) -> Bool:
        """Whether every index space is dense (no deleted elements).

        Returns:
            True for a mesh that has not lost elements since it was built
            or compressed.
        """
        return self._is_compressed

    # ===------------------------------------------------------------------=== #
    # ElementMesh, tags, and checks
    # ===------------------------------------------------------------------=== #

    def element_capacity(self, kind: Int) -> Int:
        """The number of elements of a kind a container must hold.

        Args:
            kind: An `ElementKind` value.

        Returns:
            The capacity for that kind.
        """
        if kind == ElementKind.VERTEX:
            return self.n_vertices_capacity()
        if kind == ElementKind.HALFEDGE or kind == ElementKind.CORNER:
            return self.n_halfedges_capacity()
        if kind == ElementKind.EDGE:
            return self.n_edges_capacity()
        if kind == ElementKind.FACE:
            return self.n_faces_capacity()
        return self.n_boundary_loops_capacity()

    @always_inline
    def mesh_tag(self) -> Int:
        """The mesh's tag (`MeshTag`), which its handles carry in debug
        builds.

        Returns:
            The tag.
        """
        return self._tag

    @always_inline
    def _index_size[kind: Int](self) -> Int:
        """One past the largest valid index of a kind of handle."""
        comptime if kind == ElementKind.VERTEX:
            return self._n_vertices_fill_count
        elif kind == ElementKind.HALFEDGE or kind == ElementKind.CORNER:
            return self._n_halfedges_fill_count
        elif kind == ElementKind.EDGE:
            return self._n_edges_fill_count
        elif kind == ElementKind.FACE:
            return self._n_faces_capacity_count
        else:
            return self._n_boundary_loops_fill_count

    def owns[kind: Int](self, e: ElementHandle[kind]) -> Bool:
        """Whether a handle can be used with this mesh: its index is in
        range and, in debug builds, it carries this mesh's tag or none.

        Mesh methods assert this of every handle they take (in debug
        builds); this is the same check, without stopping the program.

        Parameters:
            kind: The `ElementKind` of the handle.

        Args:
            e: The handle.

        Returns:
            True if the handle is valid here.
        """
        var ok = 0 <= e.index() < self._index_size[kind]()
        comptime if MESH_CHECKS:
            ok = ok and (e.tag() == 0 or e.tag() == self._tag)
        return ok

    @always_inline
    def _check[kind: Int](self, e: ElementHandle[kind]) -> Int:
        """The handle's index, asserting (in debug builds) that it belongs
        to this mesh."""
        comptime if MESH_CHECKS:
            debug_assert[assert_mode="safe"](
                e.tag() == 0 or e.tag() == self._tag,
                "ManifoldSurfaceMesh: element handle from another mesh, or",
                " from before compress()",
            )
            debug_assert[assert_mode="safe"](
                0 <= e.index() < self._index_size[kind](),
                "ManifoldSurfaceMesh: element index out of range",
            )
        return e.index()

    @always_inline
    def _halfedge_at(self, i: Int) -> Halfedge:
        return Halfedge(index=i, tag=self._tag)

    @always_inline
    def _vertex_at(self, i: Int) -> Vertex:
        return Vertex(index=i, tag=self._tag)

    @always_inline
    def _face_at(self, i: Int) -> Face:
        return Face(index=i, tag=self._tag)

    @always_inline
    def _face_is_boundary_loop(self, f: Int) -> Bool:
        return f >= self._n_faces_fill_count

    @always_inline
    def _element_okay[kind: Int, filter: Int](self, i: Int) -> Bool:
        """Whether index `i` is a live element an `ElementRange` yields."""
        comptime if kind == ElementKind.VERTEX:
            return Int(self._v_halfedge_arr[i]) != INVALID_IND
        elif kind == ElementKind.HALFEDGE:
            if Int(self._he_next_arr[i]) == INVALID_IND:
                return False
            comptime if filter == RangeFilter.INTERIOR:
                return not self._face_is_boundary_loop(
                    Int(self._he_face_arr[i])
                )
            elif filter == RangeFilter.EXTERIOR:
                return self._face_is_boundary_loop(Int(self._he_face_arr[i]))
            else:
                return True
        elif kind == ElementKind.CORNER:
            return Int(
                self._he_next_arr[i]
            ) != INVALID_IND and not self._face_is_boundary_loop(
                Int(self._he_face_arr[i])
            )
        elif kind == ElementKind.EDGE:
            return Int(self._he_next_arr[2 * i]) != INVALID_IND
        elif kind == ElementKind.FACE:
            return Int(self._f_halfedge_arr[i]) != INVALID_IND
        else:
            return (
                Int(self._f_halfedge_arr[self._n_faces_capacity_count - 1 - i])
                != INVALID_IND
            )

    # ===------------------------------------------------------------------=== #
    # Elements by index, and all elements of a kind
    # ===------------------------------------------------------------------=== #

    def vertex(self, index: Int) -> Vertex:
        """The vertex with an index (meaningful when the mesh is
        compressed).

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return self._vertex_at(self._check(Vertex(index)))

    def halfedge(self, index: Int) -> Halfedge:
        """The halfedge with an index.

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return self._halfedge_at(self._check(Halfedge(index)))

    def corner(self, index: Int) -> Corner:
        """The corner with an index (its halfedge's).

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return Corner(index=self._check(Corner(index)), tag=self._tag)

    def edge(self, index: Int) -> Edge:
        """The edge with an index.

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return Edge(index=self._check(Edge(index)), tag=self._tag)

    def face(self, index: Int) -> Face:
        """The face with an index (boundary loops' face indices come after
        the faces').

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return self._face_at(self._check(Face(index)))

    def boundary_loop(self, index: Int) -> BoundaryLoop:
        """The boundary loop with an index.

        Args:
            index: The index.

        Returns:
            Its handle.
        """
        return BoundaryLoop(
            index=self._check(BoundaryLoop(index)), tag=self._tag
        )

    def vertices(
        ref self,
    ) -> ElementRange[ElementKind.VERTEX, RangeFilter.ALL, origin_of(self)]:
        """Every vertex, by index.

        Returns:
            An iterator over the vertices.
        """
        return ElementRange[ElementKind.VERTEX, RangeFilter.ALL](
            Pointer(to=self), self._n_vertices_fill_count
        )

    def halfedges(
        ref self,
    ) -> ElementRange[ElementKind.HALFEDGE, RangeFilter.ALL, origin_of(self)]:
        """Every halfedge, by index.

        Returns:
            An iterator over the halfedges.
        """
        return ElementRange[ElementKind.HALFEDGE, RangeFilter.ALL](
            Pointer(to=self), self._n_halfedges_fill_count
        )

    def interior_halfedges(
        ref self,
    ) -> ElementRange[
        ElementKind.HALFEDGE, RangeFilter.INTERIOR, origin_of(self)
    ]:
        """Every halfedge in a face, by index.

        Returns:
            An iterator over the interior halfedges.
        """
        return ElementRange[ElementKind.HALFEDGE, RangeFilter.INTERIOR](
            Pointer(to=self), self._n_halfedges_fill_count
        )

    def exterior_halfedges(
        ref self,
    ) -> ElementRange[
        ElementKind.HALFEDGE, RangeFilter.EXTERIOR, origin_of(self)
    ]:
        """Every halfedge in a boundary loop, by index.

        Returns:
            An iterator over the exterior halfedges.
        """
        return ElementRange[ElementKind.HALFEDGE, RangeFilter.EXTERIOR](
            Pointer(to=self), self._n_halfedges_fill_count
        )

    def corners(
        ref self,
    ) -> ElementRange[ElementKind.CORNER, RangeFilter.ALL, origin_of(self)]:
        """Every corner, by index.

        Returns:
            An iterator over the corners.
        """
        return ElementRange[ElementKind.CORNER, RangeFilter.ALL](
            Pointer(to=self), self._n_halfedges_fill_count
        )

    def edges(
        ref self,
    ) -> ElementRange[ElementKind.EDGE, RangeFilter.ALL, origin_of(self)]:
        """Every edge, by index.

        Returns:
            An iterator over the edges.
        """
        return ElementRange[ElementKind.EDGE, RangeFilter.ALL](
            Pointer(to=self), self._n_edges_fill_count
        )

    def faces(
        ref self,
    ) -> ElementRange[ElementKind.FACE, RangeFilter.ALL, origin_of(self)]:
        """Every face (not boundary loops), by index.

        Returns:
            An iterator over the faces.
        """
        return ElementRange[ElementKind.FACE, RangeFilter.ALL](
            Pointer(to=self), self._n_faces_fill_count
        )

    def boundary_loops(
        ref self,
    ) -> ElementRange[
        ElementKind.BOUNDARY_LOOP, RangeFilter.ALL, origin_of(self)
    ]:
        """Every boundary loop, by index.

        Returns:
            An iterator over the boundary loops.
        """
        return ElementRange[ElementKind.BOUNDARY_LOOP, RangeFilter.ALL](
            Pointer(to=self), self._n_boundary_loops_fill_count
        )

    # ===------------------------------------------------------------------=== #
    # Halfedges
    # ===------------------------------------------------------------------=== #

    @always_inline
    def twin(self, h: Halfedge) -> Halfedge:
        """The other halfedge of `h`'s edge.

        Args:
            h: A halfedge.

        Returns:
            Its twin.
        """
        return self._halfedge_at(self._check(h) ^ 1)

    @always_inline
    def sibling(self, h: Halfedge) -> Halfedge:
        """The next halfedge around `h`'s edge: its twin, in a manifold
        mesh.

        Args:
            h: A halfedge.

        Returns:
            Its twin.
        """
        return self.twin(h)

    @always_inline
    def next(self, h: Halfedge) -> Halfedge:
        """The next halfedge around `h`'s face (or boundary loop).

        Args:
            h: A halfedge.

        Returns:
            The next halfedge.
        """
        return self._halfedge_at(Int(self._he_next_arr[self._check(h)]))

    @always_inline
    def vertex(self, h: Halfedge) -> Vertex:
        """The vertex `h` leaves (its tail).

        Args:
            h: A halfedge.

        Returns:
            Its tail vertex.
        """
        return self._vertex_at(Int(self._he_vertex_arr[self._check(h)]))

    @always_inline
    def tail_vertex(self, h: Halfedge) -> Vertex:
        """The vertex `h` leaves.

        Args:
            h: A halfedge.

        Returns:
            Its tail vertex.
        """
        return self.vertex(h)

    @always_inline
    def tip_vertex(self, h: Halfedge) -> Vertex:
        """The vertex `h` points to.

        Args:
            h: A halfedge.

        Returns:
            Its tip vertex.
        """
        return self.vertex(self.next(h))

    @always_inline
    def next_outgoing_neighbor(self, h: Halfedge) -> Halfedge:
        """The next halfedge leaving `h`'s tail (`next(twin(h))`).

        Args:
            h: A halfedge.

        Returns:
            The next outgoing halfedge.
        """
        return self.next(self.twin(h))

    @always_inline
    def next_incoming_neighbor(self, h: Halfedge) -> Halfedge:
        """The next halfedge arriving at `h`'s tip (`twin(next(h))`).

        Args:
            h: A halfedge.

        Returns:
            The next incoming halfedge.
        """
        return self.twin(self.next(h))

    @always_inline
    def edge(self, h: Halfedge) -> Edge:
        """The edge `h` is a side of.

        Args:
            h: A halfedge.

        Returns:
            Its edge.
        """
        return Edge(index=self._check(h) // 2, tag=self._tag)

    @always_inline
    def face(self, h: Halfedge) -> Face:
        """The face `h` is in: a boundary loop's face (see
        `is_boundary_loop`) if `h` is exterior.

        Args:
            h: A halfedge.

        Returns:
            Its face.
        """
        return self._face_at(Int(self._he_face_arr[self._check(h)]))

    @always_inline
    def corner(self, h: Halfedge) -> Corner:
        """The corner at `h`'s tail in its face (meaningful for interior
        halfedges).

        Args:
            h: A halfedge.

        Returns:
            Its corner.
        """
        return Corner(index=self._check(h), tag=self._tag)

    def is_dead(self, h: Halfedge) -> Bool:
        """Whether `h` was deleted.

        Args:
            h: A halfedge.

        Returns:
            True if it is no longer part of the mesh.
        """
        return Int(self._he_next_arr[self._check(h)]) == INVALID_IND

    @always_inline
    def is_interior(self, h: Halfedge) -> Bool:
        """Whether `h` is in a face (rather than a boundary loop).

        Args:
            h: A halfedge.

        Returns:
            True if interior.
        """
        return not self._face_is_boundary_loop(
            Int(self._he_face_arr[self._check(h)])
        )

    def orientation(self, h: Halfedge) -> Bool:
        """Whether `h` points the same way as its edge's halfedge.

        Args:
            h: A halfedge.

        Returns:
            True if `h` is its edge's halfedge.
        """
        return self._check(h) % 2 == 0

    def prev_orbit_face(self, h: Halfedge) -> Halfedge:
        """The halfedge before `h` in its face, found by walking around it.

        Args:
            h: A halfedge.

        Returns:
            The halfedge whose `next` is `h`.
        """
        var cur = h
        while True:
            var nxt = self.next(cur)
            if nxt == h:
                break
            cur = nxt
        return cur

    def prev_orbit_vertex(self, h: Halfedge) -> Halfedge:
        """The halfedge before `h` in its face, found by walking around its
        tip.

        Args:
            h: A halfedge.

        Returns:
            The halfedge whose `next` is `h`.
        """
        var cur = self.twin(h)
        while True:
            var nxt = self.next(cur)
            if nxt == h:
                break
            cur = self.twin(nxt)
        return cur

    # ===------------------------------------------------------------------=== #
    # Corners
    # ===------------------------------------------------------------------=== #

    def halfedge(self, c: Corner) -> Halfedge:
        """The halfedge leaving a corner's vertex in its face.

        Args:
            c: A corner.

        Returns:
            Its halfedge.
        """
        return self._halfedge_at(self._check(c))

    def vertex(self, c: Corner) -> Vertex:
        """A corner's vertex.

        Args:
            c: A corner.

        Returns:
            Its vertex.
        """
        return self.vertex(self.halfedge(c))

    def face(self, c: Corner) -> Face:
        """A corner's face.

        Args:
            c: A corner.

        Returns:
            Its face.
        """
        return self.face(self.halfedge(c))

    def is_dead(self, c: Corner) -> Bool:
        """Whether a corner was deleted.

        Args:
            c: A corner.

        Returns:
            True if its halfedge was.
        """
        return self.is_dead(self.halfedge(c))

    # ===------------------------------------------------------------------=== #
    # Vertices
    # ===------------------------------------------------------------------=== #

    @always_inline
    def halfedge(self, v: Vertex) -> Halfedge:
        """A halfedge leaving `v`: for a boundary vertex, the interior one
        along the boundary (whose twin is exterior).

        Args:
            v: A vertex.

        Returns:
            An outgoing halfedge.
        """
        return self._halfedge_at(Int(self._v_halfedge_arr[self._check(v)]))

    def corner(self, v: Vertex) -> Corner:
        """A corner at `v`: `halfedge(v)`'s.

        Args:
            v: A vertex.

        Returns:
            A corner.
        """
        return self.corner(self.halfedge(v))

    def is_dead(self, v: Vertex) -> Bool:
        """Whether `v` was deleted.

        Args:
            v: A vertex.

        Returns:
            True if it is no longer part of the mesh.
        """
        return Int(self._v_halfedge_arr[self._check(v)]) == INVALID_IND

    def is_boundary(self, v: Vertex) -> Bool:
        """Whether `v` is on the boundary.

        Args:
            v: A vertex.

        Returns:
            True if a boundary loop passes through it.
        """
        return not self.is_interior(self.twin(self.halfedge(v)))

    def is_manifold(self, v: Vertex) -> Bool:
        """Whether `v`'s neighborhood is a disk or half-disk: always, in a
        manifold mesh.

        Args:
            v: A vertex.

        Returns:
            True.
        """
        return True

    def degree(self, v: Vertex) -> Int:
        """The number of edges at `v`.

        Args:
            v: A vertex.

        Returns:
            Its degree.
        """
        var k = 0
        for _ in self.adjacent_edges(v):
            k += 1
        return k

    def face_degree(self, v: Vertex) -> Int:
        """The number of faces at `v` (one less than its degree on the
        boundary).

        Args:
            v: A vertex.

        Returns:
            Its face degree.
        """
        var k = 0
        for _ in self.adjacent_faces(v):
            k += 1
        return k

    def outgoing_halfedges(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexOutgoingHalfedgeNavigator, origin_of(self)]:
        """The halfedges leaving `v`, counterclockwise, from
        `halfedge(v)`.

        Args:
            v: A vertex.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self), VertexOutgoingHalfedgeNavigator(self.halfedge(v))
        )

    def incoming_halfedges(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexIncomingHalfedgeNavigator, origin_of(self)]:
        """The halfedges arriving at `v`, from the one before
        `halfedge(v)` in its face.

        Args:
            v: A vertex.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self),
            VertexIncomingHalfedgeNavigator(
                self.prev_orbit_face(self.halfedge(v))
            ),
        )

    def adjacent_vertices(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexAdjacentVertexNavigator, origin_of(self)]:
        """The neighbors of `v`, in the order of its outgoing halfedges.

        Args:
            v: A vertex.

        Returns:
            An iterator over the vertices.
        """
        return NavigationIterator(
            Pointer(to=self), VertexAdjacentVertexNavigator(self.halfedge(v))
        )

    def adjacent_corners(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexAdjacentCornerNavigator, origin_of(self)]:
        """The corners at `v`.

        Args:
            v: A vertex.

        Returns:
            An iterator over the corners.
        """
        return NavigationIterator(
            Pointer(to=self), VertexAdjacentCornerNavigator(self.halfedge(v))
        )

    def adjacent_edges(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexAdjacentEdgeNavigator, origin_of(self)]:
        """The edges at `v`.

        Args:
            v: A vertex.

        Returns:
            An iterator over the edges.
        """
        return NavigationIterator(
            Pointer(to=self), VertexAdjacentEdgeNavigator(self.halfedge(v))
        )

    def adjacent_faces(
        ref self, v: Vertex
    ) -> NavigationIterator[VertexAdjacentFaceNavigator, origin_of(self)]:
        """The faces at `v`.

        Args:
            v: A vertex.

        Returns:
            An iterator over the faces.
        """
        return NavigationIterator(
            Pointer(to=self), VertexAdjacentFaceNavigator(self.halfedge(v))
        )

    # ===------------------------------------------------------------------=== #
    # Edges
    # ===------------------------------------------------------------------=== #

    @always_inline
    def halfedge(self, e: Edge) -> Halfedge:
        """The edge's halfedge, `2e`: interior unless the edge is on the
        boundary with its face on the other side.

        Args:
            e: An edge.

        Returns:
            Its halfedge.
        """
        return self._halfedge_at(2 * self._check(e))

    def other_vertex(self, e: Edge, v: Vertex) -> Vertex:
        """The end of `e` that is not `v`.

        Args:
            e: An edge.
            v: One of its ends.

        Returns:
            The other end (the tail of `halfedge(e)` if `v` is not an end).
        """
        var h = self.halfedge(e)
        if self.tail_vertex(h) == v:
            return self.tip_vertex(h)
        return self.tail_vertex(h)

    def first_vertex(self, e: Edge) -> Vertex:
        """The tail of `halfedge(e)`.

        Args:
            e: An edge.

        Returns:
            Its first vertex.
        """
        return self.tail_vertex(self.halfedge(e))

    def second_vertex(self, e: Edge) -> Vertex:
        """The tip of `halfedge(e)`.

        Args:
            e: An edge.

        Returns:
            Its second vertex.
        """
        return self.tip_vertex(self.halfedge(e))

    def adjacent_vertices(self, e: Edge) -> EdgeVertices:
        """The two ends of `e`.

        Args:
            e: An edge.

        Returns:
            Its first and second vertices.
        """
        var h = self.halfedge(e)
        return EdgeVertices(self.tail_vertex(h), self.tip_vertex(h))

    def diamond_boundary(self, e: Edge) raises -> DiamondBoundary:
        """The four halfedges around the two triangles at `e`.

        Args:
            e: An interior edge between two triangles.

        Returns:
            `next(h)`, `next(next(h))`, `next(t)`, `next(next(t))`, for
            `h = halfedge(e)` and `t = twin(h)`.

        Raises:
            If `e` is on the boundary or a face beside it is not a
            triangle.
        """
        var h = self.halfedge(e)
        var hn = self.next(h)
        var hnn = self.next(hn)
        var ht = self.twin(h)
        var htn = self.next(ht)
        var htnn = self.next(htn)
        if not self.is_interior(ht):
            raise Error("cannot construct diamond_boundary() of boundary edge")
        if self.next(hnn) != h or self.next(htnn) != ht:
            raise Error(
                "cannot construct diamond_boundary() for non-triangular face"
            )
        return DiamondBoundary(hn, hnn, htn, htnn)

    def is_dead(self, e: Edge) -> Bool:
        """Whether `e` was deleted.

        Args:
            e: An edge.

        Returns:
            True if it is no longer part of the mesh.
        """
        return Int(self._he_next_arr[2 * self._check(e)]) == INVALID_IND

    def is_boundary(self, e: Edge) -> Bool:
        """Whether `e` is on the boundary.

        Args:
            e: An edge.

        Returns:
            True if one of its halfedges is in a boundary loop.
        """
        return not self.is_interior(self.twin(self.halfedge(e)))

    def is_manifold(self, e: Edge) -> Bool:
        """Whether at most two faces meet at `e`: always, in a manifold
        mesh.

        Args:
            e: An edge.

        Returns:
            True.
        """
        return True

    def is_oriented(self, e: Edge) -> Bool:
        """Whether the faces at `e` agree on its orientation: always, in a
        manifold mesh.

        Args:
            e: An edge.

        Returns:
            True.
        """
        return True

    def degree(self, e: Edge) -> Int:
        """The number of faces at `e`.

        Args:
            e: An edge.

        Returns:
            2, or 1 on the boundary.
        """
        var k = 0
        for _ in self.adjacent_interior_halfedges(e):
            k += 1
        return k

    def adjacent_halfedges(
        ref self, e: Edge
    ) -> NavigationIterator[EdgeAdjacentHalfedgeNavigator, origin_of(self)]:
        """The two halfedges of `e`, `halfedge(e)` first.

        Args:
            e: An edge.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self), EdgeAdjacentHalfedgeNavigator(self.halfedge(e))
        )

    def adjacent_interior_halfedges(
        ref self, e: Edge
    ) -> NavigationIterator[
        EdgeAdjacentInteriorHalfedgeNavigator, origin_of(self)
    ]:
        """The halfedges of `e` that are in faces.

        Args:
            e: An edge.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self),
            EdgeAdjacentInteriorHalfedgeNavigator(self.halfedge(e)),
        )

    def adjacent_faces(
        ref self, e: Edge
    ) -> NavigationIterator[EdgeAdjacentFaceNavigator, origin_of(self)]:
        """The faces at `e`.

        Args:
            e: An edge.

        Returns:
            An iterator over the faces.
        """
        return NavigationIterator(
            Pointer(to=self), EdgeAdjacentFaceNavigator(self.halfedge(e))
        )

    # ===------------------------------------------------------------------=== #
    # Faces
    # ===------------------------------------------------------------------=== #

    @always_inline
    def halfedge(self, f: Face) -> Halfedge:
        """A halfedge of `f`.

        Args:
            f: A face (or a boundary loop's face).

        Returns:
            One of its halfedges.
        """
        return self._halfedge_at(Int(self._f_halfedge_arr[self._check(f)]))

    def as_boundary_loop(self, f: Face) -> BoundaryLoop:
        """The boundary loop a boundary-loop face stands for.

        Args:
            f: A face for which `is_boundary_loop(f)`.

        Returns:
            The boundary loop.
        """
        debug_assert[assert_mode="safe"](
            self.is_boundary_loop(f),
            "face must be boundary loop to call as_boundary_loop()",
        )
        return BoundaryLoop(
            index=self._n_faces_capacity_count - 1 - self._check(f),
            tag=self._tag,
        )

    def is_dead(self, f: Face) -> Bool:
        """Whether `f` was deleted.

        Args:
            f: A face.

        Returns:
            True if it is no longer part of the mesh.
        """
        return Int(self._f_halfedge_arr[self._check(f)]) == INVALID_IND

    def is_boundary_loop(self, f: Face) -> Bool:
        """Whether `f` is a boundary loop's face rather than a real face.

        Args:
            f: A face.

        Returns:
            True for a boundary loop.
        """
        return self._face_is_boundary_loop(self._check(f))

    def is_triangle(self, f: Face) -> Bool:
        """Whether `f` has three sides.

        Args:
            f: A face.

        Returns:
            True for a triangle.
        """
        var h = self.halfedge(f)
        return h == self.next(self.next(self.next(h)))

    def degree(self, f: Face) -> Int:
        """The number of sides of `f`.

        Args:
            f: A face.

        Returns:
            Its degree.
        """
        var k = 0
        for _ in self.adjacent_halfedges(f):
            k += 1
        return k

    def adjacent_halfedges(
        ref self, f: Face
    ) -> NavigationIterator[FaceAdjacentHalfedgeNavigator, origin_of(self)]:
        """The halfedges around `f`, from `halfedge(f)`.

        Args:
            f: A face.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentHalfedgeNavigator(self.halfedge(f))
        )

    def adjacent_vertices(
        ref self, f: Face
    ) -> NavigationIterator[FaceAdjacentVertexNavigator, origin_of(self)]:
        """The vertices around `f`, counterclockwise.

        Args:
            f: A face.

        Returns:
            An iterator over the vertices.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentVertexNavigator(self.halfedge(f))
        )

    def adjacent_corners(
        ref self, f: Face
    ) -> NavigationIterator[FaceAdjacentCornerNavigator, origin_of(self)]:
        """The corners of `f`.

        Args:
            f: A face.

        Returns:
            An iterator over the corners.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentCornerNavigator(self.halfedge(f))
        )

    def adjacent_edges(
        ref self, f: Face
    ) -> NavigationIterator[FaceAdjacentEdgeNavigator, origin_of(self)]:
        """The edges around `f`.

        Args:
            f: A face.

        Returns:
            An iterator over the edges.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentEdgeNavigator(self.halfedge(f))
        )

    def adjacent_faces(
        ref self, f: Face
    ) -> NavigationIterator[FaceAdjacentFaceNavigator, origin_of(self)]:
        """The faces across `f`'s edges (none across a boundary edge).

        Args:
            f: A face.

        Returns:
            An iterator over the faces.
        """
        var h = self.halfedge(f)
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentFaceNavigator(h, h)
        )

    # ===------------------------------------------------------------------=== #
    # Boundary loops
    # ===------------------------------------------------------------------=== #

    def as_face(self, b: BoundaryLoop) -> Face:
        """The face a boundary loop is stored as.

        Args:
            b: A boundary loop.

        Returns:
            Its face (`is_boundary_loop` of it is true).
        """
        return self._face_at(self._n_faces_capacity_count - 1 - self._check(b))

    def halfedge(self, b: BoundaryLoop) -> Halfedge:
        """A halfedge of a boundary loop (exterior).

        Args:
            b: A boundary loop.

        Returns:
            One of its halfedges.
        """
        return self.halfedge(self.as_face(b))

    def is_dead(self, b: BoundaryLoop) -> Bool:
        """Whether a boundary loop was deleted.

        Args:
            b: A boundary loop.

        Returns:
            True if it is no longer part of the mesh.
        """
        return self.is_dead(self.as_face(b))

    def degree(self, b: BoundaryLoop) -> Int:
        """The number of edges around a boundary loop.

        Args:
            b: A boundary loop.

        Returns:
            Its degree.
        """
        return self.degree(self.as_face(b))

    def adjacent_halfedges(
        ref self, b: BoundaryLoop
    ) -> NavigationIterator[FaceAdjacentHalfedgeNavigator, origin_of(self)]:
        """The halfedges around a boundary loop.

        Args:
            b: A boundary loop.

        Returns:
            An iterator over the halfedges.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentHalfedgeNavigator(self.halfedge(b))
        )

    def adjacent_vertices(
        ref self, b: BoundaryLoop
    ) -> NavigationIterator[FaceAdjacentVertexNavigator, origin_of(self)]:
        """The vertices around a boundary loop.

        Args:
            b: A boundary loop.

        Returns:
            An iterator over the vertices.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentVertexNavigator(self.halfedge(b))
        )

    def adjacent_edges(
        ref self, b: BoundaryLoop
    ) -> NavigationIterator[FaceAdjacentEdgeNavigator, origin_of(self)]:
        """The edges around a boundary loop.

        Args:
            b: A boundary loop.

        Returns:
            An iterator over the edges.
        """
        return NavigationIterator(
            Pointer(to=self), FaceAdjacentEdgeNavigator(self.halfedge(b))
        )

    # ===------------------------------------------------------------------=== #
    # Whole-mesh queries
    # ===------------------------------------------------------------------=== #

    def has_boundary(self) -> Bool:
        """Whether the mesh has boundary (is not closed).

        Returns:
            True if it has a boundary loop.
        """
        return self._n_boundary_loops_count > 0

    def is_triangular(self) -> Bool:
        """Whether every face is a triangle (O(n)).

        Returns:
            True for a triangle mesh.
        """
        for f in self.faces():
            if not self.is_triangle(f):
                return False
        return True

    def is_manifold(self) -> Bool:
        """Whether the mesh is manifold: always, for this type.

        Returns:
            True.
        """
        return True

    def is_edge_manifold(self) -> Bool:
        """Whether every edge is manifold: always, for this type.

        Returns:
            True.
        """
        return True

    def is_oriented(self) -> Bool:
        """Whether the faces are consistently oriented: always, for this
        type.

        Returns:
            True.
        """
        return True

    def euler_characteristic(self) -> Int:
        """`V - E + F`, counting boundary loops as faces.

        Returns:
            The Euler characteristic (2 for a sphere).
        """
        return (
            self.n_vertices()
            - self.n_edges()
            + self.n_faces()
            + self.n_boundary_loops()
        )

    def genus(self) -> Int:
        """The genus, `(2 - b - chi) / 2` with `b` boundary loops (rounded
        toward zero, as geometry-central's `int` division does; meaningful
        for a connected mesh).

        Returns:
            The genus (1 for a torus).
        """
        return _c_div2(
            2 - self.n_boundary_loops() - self.euler_characteristic()
        )

    def n_connected_components(self) -> Int:
        """The number of connected components (O(n)).

        Returns:
            The count.
        """
        var vertex_index = self.get_vertex_indices()
        var sets = DisjointSets(self.n_vertices())
        for e in self.edges():
            var h = self.halfedge(e)
            sets.merge(
                vertex_index[self.vertex(h)],
                vertex_index[self.vertex(self.twin(h))],
            )
        var root_seen = List[Bool](length=self.n_vertices() + 1, fill=False)
        var n = 0
        for i in range(self.n_vertices()):
            var r = sets.find(i)
            if not root_seen[r]:
                root_seen[r] = True
                n += 1
        return n

    def connecting_edge(self, a: Vertex, b: Vertex) -> Edge:
        """The edge between two vertices.

        Args:
            a: A vertex.
            b: Another vertex.

        Returns:
            The edge from `a` to `b`, or the null `Edge()` if there is none.
        """
        for e in self.adjacent_edges(a):
            if self.other_vertex(e, a) == b:
                return e
        return Edge()

    def get_vertex_indices(self) -> VertexData[Int]:
        """A dense numbering of the vertices, in index order (their indices,
        when the mesh is compressed).

        Returns:
            Each vertex's number.
        """
        var out = VertexData[Int](self, INVALID_IND)
        var i = 0
        for v in self.vertices():
            out[v] = i
            i += 1
        return out^

    def get_interior_vertex_indices(self) -> VertexData[Int]:
        """A dense numbering of the interior vertices, in index order.

        Returns:
            Each interior vertex's number, and `INVALID_IND` for boundary
            vertices.
        """
        var out = VertexData[Int](self, INVALID_IND)
        var i = 0
        for v in self.vertices():
            if not self.is_boundary(v):
                out[v] = i
                i += 1
        return out^

    def get_halfedge_indices(self) -> HalfedgeData[Int]:
        """A dense numbering of the halfedges, in index order.

        Returns:
            Each halfedge's number.
        """
        var out = HalfedgeData[Int](self, INVALID_IND)
        var i = 0
        for h in self.halfedges():
            out[h] = i
            i += 1
        return out^

    def get_corner_indices(self) -> CornerData[Int]:
        """A dense numbering of the corners, in index order.

        Returns:
            Each corner's number.
        """
        var out = CornerData[Int](self, INVALID_IND)
        var i = 0
        for c in self.corners():
            out[c] = i
            i += 1
        return out^

    def get_edge_indices(self) -> EdgeData[Int]:
        """A dense numbering of the edges, in index order.

        Returns:
            Each edge's number.
        """
        var out = EdgeData[Int](self, INVALID_IND)
        var i = 0
        for e in self.edges():
            out[e] = i
            i += 1
        return out^

    def get_face_indices(self) -> FaceData[Int]:
        """A dense numbering of the faces, in index order.

        Returns:
            Each face's number.
        """
        var out = FaceData[Int](self, INVALID_IND)
        var i = 0
        for f in self.faces():
            out[f] = i
            i += 1
        return out^

    def get_boundary_loop_indices(self) -> BoundaryLoopData[Int]:
        """A dense numbering of the boundary loops, in index order.

        Returns:
            Each boundary loop's number.
        """
        var out = BoundaryLoopData[Int](self, INVALID_IND)
        var i = 0
        for b in self.boundary_loops():
            out[b] = i
            i += 1
        return out^

    def get_vertex_boundary_status(self) -> VertexData[Bool]:
        """Which vertices are on the boundary.

        Returns:
            `is_boundary(v)` for each vertex.
        """
        var out = VertexData[Bool](self, False)
        for v in self.vertices():
            out[v] = self.is_boundary(v)
        return out^

    def get_face_vertex_list(self) -> List[List[Int]]:
        """The faces as lists of vertex numbers (`get_vertex_indices`), the
        input format of the constructor.

        Returns:
            Each face's vertices, counterclockwise from its halfedge's
            tail.
        """
        var vertex_index = self.get_vertex_indices()
        var out = List[List[Int]](capacity=self.n_faces())
        for f in self.faces():
            var face = List[Int]()
            for v in self.adjacent_vertices(f):
                face.append(vertex_index[v])
            out.append(face^)
        return out^

    def print_statistics(self):
        """Prints the element counts, as geometry-central does."""
        print("Halfedge mesh with: ")
        print("    # verts =  ", self.n_vertices(), sep="")
        print("    # edges =  ", self.n_edges(), sep="")
        print("    # faces =  ", self.n_faces(), sep="")
        print(
            "    # halfedges =  ",
            self.n_halfedges(),
            "  (",
            self.n_interior_halfedges(),
            " interior, ",
            self.n_exterior_halfedges(),
            " exterior)",
            sep="",
        )
        print(
            "      and ",
            self.n_boundary_loops(),
            " boundary components. ",
            sep="",
        )

    # ===------------------------------------------------------------------=== #
    # Validation
    # ===------------------------------------------------------------------=== #

    def validate_connectivity(self) raises:
        """Checks the halfedge structure for consistency (geometry-central's
        `validateConnectivity`, its implicit-twin checks).

        Raises:
            On the first inconsistency found, describing it.
        """
        if (
            self.n_interior_halfedges() + self.n_exterior_halfedges()
            != self.n_halfedges()
        ):
            raise Error("nInterior + nImaginary != nTotal halfedges")
        if self._n_halfedges_count > self._n_halfedges_fill_count:
            raise Error("halfedge count > halfedge fill")
        if self._n_halfedges_fill_count > self._n_halfedges_capacity_count:
            raise Error("halfedge fill > halfedge capacity")
        if self._n_vertices_count > self._n_vertices_fill_count:
            raise Error("vertex count > vertex fill")
        if self._n_vertices_fill_count > self._n_vertices_capacity_count:
            raise Error("vertex fill > vertex capacity")
        if self._n_faces_count > self._n_faces_fill_count:
            raise Error("face count > face fill")
        if (
            self._n_faces_fill_count + self._n_boundary_loops_count
            > self._n_faces_capacity_count
        ):
            raise Error("face + bl fill > face capacity")

        # == References
        for i in range(self._n_halfedges_fill_count):
            if Int(self._he_next_arr[i]) == INVALID_IND:
                continue
            self._validate_halfedge(i ^ 1, "he.twin()")
            self._validate_halfedge(Int(self._he_next_arr[i]), "he.next()")
            self._validate_vertex(Int(self._he_vertex_arr[i]), "he.vertex()")
            self._validate_edge(i // 2, "he.edge()")
            self._validate_face(Int(self._he_face_arr[i]), "he.face()")
        for i in range(self._n_vertices_fill_count):
            if Int(self._v_halfedge_arr[i]) == INVALID_IND:
                continue
            self._validate_halfedge(
                Int(self._v_halfedge_arr[i]), "v.halfedge()"
            )
        for i in range(self._n_edges_fill_count):
            if Int(self._he_next_arr[2 * i]) == INVALID_IND:
                continue
            self._validate_halfedge(2 * i, "e.halfedge()")
        for i in range(self._n_faces_fill_count):
            if Int(self._f_halfedge_arr[i]) == INVALID_IND:
                continue
            self._validate_halfedge(
                Int(self._f_halfedge_arr[i]), "f.halfedge()"
            )

        # == Edges and twins
        for h in self.halfedges():
            var cur = h
            var count = 0
            while True:
                if self.edge(cur) != self.edge(h):
                    raise Error(
                        "(he sibling) halfedge sibling doesn't have edge =="
                        " he.edge"
                    )
                if count > self.n_halfedges():
                    raise Error(
                        "(he sibling) halfedge sibling doesn't cycle back"
                    )
                cur = self.sibling(cur)
                count += 1
                if cur == h:
                    break
        for e in self.edges():
            for h in self.adjacent_halfedges(e):
                if e != self.edge(h):
                    raise Error("edge.halfedge doesn't match halfedge.edge")

        # == Faces and next
        for f in self.faces():
            if self.face(self.halfedge(f)) != f:
                raise Error("f.halfedge().face() is not f")
            var first = self.halfedge(f)
            var cur = first
            var count = 0
            while True:
                if self.face(cur) != f:
                    raise Error("face.halfedge doesn't match halfedge.face")
                cur = self.next(cur)
                count += 1
                if count > self._n_halfedges_count:
                    raise Error("next forms non-face loop")
                if cur == first:
                    break
            if count < 2:
                raise Error("face of degree < 2")

        # == Boundary loops and next
        for b in self.boundary_loops():
            if self.as_boundary_loop(self.as_face(b)) != b:
                raise Error("b.asFace().asBoundaryLoop() is not fixed point")
            if self.face(self.halfedge(b)) != self.as_face(b):
                raise Error("bl.halfedge().face() is not bl")
            if not self.is_boundary_loop(self.face(self.halfedge(b))):
                raise Error("bl.halfedge().face() is not a boundary loop")
            var first = self.halfedge(b)
            var cur = first
            var count = 0
            while True:
                if not self.is_boundary_loop(self.face(cur)):
                    raise Error(
                        "walking around boundary loop yielded he.face() which"
                        " is not a boundary loop"
                    )
                if self.as_boundary_loop(self.face(cur)) != b:
                    raise Error(
                        "(boundary loop) face.halfedge doesn't match"
                        " halfedge.face"
                    )
                cur = self.next(cur)
                count += 1
                if count > self._n_halfedges_count:
                    raise Error("(boundary loop) next forms non-face loop")
                if cur == first:
                    break
            if count < 2:
                raise Error("(boundary loop) face of degree < 2")

        # == Halfedge neighborhoods
        var seen = List[Bool](
            length=self._n_halfedges_capacity_count, fill=False
        )
        for h in self.halfedges():
            if not self.is_interior(h):
                if h == self.halfedge(self.edge(h)):
                    raise Error("exterior halfedge is e.halfedge()")
                if not self.is_interior(self.twin(h)):
                    raise Error("he and he.twin() are both exterior")
            if self.vertex(h) != self.vertex(self.next(self.twin(h))):
                raise Error("halfedge vertices don't match")
            if seen[h.index()]:
                continue
            var cur = h
            var count = 0
            while True:
                if self.face(cur) != self.face(h):
                    raise Error("he.next.**.face doesn't match he.face")
                seen[cur.index()] = True
                cur = self.next(cur)
                count += 1
                if count > self._n_halfedges_count:
                    raise Error("next forms non-face loop")
                if cur == h:
                    break

        # == Vertex orbits
        for v in self.vertices():
            var count = 0
            for h in self.outgoing_halfedges(v):
                if count > self._n_halfedges_count:
                    raise Error("vertex outgoing halfedges has bad cycle")
                if self.vertex(h) != v:
                    raise Error("vertex.halfedge doesn't match halfedge.vertex")
                count += 1

        # == Boundary rules
        for v in self.vertices():
            var boundary_count = 0
            for h in self.outgoing_halfedges(v):
                if not self.is_interior(h):
                    boundary_count += 1
            if boundary_count > 1:
                raise Error("multiple boundaries incident on vertex")
            if boundary_count > 0:
                if not self.is_interior(self.halfedge(v)):
                    raise Error("v.halfedge() is exterior")
                if self.is_interior(self.twin(self.halfedge(v))):
                    raise Error(
                        "v.halfedge() does not border boundary on a boundary"
                        " vertex"
                    )
                if not self.is_boundary(v):
                    raise Error("computed v.isBoundary is wrong")

    def _validate_vertex(self, i: Int, msg: StaticString) raises:
        if (
            i < 0
            or i >= self._n_vertices_fill_count
            or Int(self._v_halfedge_arr[i]) == INVALID_IND
        ):
            raise Error(msg, " - bad vertex reference")

    def _validate_halfedge(self, i: Int, msg: StaticString) raises:
        if (
            i < 0
            or i >= self._n_halfedges_fill_count
            or Int(self._he_next_arr[i]) == INVALID_IND
        ):
            raise Error(msg, " - bad halfedge reference")

    def _validate_edge(self, i: Int, msg: StaticString) raises:
        if (
            i < 0
            or i >= self._n_edges_fill_count
            or Int(self._he_next_arr[2 * i]) == INVALID_IND
        ):
            raise Error(msg, " - bad edge reference")

    def _validate_face(self, i: Int, msg: StaticString) raises:
        if (
            i < 0
            or i >= self._n_faces_capacity_count
            or Int(self._f_halfedge_arr[i]) == INVALID_IND
            or (
                i >= self._n_faces_fill_count
                and i
                < self._n_faces_capacity_count - self._n_boundary_loops_count
            )
        ):
            raise Error(msg, " - bad face reference")
