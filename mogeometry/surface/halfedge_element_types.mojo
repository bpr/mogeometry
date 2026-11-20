# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Follows geometry-central's include/geometrycentral/surface/
# halfedge_element_types.h and .ipp, and utilities/element_iterators.h and
# .ipp (commit 019669d, MIT License, Copyright (c) 2017-2019 Nicholas Sharp
# and the geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Halfedge mesh elements (`mogeometry.surface.halfedge_element_types`)
====================================================================
The handles of a surface mesh's elements (`Vertex`, `Halfedge`, `Corner`,
`Edge`, `Face`, `BoundaryLoop`), the containers of values on them
(`VertexData[T]`, ...), and the iterators that `ManifoldSurfaceMesh`
returns: over all elements of a kind (`mesh.vertices()`) and around an
element (`mesh.outgoing_halfedges(v)`, `mesh.adjacent_faces(f)`).

Differences from geometry-central:

- Navigation is on the mesh, not the handle: `he.twin().vertex()` is
  `mesh.vertex(mesh.twin(he))`, `v.outgoingHalfedges()` is
  `mesh.outgoing_halfedges(v)` (see `ManifoldSurfaceMesh`).
- One iterator type, `NavigationIterator[N]`, runs every navigator `N`
  (geometry-central's `NavigationSetBase<N>`), and one, `ElementRange`,
  every element range (its `RangeSetBase<F>`). The boundary-loop
  navigators are the face ones started from the loop's halfedge.
- Only the manifold mesh exists, so the vertex navigators orbit by
  `twin().next()` alone (geometry-central's `VertexNeighborIteratorState`
  also handles nonmanifold vertices).
- `Edge::diamondBoundary()` returns a `DiamondBoundary` with named fields,
  and `Edge::adjacentVertices()` an `EdgeVertices`, instead of
  `std::array`s.
"""

from mogeometry.surface.manifold_surface_mesh import ManifoldSurfaceMesh
from mogeometry.utilities.element import ElementHandle, ElementKind
from mogeometry.utilities.mesh_data import MeshData

comptime Vertex = ElementHandle[ElementKind.VERTEX]
"""A vertex of a surface mesh."""
comptime Halfedge = ElementHandle[ElementKind.HALFEDGE]
"""A halfedge: one side of an edge, directed, in one face or boundary
loop."""
comptime Corner = ElementHandle[ElementKind.CORNER]
"""A corner of a face, where an interior halfedge leaves its vertex (the
corner and the halfedge share an index)."""
comptime Edge = ElementHandle[ElementKind.EDGE]
"""An edge: a pair of twin halfedges."""
comptime Face = ElementHandle[ElementKind.FACE]
"""A face, or a boundary loop seen as a face (`mesh.is_boundary_loop(f)`;
such a face cannot index a `FaceData`)."""
comptime BoundaryLoop = ElementHandle[ElementKind.BOUNDARY_LOOP]
"""A boundary loop: a cycle of exterior halfedges, one per hole."""

comptime VertexData[T: Copyable & Deinitable] = MeshData[Vertex, T]
"""A value of type `T` on each vertex.

Parameters:
    T: The value type.
"""
comptime HalfedgeData[T: Copyable & Deinitable] = MeshData[Halfedge, T]
"""A value of type `T` on each halfedge.

Parameters:
    T: The value type.
"""
comptime CornerData[T: Copyable & Deinitable] = MeshData[Corner, T]
"""A value of type `T` on each corner (indexed like halfedges).

Parameters:
    T: The value type.
"""
comptime EdgeData[T: Copyable & Deinitable] = MeshData[Edge, T]
"""A value of type `T` on each edge.

Parameters:
    T: The value type.
"""
comptime FaceData[T: Copyable & Deinitable] = MeshData[Face, T]
"""A value of type `T` on each face.

Parameters:
    T: The value type.
"""
comptime BoundaryLoopData[T: Copyable & Deinitable] = MeshData[BoundaryLoop, T]
"""A value of type `T` on each boundary loop.

Parameters:
    T: The value type.
"""


@fieldwise_init
struct DiamondBoundary(TrivialRegisterPassable):
    """The four halfedges around the two triangles of an interior edge,
    as `mesh.diamond_boundary(e)` returns them: for `h = halfedge(e)` and
    `t = twin(h)`, `next(h)`, `next(next(h))`, `next(t)`, `next(next(t))`.
    """

    var next: Halfedge
    """`next(h)`."""
    var next_next: Halfedge
    """`next(next(h))`."""
    var twin_next: Halfedge
    """`next(t)`."""
    var twin_next_next: Halfedge
    """`next(next(t))`."""

    def __getitem__(self, i: Int) -> Halfedge:
        """The halfedges in geometry-central's order.

        Args:
            i: 0 to 3.

        Returns:
            `next`, `next_next`, `twin_next`, or `twin_next_next`.
        """
        debug_assert[assert_mode="safe"](0 <= i < 4, "index out of range")
        if i == 0:
            return self.next
        if i == 1:
            return self.next_next
        if i == 2:
            return self.twin_next
        return self.twin_next_next


@fieldwise_init
struct EdgeVertices(TrivialRegisterPassable):
    """The two ends of an edge, as `mesh.adjacent_vertices(e)` returns
    them."""

    var first: Vertex
    """The tail of `halfedge(e)`."""
    var second: Vertex
    """The tip of `halfedge(e)`."""


# ===----------------------------------------------------------------------=== #
# Navigation around an element
# ===----------------------------------------------------------------------=== #


trait Navigator(TrivialRegisterPassable):
    """The state of a walk around an element: how to step, which states
    to yield, and what to yield (geometry-central's navigator structs).
    """

    comptime Element: Copyable
    """What the walk yields."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next state.

        Args:
            mesh: The mesh walked on.
        """
        ...

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Whether the current state yields an element.

        Args:
            mesh: The mesh walked on.

        Returns:
            True if the walk yields at this state.
        """
        ...

    def current(self, mesh: ManifoldSurfaceMesh) -> Self.Element:
        """The element the current state yields.

        Args:
            mesh: The mesh walked on.

        Returns:
            The element.
        """
        ...

    def same(self, other: Self) -> Bool:
        """Whether two states are the same (the walk ends when it returns
        to its first valid state).

        Args:
            other: The state to compare with.

        Returns:
            True if the states are equal.
        """
        ...


struct NavigationIterator[N: Navigator, origin: ImmOrigin](Copyable, Iterator):
    """Walks a navigator from its first valid state around to it again,
    yielding at each valid state (geometry-central's `NavigationSetBase`).

    Parameters:
        N: The navigator.
        origin: The origin of the mesh.
    """

    comptime Element = Self.N.Element
    """What the walk yields."""

    var _mesh: Pointer[ManifoldSurfaceMesh, Self.origin]
    """The mesh walked on."""
    var _state: Self.N
    """The state to yield next."""
    var _first: Self.N
    """The first valid state, where the walk stops."""
    var _done: Bool
    """Whether the walk is over."""

    def __init__(
        out self, mesh: Pointer[ManifoldSurfaceMesh, Self.origin], start: Self.N
    ):
        """A walk from `start`.

        Args:
            mesh: The mesh walked on.
            start: The initial state (the walk begins at the first valid
                state from it; if there is none, it yields nothing).
        """
        self._mesh = mesh
        self._done = False
        var state = start
        while not state.is_valid(mesh[]):
            state.advance(mesh[])
            if state.same(start):
                self._done = True
                break
        self._state = state
        self._first = state

    def __iter__(ref self) -> Self:
        """The iterator itself.

        Returns:
            A copy of the iterator.
        """
        return self.copy()

    def __next__(mut self) raises StopIteration -> Self.Element:
        """The next element of the walk.

        Returns:
            The element.

        Raises:
            StopIteration: When the walk is back at its first state.
        """
        if self._done:
            raise StopIteration()
        var out = self._state.current(self._mesh[])
        self._state.advance(self._mesh[])
        while not self._state.is_valid(self._mesh[]):
            self._state.advance(self._mesh[])
        if self._state.same(self._first):
            self._done = True
        return out^


@fieldwise_init
struct VertexOutgoingHalfedgeNavigator(Navigator):
    """The halfedges leaving a vertex, in order (`twin().next()`)."""

    comptime Element = Halfedge
    """Yields halfedges."""
    var he: Halfedge
    """The current outgoing halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next outgoing halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_outgoing_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every outgoing halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Halfedge:
        """The halfedge.

        Args:
            mesh: The mesh.

        Returns:
            The current outgoing halfedge.
        """
        return self.he

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct VertexIncomingHalfedgeNavigator(Navigator):
    """The halfedges arriving at a vertex, in order (`next().twin()`)."""

    comptime Element = Halfedge
    """Yields halfedges."""
    var he: Halfedge
    """The current incoming halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next incoming halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_incoming_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every incoming halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Halfedge:
        """The halfedge.

        Args:
            mesh: The mesh.

        Returns:
            The current incoming halfedge.
        """
        return self.he

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct VertexAdjacentVertexNavigator(Navigator):
    """The neighbors of a vertex: the tips of its outgoing halfedges."""

    comptime Element = Vertex
    """Yields vertices."""
    var he: Halfedge
    """The current outgoing halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next outgoing halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_outgoing_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every outgoing halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Vertex:
        """The halfedge's tip.

        Args:
            mesh: The mesh.

        Returns:
            The neighbor.
        """
        return mesh.tip_vertex(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct VertexAdjacentCornerNavigator(Navigator):
    """The corners at a vertex: its interior outgoing halfedges'."""

    comptime Element = Corner
    """Yields corners."""
    var he: Halfedge
    """The current outgoing halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next outgoing halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_outgoing_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Interior halfedges only.

        Args:
            mesh: The mesh.

        Returns:
            True if the halfedge is in a face.
        """
        return mesh.is_interior(self.he)

    def current(self, mesh: ManifoldSurfaceMesh) -> Corner:
        """The halfedge's corner.

        Args:
            mesh: The mesh.

        Returns:
            The corner.
        """
        return mesh.corner(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct VertexAdjacentEdgeNavigator(Navigator):
    """The edges at a vertex: its outgoing halfedges'."""

    comptime Element = Edge
    """Yields edges."""
    var he: Halfedge
    """The current outgoing halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next outgoing halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_outgoing_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every outgoing halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Edge:
        """The halfedge's edge.

        Args:
            mesh: The mesh.

        Returns:
            The edge.
        """
        return mesh.edge(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct VertexAdjacentFaceNavigator(Navigator):
    """The faces at a vertex: its interior outgoing halfedges'."""

    comptime Element = Face
    """Yields faces."""
    var he: Halfedge
    """The current outgoing halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next outgoing halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next_outgoing_neighbor(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Interior halfedges only.

        Args:
            mesh: The mesh.

        Returns:
            True if the halfedge is in a face.
        """
        return mesh.is_interior(self.he)

    def current(self, mesh: ManifoldSurfaceMesh) -> Face:
        """The halfedge's face.

        Args:
            mesh: The mesh.

        Returns:
            The face.
        """
        return mesh.face(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct EdgeAdjacentHalfedgeNavigator(Navigator):
    """The halfedges of an edge (`twin()`)."""

    comptime Element = Halfedge
    """Yields halfedges."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the twin.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.sibling(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Both halfedges.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Halfedge:
        """The halfedge.

        Args:
            mesh: The mesh.

        Returns:
            The current halfedge.
        """
        return self.he

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct EdgeAdjacentInteriorHalfedgeNavigator(Navigator):
    """The interior halfedges of an edge."""

    comptime Element = Halfedge
    """Yields halfedges."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the twin.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.sibling(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Interior halfedges only.

        Args:
            mesh: The mesh.

        Returns:
            True if the halfedge is in a face.
        """
        return mesh.is_interior(self.he)

    def current(self, mesh: ManifoldSurfaceMesh) -> Halfedge:
        """The halfedge.

        Args:
            mesh: The mesh.

        Returns:
            The current halfedge.
        """
        return self.he

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct EdgeAdjacentFaceNavigator(Navigator):
    """The faces on either side of an edge."""

    comptime Element = Face
    """Yields faces."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the twin.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.sibling(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Interior halfedges only.

        Args:
            mesh: The mesh.

        Returns:
            True if the halfedge is in a face.
        """
        return mesh.is_interior(self.he)

    def current(self, mesh: ManifoldSurfaceMesh) -> Face:
        """The halfedge's face.

        Args:
            mesh: The mesh.

        Returns:
            The face.
        """
        return mesh.face(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct FaceAdjacentHalfedgeNavigator(Navigator):
    """The halfedges around a face or boundary loop (`next()`)."""

    comptime Element = Halfedge
    """Yields halfedges."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Halfedge:
        """The halfedge.

        Args:
            mesh: The mesh.

        Returns:
            The current halfedge.
        """
        return self.he

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct FaceAdjacentVertexNavigator(Navigator):
    """The vertices around a face or boundary loop: its halfedges' tails."""

    comptime Element = Vertex
    """Yields vertices."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Vertex:
        """The halfedge's tail.

        Args:
            mesh: The mesh.

        Returns:
            The vertex.
        """
        return mesh.vertex(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct FaceAdjacentCornerNavigator(Navigator):
    """The corners of a face."""

    comptime Element = Corner
    """Yields corners."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Corner:
        """The halfedge's corner.

        Args:
            mesh: The mesh.

        Returns:
            The corner.
        """
        return mesh.corner(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct FaceAdjacentEdgeNavigator(Navigator):
    """The edges around a face or boundary loop."""

    comptime Element = Edge
    """Yields edges."""
    var he: Halfedge
    """The current halfedge."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the next halfedge.

        Args:
            mesh: The mesh.
        """
        self.he = mesh.next(self.he)

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Every halfedge.

        Args:
            mesh: The mesh.

        Returns:
            True.
        """
        return True

    def current(self, mesh: ManifoldSurfaceMesh) -> Edge:
        """The halfedge's edge.

        Args:
            mesh: The mesh.

        Returns:
            The edge.
        """
        return mesh.edge(self.he)

    def same(self, other: Self) -> Bool:
        """Same halfedge.

        Args:
            other: The other state.

        Returns:
            True if the halfedges are equal.
        """
        return self.he == other.he


@fieldwise_init
struct FaceAdjacentFaceNavigator(Navigator):
    """The faces across the edges of a face (a face shows up once per
    edge it shares)."""

    comptime Element = Face
    """Yields faces."""
    var he: Halfedge
    """The current halfedge of this face."""
    var other: Halfedge
    """The current halfedge of the edge, possibly `he` itself."""

    def advance(mut self, mesh: ManifoldSurfaceMesh):
        """Steps to the twin, or to the next edge once back at `he`.

        Args:
            mesh: The mesh.
        """
        self.other = mesh.sibling(self.other)
        if self.he == self.other:
            self.he = mesh.next(self.he)
            self.other = self.he

    def is_valid(self, mesh: ManifoldSurfaceMesh) -> Bool:
        """Interior halfedges across the edge.

        Args:
            mesh: The mesh.

        Returns:
            True if `other` is not `he` and is in a face.
        """
        return self.he != self.other and mesh.is_interior(self.other)

    def current(self, mesh: ManifoldSurfaceMesh) -> Face:
        """The face across.

        Args:
            mesh: The mesh.

        Returns:
            The face of `other`.
        """
        return mesh.face(self.other)

    def same(self, other: Self) -> Bool:
        """Same pair of halfedges.

        Args:
            other: The other state.

        Returns:
            True if both halfedges are equal.
        """
        return self.he == other.he and self.other == other.other


# ===----------------------------------------------------------------------=== #
# Ranges of elements
# ===----------------------------------------------------------------------=== #


struct RangeFilter:
    """Which elements of a kind an `ElementRange` yields."""

    comptime ALL = 0
    """Every live element."""
    comptime INTERIOR = 1
    """Live halfedges in faces."""
    comptime EXTERIOR = 2
    """Live halfedges in boundary loops."""


struct ElementRange[kind: Int, filter: Int, origin: ImmOrigin](
    Copyable, Iterator
):
    """The live elements of a kind, by index (geometry-central's
    `RangeSetBase`): `mesh.vertices()`, `mesh.interior_halfedges()`, ....

    Parameters:
        kind: The `ElementKind`.
        filter: The `RangeFilter`.
        origin: The origin of the mesh.
    """

    comptime Element = ElementHandle[Self.kind]
    """Yields handles of this kind."""

    var _mesh: Pointer[ManifoldSurfaceMesh, Self.origin]
    """The mesh."""
    var _i: Int
    """The next index to consider."""
    var _end: Int
    """One past the last index."""
    var _dense: Bool
    """Whether every index is a live element (a compressed mesh, and no
    filter), so none needs checking."""

    def __init__(
        out self, mesh: Pointer[ManifoldSurfaceMesh, Self.origin], end: Int
    ):
        """The range of indices `0 ..< end`.

        Args:
            mesh: The mesh.
            end: One past the last index (the mesh's index size).
        """
        self._mesh = mesh
        self._i = 0
        self._end = end
        self._dense = (
            Self.filter == RangeFilter.ALL
            and Self.kind != ElementKind.CORNER
            and mesh[].is_compressed()
        )

    def __iter__(ref self) -> Self:
        """The iterator itself.

        Returns:
            A copy of the iterator.
        """
        return self.copy()

    def __next__(mut self) raises StopIteration -> Self.Element:
        """The next live element.

        Returns:
            Its handle.

        Raises:
            StopIteration: After the last one.
        """
        if not self._dense:
            while self._i < self._end and not self._mesh[]._element_okay[
                Self.kind, Self.filter
            ](self._i):
                self._i += 1
        if self._i >= self._end:
            raise StopIteration()
        var i = self._i
        self._i += 1
        return Self.Element(index=i, tag=self._mesh[].mesh_tag())
