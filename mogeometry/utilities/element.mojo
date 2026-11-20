# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Follows geometry-central's include/geometrycentral/utilities/element.h and
# element.ipp (commit 019669d, MIT License, Copyright (c) 2017-2019 Nicholas
# Sharp and the geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Element handles (`mogeometry.utilities.element`)
================================================
`ElementHandle[kind]`, the handle type behind `Vertex`, `Halfedge`,
`Corner`, `Edge`, `Face`, and `BoundaryLoop`: the index of an element in
its mesh, as a distinct type per kind of element. The mesh navigates
(`mesh.next(h)`, `mesh.vertex(h)`); a handle holds no reference to it.

Debug builds tag handles with their mesh. A mesh returns handles whose
`Int` carries, above the index, the mesh's id and generation (see
`MeshTag`), and its methods check the tag of every handle they are given:
a handle from another mesh, or from before a `compress()`, stops the
program with an assertion. `index()` masks the tag off. A handle made
from a plain integer (`Vertex(5)`) has no tag and is not checked. Built
with `-D ASSERT=none`, handles are plain indices and nothing is checked.

Also here: the traits `MeshElement` (what `MeshData` needs of a handle)
and `ElementMesh` (what it needs of a mesh), and `MESH_INDEX`, the
integer type meshes store their connectivity in.

Differences from geometry-central:

- A handle is an index only; geometry-central's `Element` also holds a
  pointer to its mesh, and navigates itself (`he.twin().vertex()`). See
  AGENTS.md (Architecture) for why.
- `getIndex()` is `index()`; the null handle `Vertex()` has index
  `INVALID_IND`, which is -1 (geometry-central's is the largest
  `size_t`), and prints as `v_-1`.
"""

from std.hash import Hasher
from std.sys import get_defined_string

comptime MESH_INDEX = DType.int32
"""The integer type meshes store connectivity in (handles are `Int`). One
constant for every mesh; `DType.int64` lifts the limit of 2^31 halfedges
at the cost of twice the memory traffic (AGENTS.md, Architecture)."""

comptime INVALID_IND = -1
"""The index of no element (the null handle's), also stored in
connectivity arrays where there is no element."""

comptime MESH_CHECKS = get_defined_string["ASSERT", "safe"]() != "none"
"""Whether handles carry mesh tags and meshes check them: true unless
built with `-D ASSERT=none`."""

comptime INDEX_BITS = 32 if MESH_INDEX == DType.int32 else 44
"""The low bits of a tagged handle that hold its index."""

comptime INDEX_MASK = (1 << INDEX_BITS) - 1
"""Masks a tagged handle down to its index."""

comptime GENERATION_BITS = 7
"""The bits of a tag that hold the mesh's generation (its `compress()`
count, modulo 128)."""

comptime ID_BITS = 63 - INDEX_BITS - GENERATION_BITS
"""The bits of a tag that hold the mesh's id."""


struct ElementKind:
    """The kinds of mesh element, as `ElementHandle`'s parameter."""

    comptime VERTEX = 0
    """A vertex (`Vertex`)."""
    comptime HALFEDGE = 1
    """A halfedge (`Halfedge`)."""
    comptime CORNER = 2
    """A corner of a face (`Corner`): the halfedge's index, interior only."""
    comptime EDGE = 3
    """An edge (`Edge`)."""
    comptime FACE = 4
    """A face (`Face`), or a boundary loop seen as a face."""
    comptime BOUNDARY_LOOP = 5
    """A boundary loop (`BoundaryLoop`)."""

    @staticmethod
    def short_name(kind: Int) -> StaticString:
        """The prefix a handle of this kind prints with (`v` for `v_12`).

        Args:
            kind: An `ElementKind` value.

        Returns:
            The short name geometry-central gives the kind.
        """
        if kind == ElementKind.VERTEX:
            return "v"
        if kind == ElementKind.HALFEDGE:
            return "he"
        if kind == ElementKind.CORNER:
            return "c"
        if kind == ElementKind.EDGE:
            return "e"
        if kind == ElementKind.FACE:
            return "f"
        return "bl"


struct MeshTag:
    """Tags: a mesh's id and generation, as debug builds store them above
    a handle's index.

    A tag is `id | generation << ID_BITS`, nonzero for every mesh (ids
    start at 1), and 0 on a handle made from a plain integer. A mesh's id
    is a hash of the input it was built from, so a copy keeps it (a copy
    has the same elements); two meshes built from the same input share
    it, harmlessly.
    """

    @staticmethod
    def make(id: Int, generation: Int) -> Int:
        """The tag of a mesh.

        Args:
            id: The mesh's id, in `1 .. 2^ID_BITS - 1`.
            generation: The number of times the mesh has been compressed.

        Returns:
            The tag.
        """
        return id | ((generation % (1 << GENERATION_BITS)) << ID_BITS)

    @staticmethod
    def id_from_hash(hash: UInt64) -> Int:
        """A mesh id from a hash of its input.

        Args:
            hash: Any 64-bit hash.

        Returns:
            An id in `1 .. 2^ID_BITS - 1` (never 0, the untagged value).
        """
        return Int(hash % UInt64((1 << ID_BITS) - 1)) + 1


trait MeshElement(Equatable, TrivialRegisterPassable):
    """A handle that `MeshData` can index by: an element's index, its kind,
    and the tag of the mesh it came from."""

    comptime KIND: Int
    """The `ElementKind` of the handle."""

    def index(self) -> Int:
        """The element's index in its mesh.

        Returns:
            The index, or `INVALID_IND`.
        """
        ...

    def tag(self) -> Int:
        """The tag of the mesh the handle came from.

        Returns:
            The tag, or 0 if the handle is untagged (made from an integer,
            or any handle in a build with `-D ASSERT=none`).
        """
        ...


trait ElementMesh:
    """A mesh that `MeshData` can be built on: how many elements of each
    kind it has room for, and its tag."""

    def element_capacity(self, kind: Int) -> Int:
        """The number of elements of a kind a container must hold.

        Args:
            kind: An `ElementKind` value.

        Returns:
            The capacity (geometry-central's `elementCapacity<E>`).
        """
        ...

    def mesh_tag(self) -> Int:
        """The mesh's tag (`MeshTag`), which its handles carry in debug
        builds.

        Returns:
            The tag.
        """
        ...


struct ElementHandle[kind: Int](
    Comparable, Hashable, MeshElement, TrivialRegisterPassable, Writable
):
    """The handle of a mesh element: its index, with the mesh's tag above
    it in debug builds.

    Use the aliases `Vertex`, `Halfedge`, `Corner`, `Edge`, `Face`, and
    `BoundaryLoop` (`mogeometry.surface`). Handles compare, order, and
    hash by index, so handles from different meshes, or a tagged and an
    untagged one, compare equal when their indices are.

    Parameters:
        kind: The `ElementKind` of the element.
    """

    comptime KIND = Self.kind
    """The `ElementKind` of the handle."""

    var _raw: Int
    """The index, with the mesh tag above `INDEX_BITS` in debug builds;
    negative for the null handle."""

    @always_inline
    def __init__(out self):
        """The null handle, which refers to no element."""
        self._raw = INVALID_IND

    @always_inline
    def __init__(out self, index: Int):
        """An untagged handle to the element with this index.

        Args:
            index: The element's index in its mesh.
        """
        self._raw = index

    @always_inline
    def __init__(out self, *, index: Int, tag: Int):
        """A handle carrying a mesh tag (in debug builds; the tag is
        dropped in builds with `-D ASSERT=none`).

        Args:
            index: The element's index, or `INVALID_IND`.
            tag: The mesh's tag.
        """
        comptime if MESH_CHECKS:
            self._raw = index | (tag << INDEX_BITS) if index >= 0 else index
        else:
            self._raw = index

    @always_inline
    def index(self) -> Int:
        """The element's index in its mesh.

        Returns:
            The index, or `INVALID_IND` for the null handle.
        """
        comptime if MESH_CHECKS:
            return self._raw & INDEX_MASK if self._raw >= 0 else self._raw
        else:
            return self._raw

    @always_inline
    def tag(self) -> Int:
        """The tag of the mesh the handle came from.

        Returns:
            The tag, or 0 if the handle is untagged.
        """
        comptime if MESH_CHECKS:
            return self._raw >> INDEX_BITS if self._raw >= 0 else 0
        else:
            return 0

    @always_inline
    def is_valid(self) -> Bool:
        """Whether the handle refers to an element (is not the null
        handle).

        Returns:
            True unless the index is `INVALID_IND`.
        """
        return self._raw >= 0

    @always_inline
    def __eq__(self, other: Self) -> Bool:
        """Equal indices.

        Args:
            other: The handle to compare with.

        Returns:
            True if both handles have the same index.
        """
        return self.index() == other.index()

    @always_inline
    def __ne__(self, other: Self) -> Bool:
        """Different indices.

        Args:
            other: The handle to compare with.

        Returns:
            True if the handles have different indices.
        """
        return self.index() != other.index()

    @always_inline
    def __lt__(self, other: Self) -> Bool:
        """Orders by index.

        Args:
            other: The handle to compare with.

        Returns:
            True if this handle's index is smaller.
        """
        return self.index() < other.index()

    @always_inline
    def __le__(self, other: Self) -> Bool:
        """Orders by index.

        Args:
            other: The handle to compare with.

        Returns:
            True if this handle's index is smaller or equal.
        """
        return self.index() <= other.index()

    @always_inline
    def __gt__(self, other: Self) -> Bool:
        """Orders by index.

        Args:
            other: The handle to compare with.

        Returns:
            True if this handle's index is larger.
        """
        return self.index() > other.index()

    @always_inline
    def __ge__(self, other: Self) -> Bool:
        """Orders by index.

        Args:
            other: The handle to compare with.

        Returns:
            True if this handle's index is larger or equal.
        """
        return self.index() >= other.index()

    def __hash__[H: Hasher](self, mut hasher: H):
        """Hashes the index.

        Parameters:
            H: The hasher type.

        Args:
            hasher: The hasher to update.
        """
        self.index().__hash__(hasher)

    def write_to(self, mut writer: Some[Writer]):
        """Writes the handle as geometry-central prints it: the kind's
        short name and the index, `v_12` or `he_3`.

        Args:
            writer: The writer to write to.
        """
        writer.write(ElementKind.short_name(Self.kind), "_", self.index())
