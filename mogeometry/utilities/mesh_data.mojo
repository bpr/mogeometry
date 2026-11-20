# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Follows geometry-central's include/geometrycentral/utilities/mesh_data.h
# and mesh_data.ipp (commit 019669d, MIT License, Copyright (c) 2017-2019
# Nicholas Sharp and the geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Data on mesh elements (`mogeometry.utilities.mesh_data`)
========================================================
`MeshData[E, T]`, a value of type `T` for every element of kind `E` (a
vertex, a face, ...) of a mesh, indexed by handle: `data[v]`. The
aliases `VertexData[T]`, `FaceData[T]`, ... are in `mogeometry.surface`.

A container holds no reference to its mesh (AGENTS.md, Architecture).
It records the mesh's tag (`MeshTag`), and after the mesh changes,
`data.sync(mesh)` brings it up to date. In debug builds, `data[e]`
checks that `e` carries the same tag, so a handle from another mesh, or
a container not synced since the mesh was compressed, stops the program
instead of reading the wrong value.

Differences from geometry-central:

- No callbacks: geometry-central's containers register with their mesh,
  which resizes and permutes them as it changes. Here growth and
  renumbering wait for `sync()`.
- The storage is a `List[T]` (`raw()`), not an Eigen vector; there is no
  `toVector`/`fromVector` (build from a `List` instead) and no
  `reinterpretTo`.
- A default value is always given (`VertexData[Float64](mesh, 0.0)`):
  Mojo has no `T()` for every `T`.
- No elementwise arithmetic operators (`A + B`, `A *= 2`) yet; loop over
  `raw()`.
"""

from mogeometry.utilities.element import (
    MESH_CHECKS,
    ElementMesh,
    MeshElement,
)


struct MeshData[E: MeshElement, T: Copyable & Deinitable](Copyable, Sized):
    """A value of type `T` for each element of kind `E` of a mesh.

    Parameters:
        E: The element handle type (`Vertex`, `Face`, ...).
        T: The value type.
    """

    var _data: List[Self.T]
    """The values, by element index."""
    var _default: Self.T
    """The value new elements get when `sync()` grows the container."""
    var _tag: Int
    """The tag of the mesh (and generation) the container matches."""

    def __init__[M: ElementMesh](out self, mesh: M, default: Self.T):
        """A container with every value set to `default`.

        Parameters:
            M: The mesh type.

        Args:
            mesh: The mesh whose elements the values belong to.
            default: The initial value, also given to elements added later.
        """
        self._data = List[Self.T](
            length=mesh.element_capacity(Self.E.KIND), fill=default
        )
        self._default = default.copy()
        self._tag = mesh.mesh_tag()

    def __init__[
        M: ElementMesh
    ](out self, mesh: M, var values: List[Self.T], default: Self.T) raises:
        """A container holding `values`, one per element index.

        Parameters:
            M: The mesh type.

        Args:
            mesh: The mesh whose elements the values belong to.
            values: The values, by element index.
            default: The value given to elements added later.

        Raises:
            If `values` does not have one value per element index.
        """
        var n = mesh.element_capacity(Self.E.KIND)
        if len(values) != n:
            raise Error(
                String(
                    "MeshData: ", len(values), " values for ", n, " elements"
                )
            )
        self._data = values^
        self._default = default.copy()
        self._tag = mesh.mesh_tag()

    @always_inline
    def __getitem__(ref self, e: Self.E) -> ref[self._data] Self.T:
        """The value of an element.

        Args:
            e: The element's handle.

        Returns:
            A reference to its value.
        """
        comptime if MESH_CHECKS:
            debug_assert[assert_mode="safe"](
                e.tag() == 0 or e.tag() == self._tag,
                "MeshData: handle from another mesh, or the container",
                " was not synced after the mesh was compressed",
            )
        var i = e.index()
        debug_assert[assert_mode="safe"](
            0 <= i < len(self._data), "MeshData: element index out of range"
        )
        return self._data.unsafe_ptr()[unsafe_offset=i]

    @always_inline
    def __getitem__(ref self, i: Int) -> ref[self._data] Self.T:
        """The value at an element index.

        Args:
            i: The element's index.

        Returns:
            A reference to its value.
        """
        debug_assert[assert_mode="safe"](
            0 <= i < len(self._data), "MeshData: element index out of range"
        )
        return self._data.unsafe_ptr()[unsafe_offset=i]

    def __len__(self) -> Int:
        """The number of values (the mesh's capacity for `E`, which is its
        element count unless it has deleted elements).

        Returns:
            The number of values.
        """
        return len(self._data)

    def raw(ref self) -> ref[self._data] List[Self.T]:
        """The values, by element index.

        Returns:
            A reference to the list of values.
        """
        return self._data

    def fill(mut self, value: Self.T):
        """Sets every value.

        Args:
            value: The value.
        """
        for i in range(len(self._data)):
            self._data[i] = value.copy()

    def default(self) -> Self.T:
        """The value elements added to the mesh get.

        Returns:
            A copy of the default value.
        """
        return self._default.copy()

    def set_default(mut self, value: Self.T):
        """Sets the value elements added to the mesh get.

        Args:
            value: The new default.
        """
        self._default = value.copy()

    def mesh_tag(self) -> Int:
        """The tag of the mesh the container matches.

        Returns:
            The tag (`MeshTag`).
        """
        return self._tag

    def sync[M: ElementMesh](mut self, mesh: M) raises:
        """Brings the container up to date with its mesh: grows it to the
        mesh's capacity, filling new entries with the default value.

        Parameters:
            M: The mesh type.

        Args:
            mesh: The mesh the container was built on, or a copy of it.

        Raises:
            If `mesh` is a different mesh, or has been compressed since the
            container was built (renumbering arrives with `compress()`).
        """
        var tag = mesh.mesh_tag()
        if tag != self._tag:
            raise Error(
                "MeshData.sync: not the mesh this container was built on,"
                " or the mesh was compressed (not yet supported)"
            )
        var n = mesh.element_capacity(Self.E.KIND)
        while len(self._data) < n:
            self._data.append(self._default.copy())
