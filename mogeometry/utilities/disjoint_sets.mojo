# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's include/geometrycentral/utilities/
# disjoint_sets.h and src/utilities/disjoint_sets.cpp (commit 019669d, MIT
# License, Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
# contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Disjoint sets (`mogeometry.utilities.disjoint_sets`)
====================================================
Union-find over the integers `0 .. n`, with union by rank and path
compression: `DisjointSets`, and `MarkedDisjointSets`, which also tracks
whether each set is marked.

Differences from geometry-central: `find` is a loop rather than a
recursion (the same path compression), so long chains cannot overflow
the stack.
"""


struct DisjointSets(Copyable, Movable):
    """Disjoint sets of the integers `0 .. n` (inclusive, as
    geometry-central's), each starting alone."""

    var _parent: List[Int]
    """Each element's parent; a root is its own parent."""
    var _rank: List[Int]
    """An upper bound on the height of each root's tree."""

    def __init__(out self, n: Int):
        """Sets `{0}, {1}, ..., {n}`.

        Args:
            n: The largest element.
        """
        self._parent = List[Int](capacity=n + 1)
        for i in range(n + 1):
            self._parent.append(i)
        self._rank = List[Int](length=n + 1, fill=0)

    def find(mut self, x: Int) -> Int:
        """The representative of the set containing `x`; compresses the
        path from `x` to it.

        Args:
            x: An element.

        Returns:
            The root of `x`'s set.
        """
        var root = x
        while self._parent[root] != root:
            root = self._parent[root]
        var y = x
        while self._parent[y] != root:
            var next = self._parent[y]
            self._parent[y] = root
            y = next
        return root

    def merge(mut self, x: Int, y: Int):
        """Unites the sets containing `x` and `y`, by rank.

        Args:
            x: An element of one set.
            y: An element of the other.
        """
        var rx = self.find(x)
        var ry = self.find(y)
        if self._rank[rx] > self._rank[ry]:
            self._parent[ry] = rx
        else:
            self._parent[rx] = ry
        if self._rank[rx] == self._rank[ry]:
            self._rank[ry] += 1


struct MarkedDisjointSets(Copyable, Movable):
    """Disjoint sets of the integers `0 .. n` that can be marked; a union
    is marked if either part was."""

    var _sets: DisjointSets
    """The sets."""
    var _marked: List[Bool]
    """Whether each root's set is marked."""

    def __init__(out self, n: Int):
        """Sets `{0}, {1}, ..., {n}`, none marked.

        Args:
            n: The largest element.
        """
        self._sets = DisjointSets(n)
        self._marked = List[Bool](length=n + 1, fill=False)

    def find(mut self, x: Int) -> Int:
        """The representative of the set containing `x`.

        Args:
            x: An element.

        Returns:
            The root of `x`'s set.
        """
        return self._sets.find(x)

    def merge(mut self, x: Int, y: Int):
        """Unites the sets containing `x` and `y`; the union is marked if
        either was.

        Args:
            x: An element of one set.
            y: An element of the other.
        """
        var rx = self._sets.find(x)
        var ry = self._sets.find(y)
        self._sets.merge(rx, ry)
        if self._marked[rx] or self._marked[ry]:
            self._marked[rx] = True
            self._marked[ry] = True

    def mark(mut self, x: Int):
        """Marks the set containing `x`.

        Args:
            x: An element.
        """
        self._marked[self._sets.find(x)] = True

    def unmark(mut self, x: Int):
        """Unmarks the set containing `x`.

        Args:
            x: An element.
        """
        self._marked[self._sets.find(x)] = False

    def is_marked(mut self, x: Int) -> Bool:
        """Whether the set containing `x` is marked.

        Args:
            x: An element.

        Returns:
            True if the set is marked.
        """
        return self._marked[self._sets.find(x)]
