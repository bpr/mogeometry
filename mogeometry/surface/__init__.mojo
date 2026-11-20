"""
Surfaces (`mogeometry.surface`)
===============================
Halfedge meshes, the data stored on their elements, their geometry, and
the algorithms on them: geometry-central's `surface/`.
"""

from mogeometry.surface.halfedge_element_types import (
    BoundaryLoop,
    BoundaryLoopData,
    Corner,
    CornerData,
    DiamondBoundary,
    Edge,
    EdgeData,
    EdgeVertices,
    ElementRange,
    Face,
    FaceData,
    Halfedge,
    HalfedgeData,
    NavigationIterator,
    Navigator,
    RangeFilter,
    Vertex,
    VertexData,
)
from mogeometry.surface.manifold_surface_mesh import ManifoldSurfaceMesh
