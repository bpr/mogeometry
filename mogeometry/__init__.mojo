# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# A port of geometry-central (https://github.com/nmwsharp/geometry-central,
# MIT License, Copyright (c) 2017-2019 Nicholas Sharp and the
# geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
The Top-Level Package (`mogeometry`).
=====================================
Geometry processing on surface meshes and point clouds: halfedge meshes,
their intrinsic and embedded geometry, and the algorithms built on them
(geodesic distance, vector transport, direction fields, intrinsic
triangulations, parameterization, remeshing). A Mojo port of
geometry-central. Sparse linear algebra comes from molina.

Ported so far: the small vectors and elementary geometry of
`mogeometry.utilities`, and `ManifoldSurfaceMesh` with its element
handles, navigation, and data containers (`mogeometry.surface`). `docs/plan.md` gives the order of the rest.
"""

from mogeometry.surface import (
    BoundaryLoop,
    BoundaryLoopData,
    Corner,
    CornerData,
    DiamondBoundary,
    Edge,
    EdgeData,
    EdgeVertices,
    Face,
    FaceData,
    Halfedge,
    HalfedgeData,
    ManifoldSurfaceMesh,
    Vertex,
    VertexData,
)
from mogeometry.utilities import (
    INVALID_IND,
    MESH_CHECKS,
    MESH_INDEX,
    PI,
    DisjointSets,
    ElementHandle,
    ElementKind,
    MarkedDisjointSets,
    MeshData,
    RayRayIntersectionResult2D,
    RaySegmentIntersectionResult2D,
    SegmentSegmentIntersectionResult2D,
    TangentBasis,
    Vector2,
    Vector3,
    clamp,
    fmax,
    fmin,
    in_circle_test,
    layout_triangle_vertex,
    libm_exp,
    libm_pow,
    point_line_segment_distance,
    point_line_segment_nearest_location,
    ray_ray_intersection,
    ray_segment_intersection,
    regularize_angle,
    segment_segment_intersection,
    tet_fourth_point,
    triangle_area,
    triangle_is_valid,
)
