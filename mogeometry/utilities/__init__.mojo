"""
Utilities (`mogeometry.utilities`)
==================================
Small vectors (`Vector2`, `Vector3`), quaternions, elementary geometry,
disjoint sets, nearest neighbors, and the element handles and data
containers that meshes share: geometry-central's `utilities/`.
"""

from mogeometry.utilities.disjoint_sets import DisjointSets, MarkedDisjointSets
from mogeometry.utilities.element import (
    INVALID_IND,
    MESH_CHECKS,
    MESH_INDEX,
    ElementHandle,
    ElementKind,
    ElementMesh,
    MeshElement,
    MeshTag,
)
from mogeometry.utilities.elementary_geometry import (
    RayRayIntersectionResult2D,
    RaySegmentIntersectionResult2D,
    SegmentSegmentIntersectionResult2D,
    in_circle_test,
    layout_triangle_vertex,
    point_line_segment_distance,
    point_line_segment_nearest_location,
    ray_ray_intersection,
    ray_segment_intersection,
    segment_segment_intersection,
    tet_fourth_point,
    triangle_area,
    triangle_is_valid,
)
from mogeometry.utilities.mesh_data import MeshData
from mogeometry.utilities.utilities import (
    PI,
    clamp,
    fmax,
    fmin,
    libm_exp,
    libm_pow,
    regularize_angle,
)
from mogeometry.utilities.vector2 import Vector2
from mogeometry.utilities.vector3 import TangentBasis, Vector3
