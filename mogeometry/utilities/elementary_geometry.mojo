# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's
# include/geometrycentral/utilities/elementary_geometry.h and .ipp, and
# src/utilities/elementary_geometry.cpp (commit 019669d, MIT License,
# Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
# contributors); see NOTICE. `in_circle_test` evaluates the determinant by
# the cofactor formula of Eigen 3.3's `determinant_impl<Derived, 4>`
# (src/LU/Determinant.h, MPL2, Copyright (C) 2008 Benoit Jacob), which
# geometry-central calls.
# ===----------------------------------------------------------------------=== #
"""
Elementary Geometry (`mogeometry.utilities.elementary_geometry`)
================================================================
Triangles, segments, and rays: triangle areas from edge lengths, laying
out a triangle in the plane, distances to segments, line intersections,
the incircle test, and the fourth vertex of a tetrahedron from its edge
lengths.

Differences from geometry-central:

- `pointLineSegmentNeaestLocation` is spelled
  `point_line_segment_nearest_location`.
- `triTriIntersection` (`src/utilities/tri_tri_intersect.cpp`) is not
  ported yet; it has no users before the mesh-changing algorithms.

None of these functions guards against degenerate input: as in
geometry-central, colinear lines give infinite or NaN parameters, and a
triangle that violates the triangle inequality has area 0.
"""

from std.math import sqrt

from mogeometry.utilities.utilities import clamp
from mogeometry.utilities.vector2 import Vector2
from mogeometry.utilities.vector3 import Vector3


def _std_max(a: Float64, b: Float64) -> Float64:
    # std::max(a, b): b only if a < b, so a NaN b gives a
    return b if a < b else a


# ===----------------------------------------------------------------------=== #
# Two dimensions
# ===----------------------------------------------------------------------=== #


def triangle_area(l_ab: Float64, l_bc: Float64, l_ca: Float64) -> Float64:
    """The area of a triangle from its edge lengths, by Heron's formula.

    Args:
        l_ab: The length of edge AB.
        l_bc: The length of edge BC.
        l_ca: The length of edge CA.

    Returns:
        The area; 0 if the lengths violate the triangle inequality.
    """
    var s = (l_ab + l_bc + l_ca) / 2.0
    var arg = _std_max(0.0, s * (s - l_ab) * (s - l_bc) * (s - l_ca))
    return sqrt(arg)


def triangle_is_valid(
    a: Float64, b: Float64, c: Float64, tolerance: Float64 = 0.0
) -> Bool:
    """Whether three lengths satisfy the triangle inequalities, up to
    `tolerance`.

    Args:
        a: The first length.
        b: The second length.
        c: The third length.
        tolerance: How far each inequality may be violated.

    Returns:
        True if each length is at most the sum of the other two plus
        `tolerance`.
    """
    return (
        a + b >= c - tolerance
        and b + c >= a - tolerance
        and c + a >= b - tolerance
    )


def layout_triangle_vertex(
    p_a: Vector2, p_b: Vector2, l_bc: Float64, l_ca: Float64
) -> Vector2:
    """The third vertex `C` of triangle ABC, given `A`, `B`, and the
    lengths of the other two edges.

    Args:
        p_a: The position of A.
        p_b: The position of B.
        l_bc: The length of edge BC.
        l_ca: The length of edge CA.

    Returns:
        The position of C on the side of AB that winds ABC
        counterclockwise.
    """
    var l_ab = (p_b - p_a).norm()
    var t_area = triangle_area(l_ab, l_bc, l_ca)

    # signed width and height of the right triangle formed by the altitude
    # from C
    var h = 2.0 * t_area / l_ab
    var w = (l_ab * l_ab - l_bc * l_bc + l_ca * l_ca) / (2.0 * l_ab)

    var v_ab_n = (p_b - p_a) / l_ab
    var v_ab_n_perp = Vector2(-v_ab_n.y, v_ab_n.x)
    return p_a + w * v_ab_n + h * v_ab_n_perp


def point_line_segment_distance(
    p: Vector2, line_a: Vector2, line_b: Vector2
) -> Float64:
    """The distance from a point to a segment in the plane.

    Args:
        p: The point.
        line_a: One end of the segment.
        line_b: The other end.

    Returns:
        The distance from `p` to the nearest point of the segment.
    """
    var len2 = (line_a - line_b).norm2()
    if len2 == 0.0:
        return (line_a - p).norm()  # degenerate segment
    var t = clamp((p - line_a).dot(line_b - line_a) / len2, 0.0, 1.0)
    var proj = line_a + t * (line_b - line_a)
    return (p - proj).norm()


@fieldwise_init
struct SegmentSegmentIntersectionResult2D(TrivialRegisterPassable):
    """Where two segments' lines cross, as `segment_segment_intersection`
    returns it."""

    var t_a: Float64
    """The crossing's parameter along the first segment, 0 at its start
    and 1 at its end."""
    var t_b: Float64
    """The crossing's parameter along the second segment."""
    var hit: Bool
    """Whether both parameters are in `[0, 1]`: the segments intersect."""


@fieldwise_init
struct RaySegmentIntersectionResult2D(TrivialRegisterPassable):
    """Where a ray's line crosses a segment's, as `ray_segment_intersection`
    returns it."""

    var t_ray: Float64
    """The crossing's distance along the ray, in units of its direction;
    infinite for parallel lines."""
    var t_line: Float64
    """The crossing's parameter along the segment, 0 at its start and 1
    at its end."""


@fieldwise_init
struct RayRayIntersectionResult2D(TrivialRegisterPassable):
    """Where two rays' lines cross, as `ray_ray_intersection` returns
    it."""

    var t_ray1: Float64
    """The crossing's distance along the first ray, in units of its
    direction; infinite for parallel lines."""
    var t_ray2: Float64
    """The crossing's distance along the second ray."""


def segment_segment_intersection(
    p_a_start: Vector2, p_a_end: Vector2, p_b_start: Vector2, p_b_end: Vector2
) -> SegmentSegmentIntersectionResult2D:
    """Where the lines of two segments cross, and whether the segments do.

    No special care is taken for nearly colinear segments.

    Args:
        p_a_start: The start of the first segment.
        p_a_end: The end of the first segment.
        p_b_start: The start of the second segment.
        p_b_end: The end of the second segment.

    Returns:
        The crossing's parameters along both segments, and whether both
        are in `[0, 1]`.
    """
    var len_a = (p_a_end - p_a_start).norm()
    var ray = ray_segment_intersection(
        p_a_start, (p_a_end - p_a_start) / len_a, p_b_start, p_b_end
    )
    var t_a = ray.t_ray / len_a
    var t_b = ray.t_line
    var hit = t_a >= 0 and t_a <= 1 and t_b >= 0 and t_b <= 1
    return SegmentSegmentIntersectionResult2D(t_a, t_b, hit)


def ray_segment_intersection(
    ray_start: Vector2, ray_dir: Vector2, line_start: Vector2, line_end: Vector2
) -> RaySegmentIntersectionResult2D:
    """Where a ray's line crosses a segment's line.

    The parameters are not limited to the ray or the segment; check
    `t_ray >= 0` and `0 <= t_line <= 1` for an actual hit. No special
    care is taken for nearly parallel lines.

    Args:
        ray_start: The ray's origin.
        ray_dir: The ray's direction.
        line_start: The start of the segment.
        line_end: The end of the segment.

    Returns:
        The crossing's parameters along the ray and the segment.
    """
    var v1 = ray_start - line_start
    var v2 = line_end - line_start
    var v3 = Vector2(-ray_dir.y, ray_dir.x)
    var cross21 = v2.cross(v1)
    var t_ray = cross21 / v2.dot(v3)
    var t_line = v1.dot(v3) / v2.dot(v3)
    return RaySegmentIntersectionResult2D(t_ray, t_line)


def ray_ray_intersection(
    ray1_start: Vector2,
    ray1_dir: Vector2,
    ray2_start: Vector2,
    ray2_dir: Vector2,
) -> RayRayIntersectionResult2D:
    """Where the lines of two rays cross.

    The parameters may be negative; no special care is taken for nearly
    parallel rays.

    Args:
        ray1_start: The first ray's origin.
        ray1_dir: The first ray's direction.
        ray2_start: The second ray's origin.
        ray2_dir: The second ray's direction.

    Returns:
        The crossing's parameters along both rays.
    """
    var v1 = ray1_start - ray2_start
    var v2 = ray2_dir
    var v3 = Vector2(-ray1_dir.y, ray1_dir.x)
    var cross21 = v2.cross(v1)
    var t_ray1 = cross21 / v2.dot(v3)
    var t_ray2 = v1.dot(v3) / v2.dot(v3)
    return RayRayIntersectionResult2D(t_ray1, t_ray2)


def _det4_helper(
    m: SIMD[DType.float64, 16], j: Int, k: Int, r: Int, s: Int
) -> Float64:
    # Eigen's bruteforce_det4_helper, m row major
    return (m[4 * j + 0] * m[4 * k + 1] - m[4 * k + 0] * m[4 * j + 1]) * (
        m[4 * r + 2] * m[4 * s + 3] - m[4 * s + 2] * m[4 * r + 3]
    )


def in_circle_test(
    p_a: Vector2, p_b: Vector2, p_c: Vector2, p_test: Vector2
) -> Bool:
    """Whether `p_test` is inside the circumcircle of the counterclockwise
    triangle ABC.

    A plain floating-point determinant: no exact predicates or epsilons,
    so the answer near the circle is not reliable.

    Args:
        p_a: The first vertex of the triangle.
        p_b: The second vertex.
        p_c: The third vertex.
        p_test: The point to test.

    Returns:
        True if the incircle determinant is positive.
    """
    # fmt: off
    var m = SIMD[DType.float64, 16](
        p_a.x, p_a.y, p_a.norm2(), 1.0,
        p_b.x, p_b.y, p_b.norm2(), 1.0,
        p_c.x, p_c.y, p_c.norm2(), 1.0,
        p_test.x, p_test.y, p_test.norm2(), 1.0,
    )
    # fmt: on
    var det = (
        _det4_helper(m, 0, 1, 2, 3)
        - _det4_helper(m, 0, 2, 1, 3)
        + _det4_helper(m, 0, 3, 1, 2)
        + _det4_helper(m, 1, 2, 0, 3)
        - _det4_helper(m, 1, 3, 0, 2)
        + _det4_helper(m, 2, 3, 0, 1)
    )
    return det > 0.0


# ===----------------------------------------------------------------------=== #
# Three dimensions
# ===----------------------------------------------------------------------=== #


def point_line_segment_nearest_location(
    p: Vector3, line_a: Vector3, line_b: Vector3
) -> Float64:
    """The parameter of the point of a segment nearest to `p`.

    Args:
        p: The point.
        line_a: One end of the segment (parameter 0).
        line_b: The other end (parameter 1).

    Returns:
        The parameter in `[0, 1]`; 0 for a degenerate segment.
    """
    var len2 = (line_a - line_b).norm2()
    if len2 == 0.0:
        return 0.0  # degenerate segment
    return clamp((p - line_a).dot(line_b - line_a) / len2, 0.0, 1.0)


def point_line_segment_distance(
    p: Vector3, line_a: Vector3, line_b: Vector3
) -> Float64:
    """The distance from a point to a segment in space.

    Args:
        p: The point.
        line_a: One end of the segment.
        line_b: The other end.

    Returns:
        The distance from `p` to the nearest point of the segment.
    """
    var t = point_line_segment_nearest_location(p, line_a, line_b)
    var proj = line_a + t * (line_b - line_a)
    return (p - proj).norm()


def tet_fourth_point(
    a: Vector3,
    b: Vector3,
    c: Vector3,
    l_ad: Float64,
    l_bd: Float64,
    l_cd: Float64,
) -> Vector3:
    """The fourth vertex `d` of a tetrahedron, given the other three and
    the lengths of the edges to it.

    Args:
        a: The first vertex.
        b: The second vertex.
        c: The third vertex.
        l_ad: The length of edge ad.
        l_bd: The length of edge bd.
        l_cd: The length of edge cd.

    Returns:
        The position of `d` that gives the tetrahedron `(a, b, c, d)`
        positive orientation. If the lengths admit no tetrahedron, `d` is
        placed in the plane of `abc`.
    """
    # a local orthonormal frame for the base triangle, based at a
    var e1 = (b - a) / (b - a).norm()
    var e2 = c - a
    e2 -= e2.dot(e1) * e1
    e2 = e2 / e2.norm()
    var e3 = e1.cross(e2)

    # in that frame a, b, c are (0, 0), (x2, 0), (x3, y3)
    var x2 = (b - a).dot(e1)
    var x3 = (c - a).dot(e1)
    var y3 = (c - a).dot(e2)

    # solve the three quadratic equations for d's coordinates
    var x4 = (l_ad * l_ad - l_bd * l_bd + x2 * x2) / (2.0 * x2)
    var y4 = (l_ad * l_ad - l_cd * l_cd + x3 * x3 + y3 * y3 - 2.0 * x4 * x3) / (
        2.0 * y3
    )
    var z4 = sqrt(_std_max(0.0, l_ad * l_ad - x4 * x4 - y4 * y4))

    # back to world coordinates: the two solutions, reflected across abc
    var d1 = a + x4 * e1 + y4 * e2 + z4 * e3
    var d2 = a + x4 * e1 + y4 * e2 - z4 * e3

    # return the positively oriented one
    var vol1 = (d1 - a).dot((b - a).cross(c - a))
    if vol1 < 0.0:
        return d1
    return d2
