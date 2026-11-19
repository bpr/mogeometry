# Vector2, Vector3, elementary geometry, and the scalar helpers, against
# geometry-central: tests/data/utilities/reference.txt is written by
# tests/data/utilities/reference.cpp (`pixi run -e ref ref-data`), which
# calls geometry-central on random inputs. Hand-checked cases cover what
# random inputs do not reach (special values, degenerate input).

from std.math import isnan, sqrt
from std.memory import bitcast
from std.testing import (
    TestSuite,
    assert_almost_equal,
    assert_equal,
    assert_false,
    assert_true,
)

from mogeometry import (
    PI,
    fmax,
    fmin,
    libm_exp,
    libm_pow,
    Vector2,
    Vector3,
    clamp,
    in_circle_test,
    layout_triangle_vertex,
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

comptime REFERENCE = "tests/data/utilities/reference.txt"


# ===----------------------------------------------------------------------=== #
# Reading the reference data
# ===----------------------------------------------------------------------=== #


@fieldwise_init
struct Case(Copyable, Movable):
    """One line of reference data: inputs and geometry-central's outputs."""

    var name: String
    var inputs: List[Float64]
    var outputs: List[Float64]


def _bits_to_float(word: StringSlice) raises -> Float64:
    return bitcast[DType.float64](Int64(atol(word)))


def _cases(name: String) raises -> List[Case]:
    var out = List[Case]()
    with open(REFERENCE, "r") as f:
        for line in f.read().split("\n"):
            var words = line.split()
            if len(words) == 0 or words[0] != name:
                continue
            var inputs = List[Float64]()
            var outputs = List[Float64]()
            var seen_arrow = False
            for i in range(1, len(words)):
                if words[i] == "->":
                    seen_arrow = True
                elif seen_arrow:
                    outputs.append(_bits_to_float(words[i]))
                else:
                    inputs.append(_bits_to_float(words[i]))
            out.append(Case(String(name), inputs^, outputs^))
    if len(out) == 0:
        raise Error("no reference cases named " + name)
    return out^


struct Inputs:
    """Reads a case's inputs in order."""

    var values: List[Float64]
    var pos: Int

    def __init__(out self, values: List[Float64]):
        self.values = values.copy()
        self.pos = 0

    def f(mut self) -> Float64:
        self.pos += 1
        return self.values[self.pos - 1]

    def v2(mut self) -> Vector2:
        var x = self.f()
        return Vector2(x, self.f())

    def v3(mut self) -> Vector3:
        var x = self.f()
        var y = self.f()
        return Vector3(x, y, self.f())


struct Outputs:
    """Collects the port's outputs in the order the reference lists them."""

    var values: List[Float64]

    def __init__(out self):
        self.values = List[Float64]()

    def f(mut self, x: Float64):
        self.values.append(x)

    def b(mut self, x: Bool):
        self.values.append(1.0 if x else 0.0)

    def v2(mut self, v: Vector2):
        self.f(v.x)
        self.f(v.y)

    def v3(mut self, v: Vector3):
        self.f(v.x)
        self.f(v.y)
        self.f(v.z)


comptime EPS = 2.220446049250313e-16
"""Double precision machine epsilon, `2^-52`."""


def _compare(cases: List[Case], got: List[List[Float64]], ulps: Float64) raises:
    # every output within `ulps` roundoffs of geometry-central's, a
    # roundoff being EPS * max(1, |value|), and NaN matching NaN. Mojo
    # contracts a * b + c into fused multiply-adds where the reference
    # (built with -ffp-contract=off) rounds twice, so results differ in
    # the last bits, by more where the problem is ill-conditioned. The
    # worst case is found first and printed, so a pass shows the margin.
    # The driver keeps its random inputs well-conditioned, so the observed
    # worst is under 12 roundoffs; ill-conditioned input is tested by hand.
    var worst = 0.0
    var worst_k = 0
    var worst_i = 0
    for k in range(len(cases)):
        ref want = cases[k].outputs
        assert_equal(len(got[k]), len(want), cases[k].name)
        for i in range(len(want)):
            if isnan(want[i]) or isnan(got[k][i]):
                assert_true(
                    isnan(want[i]) and isnan(got[k][i]),
                    String(cases[k].name, " case ", k, " output ", i, ": NaN"),
                )
                continue
            if got[k][i] == want[i]:
                continue
            var err = abs(got[k][i] - want[i]) / (EPS * max(1.0, abs(want[i])))
            if err > worst:
                worst = err
                worst_k = k
                worst_i = i
    print("   ", cases[0].name, "worst difference in roundoffs:", worst)
    ref want = cases[worst_k].outputs[worst_i]
    assert_almost_equal(
        got[worst_k][worst_i],
        want,
        atol=ulps * EPS * max(1.0, abs(want)),
        rtol=0,
        msg=String(cases[worst_k].name, " case ", worst_k, " output ", worst_i),
    )


# ===----------------------------------------------------------------------=== #
# Against geometry-central
# ===----------------------------------------------------------------------=== #


def test_vector3_matches_reference() raises:
    var cases = _cases("vector3")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var u = r.v3()
        var v = r.v3()
        var w = Outputs()
        w.f(u.norm())
        w.f(u.norm2())
        w.v3(u.normalize())
        w.v3(u.normalize_cutoff(1.5))
        w.v3(u.cross(v))
        w.f(u.dot(v))
        w.f(u.sum())
        w.f(u.angle(v))
        w.v3(u / 3.0)
        w.v3(3.0 * u - v)
        var q = u
        q /= 3.0
        w.v3(q)
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_vector3_frames_match_reference() raises:
    var cases = _cases("vector3_frame")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var u = r.v3()
        var v = r.v3()
        var n = r.v3()
        var theta = r.f()
        var n_hat = r.v3()
        var w = Outputs()
        w.f(u.angle_in_plane(v, n))
        w.v3(u.rotate_around(n, theta))
        w.v3(u.remove_component(n_hat))
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_tangent_basis_matches_reference() raises:
    var cases = _cases("vector3_tangent_basis")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var basis = r.v3().build_tangent_basis()
        var w = Outputs()
        w.v3(basis.basis_x)
        w.v3(basis.basis_y)
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_vector2_matches_reference() raises:
    var cases = _cases("vector2")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var u = r.v2()
        var v = r.v2()
        var theta = r.f()
        var w = Outputs()
        w.v2(u * v)
        w.v2(u / v)
        w.v2(u.inv())
        w.f(u.norm())
        w.v2(u.normalize())
        w.v2(u.normalize_cutoff(1.5))
        w.f(u.arg())
        w.f(u.angle(v))
        w.f(u.oriented_angle(v))
        w.f(u.dot(v))
        w.f(u.cross(v))
        w.v2(u.rotate_cw(theta))
        w.v2(u.rotate_ccw(theta))
        w.v2(Vector2.from_angle(theta))
        w.v2(u / 3.0)
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_vector2_pow_matches_reference() raises:
    var cases = _cases("vector2_pow")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var z = r.v2()
        var p = r.f()
        var q = r.v2()
        var w = Outputs()
        w.v2(z.pow(p))
        w.v2(z.pow(q))
        got.append(w.values.copy())
    # exp(e) multiplies the error in e by |e| (up to ~9 here)
    _compare(cases, got, ulps=32)


def test_regularize_angle_matches_reference() raises:
    var cases = _cases("regularize_angle")
    var got = List[List[Float64]]()
    for c in cases:
        var w = Outputs()
        w.f(regularize_angle(c.inputs[0]))
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_triangles_match_reference() raises:
    var cases = _cases("triangle")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var p_a = r.v2()
        var p_b = r.v2()
        var l_ab = r.f()
        var l_bc = r.f()
        var l_ca = r.f()
        var w = Outputs()
        w.f(triangle_area(l_ab, l_bc, l_ca))
        w.b(triangle_is_valid(l_ab, l_bc, l_ca))
        w.b(triangle_is_valid(l_ab, l_bc, l_ca, 0.1))
        w.v2(layout_triangle_vertex(p_a, p_b, l_bc, l_ca))
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_segments_match_reference() raises:
    var cases = _cases("segments2")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var p = r.v2()
        var a = r.v2()
        var b = r.v2()
        var d = r.v2()
        var s = segment_segment_intersection(p, a, b, d)
        var ray = ray_segment_intersection(p, a, b, d)
        var rr = ray_ray_intersection(p, a, b, d)
        var w = Outputs()
        w.f(point_line_segment_distance(p, a, b))
        w.f(s.t_a)
        w.f(s.t_b)
        w.b(s.hit)
        w.f(ray.t_ray)
        w.f(ray.t_line)
        w.f(rr.t_ray1)
        w.f(rr.t_ray2)
        w.b(in_circle_test(a, b, d, p))
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


def test_space_matches_reference() raises:
    var cases = _cases("space3")
    var got = List[List[Float64]]()
    for c in cases:
        var r = Inputs(c.inputs)
        var p = r.v3()
        var a = r.v3()
        var b = r.v3()
        var d = r.v3()
        var l_ad = r.f()
        var l_bd = r.f()
        var l_cd = r.f()
        var w = Outputs()
        w.f(point_line_segment_nearest_location(p, a, b))
        w.f(point_line_segment_distance(p, a, b))
        w.v3(tet_fourth_point(a, b, d, l_ad, l_bd, l_cd))
        got.append(w.values.copy())
    _compare(cases, got, ulps=16)


# ===----------------------------------------------------------------------=== #
# Special values and degenerate input
# ===----------------------------------------------------------------------=== #


def test_vector3_basics() raises:
    var v = Vector3(1, 2, 3)
    assert_equal(v[0], 1)
    assert_equal(v[2], 3)
    v[1] = 5
    assert_true(v == Vector3(1, 5, 3))
    assert_true(v != Vector3(1, 2, 3))
    assert_true(Vector3.zero() == Vector3.constant(0))
    assert_true(-v + v == Vector3.zero())
    assert_equal(String(Vector3(1, -2, 0.5)), "<1.0, -2.0, 0.5>")
    assert_true(v.is_finite() and v.is_defined())
    assert_false(Vector3.infinity().is_finite())
    assert_true(Vector3.infinity().is_defined())
    assert_false(Vector3.undefined().is_defined())
    var w = Vector3(1, 2, 3)
    w += Vector3(1, 1, 1)
    w -= Vector3(0, 1, 0)
    w *= 2.0
    assert_true(w == Vector3(4, 4, 8))


def test_vector3_normalize_cutoff_keeps_short_vectors() raises:
    var short = Vector3(0.1, 0.2, 0.2)  # length 0.3
    assert_true(short.normalize_cutoff(0.5) == short)
    assert_almost_equal(
        short.normalize_cutoff(0.1).norm(), 1, atol=4 * EPS, rtol=0
    )
    # the zero vector: NaNs from normalize, itself from normalize_cutoff
    assert_false(Vector3.zero().normalize().is_defined())
    assert_true(Vector3.zero().normalize_cutoff() == Vector3.zero())


def test_vector3_componentwise() raises:
    var nan = Vector3.undefined().x
    var u = Vector3(1, nan, -3)
    var v = Vector3(2, 0, -4)
    assert_true(u.componentwise_min(v) == Vector3(1, 0, -4))
    assert_true(u.componentwise_max(v) == Vector3(2, 0, -3))
    var c = Vector3(-5, 0.5, 9).clamp(Vector3.constant(0), Vector3.constant(1))
    assert_true(c == Vector3(0, 0.5, 1))
    assert_equal(Vector3(1, 2, 3).sum(), 6)


def test_tangent_basis_is_orthonormal() raises:
    # both branches: along x (test vector y) and elsewhere (test vector x)
    var normals: List[Vector3] = [Vector3(3, 0.1, -0.2), Vector3(0.2, -1, 2)]
    for n in normals:
        var b = n.build_tangent_basis()
        var n_hat = n.unit()
        assert_almost_equal(b.basis_x.norm(), 1, atol=2 * EPS, rtol=0)
        assert_almost_equal(b.basis_y.norm(), 1, atol=2 * EPS, rtol=0)
        assert_almost_equal(b.basis_x.dot(b.basis_y), 0, atol=2 * EPS, rtol=0)
        assert_almost_equal(b.basis_x.dot(n_hat), 0, atol=2 * EPS, rtol=0)
        var z = b.basis_x.cross(b.basis_y)
        for i in range(3):
            assert_almost_equal(z[i], n_hat[i], atol=4 * EPS, rtol=0)


def test_rotate_around() raises:
    # a quarter turn about z takes x to y; a vector along the axis stays
    var r = Vector3(1, 0, 0).rotate_around(Vector3(0, 0, 2), PI / 2)
    assert_almost_equal(r.x, 0, atol=2 * EPS, rtol=0)
    assert_almost_equal(r.y, 1, atol=2 * EPS, rtol=0)
    assert_equal(r.z, 0)
    assert_true(
        Vector3(0, 0, 3).rotate_around(Vector3(0, 0, 1), 1.0)
        == Vector3(0, 0, 3)
    )


def test_vector2_complex_arithmetic() raises:
    var i = Vector2(0, 1)
    assert_true(i * i == Vector2(-1, 0))
    assert_true(Vector2(1, 2) * Vector2(3, 4) == Vector2(-5, 10))
    assert_true(Vector2(-5, 10) / Vector2(3, 4) == Vector2(1, 2))
    assert_true(Vector2(0, 2).inv() == Vector2(0, -0.5))
    assert_true(Vector2(1, 2).conj() == Vector2(1, -2))
    assert_true(Vector2(1, 2).rotate90() == i * Vector2(1, 2))
    var z = Vector2(1, 2)
    z *= Vector2(3, 4)
    assert_true(z == Vector2(-5, 10))
    z /= Vector2(3, 4)
    assert_true(z == Vector2(1, 2))
    z *= 2.0
    assert_true(z == Vector2(2, 4))
    assert_true(
        Vector2.from_complex(Vector2(3, -1).to_complex()) == Vector2(3, -1)
    )
    assert_true(Vector2(1, 0).cross3(Vector2(0, 1)) == Vector3(0, 0, 1))
    assert_equal(String(Vector2(1, -2)), "<1.0, -2.0>")


def test_vector2_pow_branches() raises:
    # a positive real base takes the real pow: exact for exact results
    assert_true(Vector2(4, 0).pow(0.5) == Vector2(2, 0))
    assert_true(Vector2(2, 0).pow(-2.0) == Vector2(0.25, 0))
    # the principal square root of -1 is i
    var r = Vector2(-1, 0).pow(0.5)
    assert_almost_equal(r.x, 0, atol=EPS, rtol=0)
    assert_almost_equal(r.y, 1, atol=EPS, rtol=0)
    # i^i = exp(-pi/2), real
    var ii = Vector2(0, 1).pow(Vector2(0, 1))
    assert_almost_equal(ii.x, 0.20787957635076193, atol=2 * EPS, rtol=0)
    assert_almost_equal(ii.y, 0, atol=EPS, rtol=0)


def test_vector2_angles() raises:
    assert_almost_equal(
        Vector2(1, 0).oriented_angle(Vector2(0, 1)), PI / 2, atol=EPS, rtol=0
    )
    assert_almost_equal(
        Vector2(0, 1).oriented_angle(Vector2(1, 0)), -PI / 2, atol=EPS, rtol=0
    )
    assert_almost_equal(
        Vector2(1, 0).angle(Vector2(-1, 0)), PI, atol=EPS, rtol=0
    )
    # parallel vectors: the cosine of unit vectors rounds to within an ulp
    # of 1, and acos's slope there is infinite, so the angle is only good
    # to sqrt(2 eps) ~ 2e-8; the clamp keeps a cosine past 1 from giving NaN
    var pairs: List[Vector3] = [
        Vector3(0.1, 0.3, 0.7),
        Vector3(1, 1, 0),
        Vector3(3, -7, 2),
        Vector3(0.3, 0.3, 0.3),
    ]
    for v in pairs:
        var a3 = v.angle(v * 3.0)
        assert_false(isnan(a3))
        assert_almost_equal(a3, 0, atol=sqrt(4 * EPS), rtol=0)
        var w = Vector2(v.x, v.y)
        var a2 = w.angle(w * 3.0)
        assert_false(isnan(a2))
        assert_almost_equal(a2, 0, atol=sqrt(4 * EPS), rtol=0)
    var r = Vector2(1, 0).rotate_ccw(PI / 2)
    assert_almost_equal(r.x, 0, atol=EPS, rtol=0)
    assert_almost_equal(r.y, 1, atol=EPS, rtol=0)


def test_scalar_helpers() raises:
    var nan = Vector3.undefined().x
    assert_equal(clamp(5, 0, 1), 1)
    assert_equal(clamp(-5, 0, 1), 0)
    assert_true(isnan(clamp(nan, 0, 1)))
    assert_equal(regularize_angle(2 * PI), 0)
    assert_almost_equal(
        regularize_angle(-PI / 2), 3 * PI / 2, atol=4 * EPS, rtol=0
    )
    # libm's pow is correctly rounded here; Mojo's ** is not
    assert_equal(libm_pow(2, 0.5), sqrt(2.0))
    assert_equal(libm_exp(0), 1)
    assert_equal(fmin(nan, 2), 2)
    assert_equal(fmax(1, nan), 1)


def test_triangle_degenerate() raises:
    assert_equal(triangle_area(3, 4, 5), 6)
    # lengths that violate the triangle inequality give area 0
    assert_equal(triangle_area(1, 1, 3), 0)
    assert_false(triangle_is_valid(1, 1, 3))
    assert_true(triangle_is_valid(1, 1, 2.05, tolerance=0.1))
    # the 3-4-5 triangle, counterclockwise
    var c = layout_triangle_vertex(Vector2(0, 0), Vector2(3, 0), 4, 5)
    assert_true(c == Vector2(3, 4))


def test_segments_degenerate() raises:
    # a degenerate segment is a point
    assert_equal(
        point_line_segment_distance(
            Vector2(3, 4), Vector2(0, 0), Vector2(0, 0)
        ),
        5,
    )
    assert_equal(
        point_line_segment_nearest_location(
            Vector3(1, 1, 1), Vector3(0, 0, 0), Vector3(0, 0, 0)
        ),
        0,
    )
    # beyond an end, the distance is to that end
    assert_equal(
        point_line_segment_distance(
            Vector3(5, 0, 0), Vector3(0, 0, 0), Vector3(1, 0, 0)
        ),
        4,
    )
    # crossing diagonals of the unit square
    var s = segment_segment_intersection(
        Vector2(0, 0), Vector2(1, 1), Vector2(0, 1), Vector2(1, 0)
    )
    assert_equal(s.t_a, 0.5)
    assert_equal(s.t_b, 0.5)
    assert_true(s.hit)
    var miss = segment_segment_intersection(
        Vector2(0, 0), Vector2(1, 1), Vector2(3, 0), Vector2(3, 1)
    )
    assert_false(miss.hit)
    # parallel lines: no crossing, an infinite or NaN parameter, no hit
    var par = ray_segment_intersection(
        Vector2(0, 0), Vector2(1, 0), Vector2(0, 1), Vector2(1, 1)
    )
    assert_false(Vector2(par.t_ray, 0).is_finite())
    var rr = ray_ray_intersection(
        Vector2(0, 0), Vector2(1, 0), Vector2(0, 1), Vector2(2, 0)
    )
    assert_false(Vector2(rr.t_ray1, rr.t_ray2).is_finite())


def test_in_circle() raises:
    var a = Vector2(1, 0)
    var b = Vector2(0, 1)
    var c = Vector2(-1, 0)
    assert_true(in_circle_test(a, b, c, Vector2(0.1, -0.2)))
    assert_false(in_circle_test(a, b, c, Vector2(0, -1.5)))
    # clockwise triangles flip the sign
    assert_false(in_circle_test(c, b, a, Vector2(0.1, -0.2)))


def test_tet_fourth_point() raises:
    # the corner of the unit cube's corner tetrahedron, on geometry-central's
    # "positive" side, which is -z for this base
    var d = tet_fourth_point(
        Vector3(0, 0, 0),
        Vector3(1, 0, 0),
        Vector3(0, 1, 0),
        1,
        sqrt(2.0),
        sqrt(2.0),
    )
    assert_almost_equal(d.x, 0, atol=2 * EPS, rtol=0)
    assert_almost_equal(d.y, 0, atol=2 * EPS, rtol=0)
    assert_almost_equal(d.z, -1, atol=2 * EPS, rtol=0)
    # lengths that admit no tetrahedron put d in the base's plane
    var flat = tet_fourth_point(
        Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0), 0.1, 2, 2
    )
    assert_equal(flat.z, 0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
