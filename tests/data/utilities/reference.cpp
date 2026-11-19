// Reference values for tests/test_utilities.mojo: geometry-central's
// Vector2, Vector3, elementary geometry, and regularizeAngle on random
// inputs. Build and run with `pixi run -e ref ref-data`.
//
// Output: one case per line, `name in... -> out...`, every number a double
// written as its IEEE bit pattern (a signed 64-bit integer), so it reads
// back exactly. Booleans are written as the doubles 0 and 1.

#include "geometrycentral/utilities/elementary_geometry.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace geometrycentral;

static std::mt19937_64 rng(20261118);
static double uniform(double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(rng); }
static Vector2 rand2(double r = 2.) { return Vector2{uniform(-r, r), uniform(-r, r)}; }
static Vector3 rand3(double r = 2.) { return Vector3{uniform(-r, r), uniform(-r, r), uniform(-r, r)}; }

struct Line {
  std::vector<double> in, out;
  void i(double x) { in.push_back(x); }
  void i(Vector2 v) { i(v.x), i(v.y); }
  void i(Vector3 v) { i(v.x), i(v.y), i(v.z); }
  void o(double x) { out.push_back(x); }
  void o(Vector2 v) { o(v.x), o(v.y); }
  void o(Vector3 v) { o(v.x), o(v.y), o(v.z); }
};

static void put(double x) {
  int64_t b;
  std::memcpy(&b, &x, 8);
  std::printf(" %lld", (long long)b);
}

static void emit(const char* name, const Line& l) {
  std::printf("%s", name);
  for (double x : l.in) put(x);
  std::printf(" ->");
  for (double x : l.out) put(x);
  std::printf("\n");
}

const int N = 40;

int main() {
  // === Vector3

  for (int k = 0; k < N; k++) {
    Line l;
    Vector3 u = rand3(), v = rand3();
    l.i(u), l.i(v);
    l.o(norm(u)), l.o(norm2(u)), l.o(normalize(u)), l.o(normalizeCutoff(u, 1.5));
    l.o(cross(u, v)), l.o(dot(u, v)), l.o(sum(u)), l.o(angle(u, v));
    l.o(u / 3.), l.o(3. * u - v);
    Vector3 w = u;
    w /= 3.;
    l.o(w);
    emit("vector3", l);
  }
  for (int k = 0; k < N; k++) {
    Line l;
    Vector3 u = rand3(), v = rand3(), n = rand3();
    double theta = uniform(-7., 7.);
    Vector3 nHat = unit(n);
    l.i(u), l.i(v), l.i(n), l.i(theta), l.i(nHat);
    l.o(angleInPlane(u, v, n)), l.o(u.rotateAround(n, theta)), l.o(u.removeComponent(nHat));
    emit("vector3_frame", l);
  }
  for (int k = 0; k < N; k++) {
    Line l;
    // half nearly along x, to take buildTangentBasis's other branch
    Vector3 u = k % 2 ? rand3() : Vector3{uniform(1., 2.), uniform(-.2, .2), uniform(-.2, .2)};
    l.i(u);
    std::array<Vector3, 2> b = u.buildTangentBasis();
    l.o(b[0]), l.o(b[1]);
    emit("vector3_tangent_basis", l);
  }

  // === Vector2

  for (int k = 0; k < N; k++) {
    Line l;
    Vector2 u = rand2(), v = rand2();
    double theta = uniform(-7., 7.);
    l.i(u), l.i(v), l.i(theta);
    l.o(u * v), l.o(u / v), l.o(u.inv()), l.o(norm(u)), l.o(normalize(u)), l.o(u.normalizeCutoff(1.5));
    l.o(arg(u)), l.o(angle(u, v)), l.o(orientedAngle(u, v)), l.o(dot(u, v)), l.o(cross(u, v));
    l.o(u.rotateCW(theta)), l.o(u.rotateCCW(theta)), l.o(Vector2::fromAngle(theta)), l.o(u / 3.);
    emit("vector2", l);
  }
  for (int k = 0; k < N; k++) {
    Line l;
    Vector2 z = rand2();
    if (k % 4 == 0) z.y = 0.;  // the real axis: positive takes std::pow's real branch
    double p = uniform(-3., 3.);
    Vector2 q = rand2(1.5);
    l.i(z), l.i(p), l.i(q);
    l.o(z.pow(p)), l.o(z.pow(q));
    emit("vector2_pow", l);
  }

  // === utilities.h

  for (int k = 0; k < N; k++) {
    Line l;
    double theta = uniform(-20., 20.);
    l.i(theta);
    l.o(regularizeAngle(theta));
    emit("regularize_angle", l);
  }

  // === Elementary geometry

  for (int k = 0; k < N; k++) {
    Line l;
    Vector2 pA = rand2(), pB = rand2(), pC = rand2();
    double lAB = norm(pB - pA), lBC = norm(pC - pB), lCA = norm(pA - pC);
    if (k % 5 == 0) lCA = lAB + lBC + uniform(0., .3);  // violates the triangle inequality
    l.i(pA), l.i(pB), l.i(lAB), l.i(lBC), l.i(lCA);
    l.o(triangleArea(lAB, lBC, lCA)), l.o(triangleIsValid(lAB, lBC, lCA)), l.o(triangleIsValid(lAB, lBC, lCA, .1));
    l.o(layoutTriangleVertex(pA, pB, lBC, lCA));
    emit("triangle", l);
  }
  for (int k = 0; k < N; k++) {
    Line l;
    // the lines meet at 15 degrees or more, so the parameters are
    // well-conditioned (nearly parallel lines lose digits to cancellation)
    Vector2 p, a, b, c;
    do {
      p = rand2(), a = rand2(), b = rand2(), c = rand2();
    } while (norm(a) < .25 || std::fabs(cross(unit(a), unit(c - b))) < .25 || std::fabs(cross(unit(a), unit(c))) < .25);
    l.i(p), l.i(a), l.i(b), l.i(c);
    SegmentSegmentIntersectionResult2D s = segmentSegmentIntersection(p, a, b, c);
    RaySegmentIntersectionResult2D r = raySegmentIntersection(p, a, b, c);
    RayRayIntersectionResult2D rr = rayRayIntersection(p, a, b, c);
    l.o(pointLineSegmentDistance(p, a, b)), l.o(s.tA), l.o(s.tB), l.o(s.hit), l.o(r.tRay), l.o(r.tLine);
    l.o(rr.tRay1), l.o(rr.tRay2), l.o(inCircleTest(a, b, c, p));
    emit("segments2", l);
  }
  for (int k = 0; k < N; k++) {
    Line l;
    // no sliver bases or flat tetrahedra: tetFourthPoint's height is the
    // square root of a difference that cancels as d nears the plane of abc
    Vector3 p, a, b, c, d;
    do {
      p = rand3(), a = rand3(), b = rand3(), c = rand3(), d = rand3();
    } while (norm(cross(unit(b - a), unit(c - a))) < .25 || std::fabs(dot(d - a, unit(cross(b - a, c - a)))) < .3);
    l.i(p), l.i(a), l.i(b), l.i(c), l.i(norm(d - a)), l.i(norm(d - b)), l.i(norm(d - c));
    l.o(pointLineSegmentNeaestLocation(p, a, b)), l.o(pointLineSegmentDistance(p, a, b));
    l.o(tetFourthPoint(a, b, c, norm(d - a), norm(d - b), norm(d - c)));
    emit("space3", l);
  }
  return 0;
}
