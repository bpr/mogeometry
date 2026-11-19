# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's include/geometrycentral/utilities/vector3.h
# and vector3.ipp (commit 019669d, MIT License, Copyright (c) 2017-2019
# Nicholas Sharp and the geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Vectors in 3D (`mogeometry.utilities.vector3`)
==============================================
`Vector3`, a point or direction in space: three `Float64` fields, with the
arithmetic, norms, products, and angles that geometry code needs.

Differences from geometry-central:

- Free functions are methods: `dot(u, v)` is `u.dot(v)`, `norm(v)` is
  `v.norm()`, `cross(u, v)` is `u.cross(v)`, `isfinite(v)` is
  `v.is_finite()`, `angleInPlane(u, v, n)` is `u.angle_in_plane(v, n)`,
  `componentwiseMin(u, v)` is `u.componentwise_min(v)`, and
  `clamp(v, lo, hi)` is `v.clamp(lo, hi)`. Mojo does not merge
  overloads imported from different modules, so a free `dot` for
  `Vector2` and one for `Vector3` could not both be imported.
- `buildTangentBasis()` returns a `TangentBasis` with named fields
  instead of a `std::array`.
- No stream input, `std::hash`, or `to_string`; `Vector3` is `Writable`
  and prints as geometry-central's `operator<<` does, `<x, y, z>`.

The layout (three named fields rather than a `SIMD[DType.float64, 4]`)
is provisional: `docs/plan.md` (Performance experiments) records the
SIMD comparison to run.
"""

from std.math import acos, atan2, cos, inf, isfinite, isnan, nan, sin, sqrt

from mogeometry.utilities.utilities import clamp, fmax, fmin


@fieldwise_init
struct TangentBasis(TrivialRegisterPassable, Writable):
    """Two unit vectors that, with a given normal, form an orthonormal
    frame, as `Vector3.build_tangent_basis()` returns them."""

    var basis_x: Vector3
    """The first tangent direction."""
    var basis_y: Vector3
    """The second tangent direction: the normal crossed with `basis_x`."""

    def write_to(self, mut writer: Some[Writer]):
        """Writes the two vectors as `(<x>, <y>)`.

        Args:
            writer: The writer to write to.
        """
        writer.write("(", self.basis_x, ", ", self.basis_y, ")")


@fieldwise_init
struct Vector3(Equatable, TrivialRegisterPassable, Writable):
    """A vector in 3D, `(x, y, z)`.

    Arithmetic is componentwise; `*` and `/` take a scalar. Dividing by a
    scalar multiplies by its reciprocal, as geometry-central does, except
    for `/=`, which divides each component.
    """

    var x: Float64
    """The first component."""
    var y: Float64
    """The second component."""
    var z: Float64
    """The third component."""

    # ===------------------------------------------------------------------=== #
    # Constructors
    # ===------------------------------------------------------------------=== #

    @staticmethod
    def zero() -> Self:
        """The zero vector.

        Returns:
            `(0, 0, 0)`.
        """
        return Self(0.0, 0.0, 0.0)

    @staticmethod
    def constant(c: Float64) -> Self:
        """A vector with every component equal to `c`.

        Args:
            c: The value of each component.

        Returns:
            `(c, c, c)`.
        """
        return Self(c, c, c)

    @staticmethod
    def infinity() -> Self:
        """The vector of positive infinities.

        Returns:
            `(inf, inf, inf)`.
        """
        var i = inf[DType.float64]()
        return Self(i, i, i)

    @staticmethod
    def undefined() -> Self:
        """The vector of NaNs, marking a value that was not computed.

        Returns:
            `(nan, nan, nan)`.
        """
        var n = nan[DType.float64]()
        return Self(n, n, n)

    # ===------------------------------------------------------------------=== #
    # Indexing and operators
    # ===------------------------------------------------------------------=== #

    def __getitem__(self, index: Int) -> Float64:
        """The component at `index`: 0 is `x`, 1 is `y`, 2 is `z`.

        Args:
            index: 0, 1, or 2.

        Returns:
            The component.
        """
        debug_assert(0 <= index < 3, "Vector3 index out of range")
        if index == 0:
            return self.x
        if index == 1:
            return self.y
        return self.z

    def __setitem__(mut self, index: Int, value: Float64):
        """Sets the component at `index`: 0 is `x`, 1 is `y`, 2 is `z`.

        Args:
            index: 0, 1, or 2.
            value: The new component.
        """
        debug_assert(0 <= index < 3, "Vector3 index out of range")
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        else:
            self.z = value

    def __add__(self, v: Self) -> Self:
        """The componentwise sum.

        Args:
            v: The vector to add.

        Returns:
            `self + v`.
        """
        return Self(self.x + v.x, self.y + v.y, self.z + v.z)

    def __sub__(self, v: Self) -> Self:
        """The componentwise difference.

        Args:
            v: The vector to subtract.

        Returns:
            `self - v`.
        """
        return Self(self.x - v.x, self.y - v.y, self.z - v.z)

    def __mul__(self, s: Float64) -> Self:
        """The vector scaled by `s`.

        Args:
            s: The scale factor.

        Returns:
            `self * s`.
        """
        return Self(self.x * s, self.y * s, self.z * s)

    def __rmul__(self, s: Float64) -> Self:
        """The vector scaled by `s`, written `s * v`.

        Args:
            s: The scale factor.

        Returns:
            `s * self`.
        """
        return Self(s * self.x, s * self.y, s * self.z)

    def __truediv__(self, s: Float64) -> Self:
        """The vector divided by `s`, as a multiplication by `1 / s`.

        Args:
            s: The divisor.

        Returns:
            `self * (1 / s)`.
        """
        var r = 1.0 / s
        return Self(self.x * r, self.y * r, self.z * r)

    def __neg__(self) -> Self:
        """The negated vector.

        Returns:
            `-self`.
        """
        return Self(-self.x, -self.y, -self.z)

    def __iadd__(mut self, v: Self):
        """Adds `v` in place.

        Args:
            v: The vector to add.
        """
        self.x += v.x
        self.y += v.y
        self.z += v.z

    def __isub__(mut self, v: Self):
        """Subtracts `v` in place.

        Args:
            v: The vector to subtract.
        """
        self.x -= v.x
        self.y -= v.y
        self.z -= v.z

    def __imul__(mut self, s: Float64):
        """Scales by `s` in place.

        Args:
            s: The scale factor.
        """
        self.x *= s
        self.y *= s
        self.z *= s

    def __itruediv__(mut self, s: Float64):
        """Divides each component by `s` in place (a true division, unlike
        `/`).

        Args:
            s: The divisor.
        """
        self.x /= s
        self.y /= s
        self.z /= s

    def __eq__(self, v: Self) -> Bool:
        """Whether every component is equal.

        Args:
            v: The vector to compare with.

        Returns:
            True if `self` and `v` are equal componentwise.
        """
        return self.x == v.x and self.y == v.y and self.z == v.z

    def __ne__(self, v: Self) -> Bool:
        """Whether some component differs.

        Args:
            v: The vector to compare with.

        Returns:
            True if `self` and `v` differ in some component.
        """
        return not (self == v)

    def write_to(self, mut writer: Some[Writer]):
        """Writes the vector as `<x, y, z>`.

        Args:
            writer: The writer to write to.
        """
        writer.write("<", self.x, ", ", self.y, ", ", self.z, ">")

    # ===------------------------------------------------------------------=== #
    # Norms and normalization
    # ===------------------------------------------------------------------=== #

    def norm(self) -> Float64:
        """The Euclidean length.

        Returns:
            `|self|`.
        """
        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def norm2(self) -> Float64:
        """The squared Euclidean length.

        Returns:
            `|self|^2`.
        """
        return self.x * self.x + self.y * self.y + self.z * self.z

    def normalize(self) -> Self:
        """The unit vector in the same direction.

        The zero vector gives NaNs.

        Returns:
            `self / |self|`.
        """
        var r = 1.0 / sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        return self * r

    def normalize_cutoff(self, mag: Float64 = 0.0) -> Self:
        """The unit vector in the same direction, or `self` unchanged if its
        length is at most `mag`.

        Args:
            mag: Vectors no longer than this are returned unscaled.

        Returns:
            `self / |self|`, or `self`.
        """
        var len = sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        if len <= mag:
            len = 1.0
        var r = 1.0 / len
        return self * r

    def unit(self) -> Self:
        """The unit vector in the same direction (`normalize`).

        Returns:
            `self / |self|`.
        """
        return self.normalize()

    # ===------------------------------------------------------------------=== #
    # Products and angles
    # ===------------------------------------------------------------------=== #

    def dot(self, v: Self) -> Float64:
        """The dot product.

        Args:
            v: The other vector.

        Returns:
            `self . v`.
        """
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cross(self, v: Self) -> Self:
        """The cross product.

        Args:
            v: The other vector.

        Returns:
            `self x v`.
        """
        var x = self.y * v.z - self.z * v.y
        var y = self.z * v.x - self.x * v.z
        var z = self.x * v.y - self.y * v.x
        return Self(x, y, z)

    def sum(self) -> Float64:
        """The sum of the components.

        Returns:
            `x + y + z`.
        """
        return self.x + self.y + self.z

    def angle(self, v: Self) -> Float64:
        """The unsigned angle between `self` and `v`.

        Args:
            v: The other vector.

        Returns:
            The angle in `[0, pi]`.
        """
        return acos(fmax(-1.0, fmin(1.0, self.unit().dot(v.unit()))))

    def angle_in_plane(self, v: Self, normal: Self) -> Float64:
        """The signed angle from `self` to `v`, measured counterclockwise
        about `normal` in the plane orthogonal to it.

        Args:
            v: The vector the angle is measured to.
            normal: The axis of the plane; need not be unit.

        Returns:
            The angle in `[-pi, pi]`.
        """
        # put self in the plane with the normal
        var n = normal.unit()
        var u_plane = (self - self.dot(n) * n).unit()
        var basis_y = normal.cross(u_plane).unit()
        var x_comp = v.dot(u_plane)
        var y_comp = v.dot(basis_y)
        return atan2(y_comp, x_comp)

    # ===------------------------------------------------------------------=== #
    # Other operations
    # ===------------------------------------------------------------------=== #

    def rotate_around(self, axis: Self, theta: Float64) -> Self:
        """The vector rotated counterclockwise by `theta` about `axis`.

        Args:
            axis: The axis of rotation; need not be unit.
            theta: The angle in radians.

        Returns:
            The rotated vector.
        """
        var axis_n = axis.normalize()
        var parallel_comp = axis_n * self.dot(axis_n)
        var tangent_comp = self - parallel_comp
        if tangent_comp.norm2() > 0.0:
            var basis_x = tangent_comp.normalize()
            var basis_y = axis_n.cross(basis_x)
            var tangent_mag = tangent_comp.norm()
            var rotated = tangent_mag * (
                cos(theta) * basis_x + sin(theta) * basis_y
            )
            return rotated + parallel_comp
        return parallel_comp

    def remove_component(self, unit_dir: Self) -> Self:
        """The vector with its component along `unit_dir` removed.

        Args:
            unit_dir: A unit vector.

        Returns:
            `self - (unit_dir . self) unit_dir`.
        """
        return self - unit_dir * unit_dir.dot(self)

    def build_tangent_basis(self) -> TangentBasis:
        """Two unit vectors orthogonal to `self` and to each other.

        `self` need not be unit. With `n = self.unit()`, the result
        `(basis_x, basis_y)` satisfies `basis_x x basis_y = n`.

        Returns:
            The tangent basis.
        """
        var unit_dir = self.normalize()
        var test_vec = Self(1.0, 0.0, 0.0)
        if abs(test_vec.dot(unit_dir)) > 0.9:
            test_vec = Self(0.0, 1.0, 0.0)
        var basis_x = test_vec.cross(unit_dir).normalize()
        var basis_y = unit_dir.cross(basis_x).normalize()
        return TangentBasis(basis_x, basis_y)

    def is_finite(self) -> Bool:
        """Whether every component is finite.

        Returns:
            True if no component is infinite or NaN.
        """
        return isfinite(self.x) and isfinite(self.y) and isfinite(self.z)

    def is_defined(self) -> Bool:
        """Whether no component is NaN.

        Returns:
            True if no component is NaN.
        """
        return not isnan(self.x) and not isnan(self.y) and not isnan(self.z)

    def clamp(self, low: Self, high: Self) -> Self:
        """The vector clamped componentwise to `[low, high]`.

        Args:
            low: The lower bounds.
            high: The upper bounds.

        Returns:
            Each component limited to its interval.
        """
        return Self(
            clamp(self.x, low.x, high.x),
            clamp(self.y, low.y, high.y),
            clamp(self.z, low.z, high.z),
        )

    def componentwise_min(self, v: Self) -> Self:
        """The componentwise minimum, ignoring NaNs (C's `fmin`).

        Args:
            v: The other vector.

        Returns:
            The smaller of each pair of components.
        """
        return Self(fmin(self.x, v.x), fmin(self.y, v.y), fmin(self.z, v.z))

    def componentwise_max(self, v: Self) -> Self:
        """The componentwise maximum, ignoring NaNs (C's `fmax`).

        Args:
            v: The other vector.

        Returns:
            The larger of each pair of components.
        """
        return Self(fmax(self.x, v.x), fmax(self.y, v.y), fmax(self.z, v.z))
