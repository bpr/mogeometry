# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's include/geometrycentral/utilities/vector2.h
# and vector2.ipp (commit 019669d, MIT License, Copyright (c) 2017-2019
# Nicholas Sharp and the geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Vectors in 2D (`mogeometry.utilities.vector2`)
==============================================
`Vector2`, a vector in the plane that is also a complex number `x + iy`:
`*` and `/` between two vectors are complex multiplication and division,
which is how geometry-central represents tangent vectors and rotations
(multiplying by a unit vector rotates by its angle).

Differences from geometry-central:

- Free functions are methods, as for `Vector3` (see
  `mogeometry.utilities.vector3`): `dot(u, v)` is `u.dot(v)`,
  `cross(u, v)` is `u.cross(v)`, `cross3(u, v)` is `u.cross3(v)`,
  `orientedAngle(u, v)` is `u.oriented_angle(v)`, and so on.
- `rotate(theta)`, deprecated in geometry-central because it rotates
  clockwise, is not ported; use `rotate_cw` or `rotate_ccw`.
- The conversions to and from `std::complex<double>` are
  `Vector2.from_complex` and `to_complex`, with Mojo's `ComplexFloat64`.
- `pow` follows libstdc++'s `std::pow` for a real exponent and glibc's
  `cpow` (`exp(p log z)`) for a complex one, with the C library's `exp`
  and `pow` (`libm_exp`, `libm_pow`), as Mojo's `exp` is inaccurate.
  `log |z|` is computed as `log(hypot(x, y))`, which can differ from
  glibc's `clog` in the last bit when `|z|` is near 1.
- No stream input, `std::hash`, or `to_string`; `Vector2` is `Writable`
  and prints as `<x, y>`.
"""

from std.complex import ComplexFloat64
from std.math import (
    acos,
    atan2,
    cos,
    hypot,
    inf,
    isfinite,
    isnan,
    log,
    nan,
    sin,
    sqrt,
)

from mogeometry.utilities.utilities import (
    clamp,
    fmax,
    fmin,
    libm_exp,
    libm_pow,
)
from mogeometry.utilities.vector3 import Vector3


@fieldwise_init
struct Vector2(Equatable, TrivialRegisterPassable, Writable):
    """A vector in 2D, `(x, y)`, which is also the complex number `x + iy`.

    `+` and `-` are componentwise; `*` and `/` by a scalar scale, and by
    another `Vector2` multiply or divide as complex numbers.
    """

    var x: Float64
    """The first component (the real part)."""
    var y: Float64
    """The second component (the imaginary part)."""

    # ===------------------------------------------------------------------=== #
    # Constructors
    # ===------------------------------------------------------------------=== #

    @staticmethod
    def zero() -> Self:
        """The zero vector.

        Returns:
            `(0, 0)`.
        """
        return Self(0.0, 0.0)

    @staticmethod
    def constant(c: Float64) -> Self:
        """A vector with both components equal to `c`.

        Args:
            c: The value of each component.

        Returns:
            `(c, c)`.
        """
        return Self(c, c)

    @staticmethod
    def from_angle(theta: Float64) -> Self:
        """The unit vector at angle `theta` counterclockwise from the x
        axis.

        Args:
            theta: The angle in radians.

        Returns:
            `(cos theta, sin theta)`.
        """
        return Self(cos(theta), sin(theta))

    @staticmethod
    def from_complex(c: ComplexFloat64) -> Self:
        """The vector of a complex number's real and imaginary parts.

        Args:
            c: The complex number.

        Returns:
            `(c.re, c.im)`.
        """
        return Self(c.re, c.im)

    @staticmethod
    def infinity() -> Self:
        """The vector of positive infinities.

        Returns:
            `(inf, inf)`.
        """
        var i = inf[DType.float64]()
        return Self(i, i)

    @staticmethod
    def undefined() -> Self:
        """The vector of NaNs, marking a value that was not computed.

        Returns:
            `(nan, nan)`.
        """
        var n = nan[DType.float64]()
        return Self(n, n)

    def to_complex(self) -> ComplexFloat64:
        """The vector as a complex number.

        Returns:
            `x + iy`.
        """
        return ComplexFloat64(self.x, self.y)

    # ===------------------------------------------------------------------=== #
    # Indexing and operators
    # ===------------------------------------------------------------------=== #

    def __getitem__(self, index: Int) -> Float64:
        """The component at `index`: 0 is `x`, 1 is `y`.

        Args:
            index: 0 or 1.

        Returns:
            The component.
        """
        debug_assert(0 <= index < 2, "Vector2 index out of range")
        if index == 0:
            return self.x
        return self.y

    def __setitem__(mut self, index: Int, value: Float64):
        """Sets the component at `index`: 0 is `x`, 1 is `y`.

        Args:
            index: 0 or 1.
            value: The new component.
        """
        debug_assert(0 <= index < 2, "Vector2 index out of range")
        if index == 0:
            self.x = value
        else:
            self.y = value

    def __add__(self, v: Self) -> Self:
        """The componentwise sum.

        Args:
            v: The vector to add.

        Returns:
            `self + v`.
        """
        return Self(self.x + v.x, self.y + v.y)

    def __sub__(self, v: Self) -> Self:
        """The componentwise difference.

        Args:
            v: The vector to subtract.

        Returns:
            `self - v`.
        """
        return Self(self.x - v.x, self.y - v.y)

    def __mul__(self, v: Self) -> Self:
        """The complex product.

        Args:
            v: The other factor.

        Returns:
            `self * v` as complex numbers.
        """
        return Self(self.x * v.x - self.y * v.y, self.x * v.y + self.y * v.x)

    def __truediv__(self, v: Self) -> Self:
        """The complex quotient.

        Args:
            v: The divisor; zero gives infinities or NaNs.

        Returns:
            `self / v` as complex numbers.
        """
        var denom = v.x * v.x + v.y * v.y
        return (
            Self(self.x * v.x + self.y * v.y, self.y * v.x - self.x * v.y)
            / denom
        )

    def __mul__(self, s: Float64) -> Self:
        """The vector scaled by `s`.

        Args:
            s: The scale factor.

        Returns:
            `self * s`.
        """
        return Self(self.x * s, self.y * s)

    def __rmul__(self, s: Float64) -> Self:
        """The vector scaled by `s`, written `s * v`.

        Args:
            s: The scale factor.

        Returns:
            `s * self`.
        """
        return Self(s * self.x, s * self.y)

    def __truediv__(self, s: Float64) -> Self:
        """The vector divided by `s`, as a multiplication by `1 / s`.

        Args:
            s: The divisor.

        Returns:
            `self * (1 / s)`.
        """
        var r = 1.0 / s
        return Self(self.x * r, self.y * r)

    def __neg__(self) -> Self:
        """The negated vector.

        Returns:
            `-self`.
        """
        return Self(-self.x, -self.y)

    def __iadd__(mut self, v: Self):
        """Adds `v` in place.

        Args:
            v: The vector to add.
        """
        self.x += v.x
        self.y += v.y

    def __isub__(mut self, v: Self):
        """Subtracts `v` in place.

        Args:
            v: The vector to subtract.
        """
        self.x -= v.x
        self.y -= v.y

    def __imul__(mut self, v: Self):
        """Multiplies by `v` in place, as complex numbers.

        Args:
            v: The other factor.
        """
        self = self * v

    def __itruediv__(mut self, v: Self):
        """Divides by `v` in place, as complex numbers.

        Args:
            v: The divisor.
        """
        self = self / v

    def __imul__(mut self, s: Float64):
        """Scales by `s` in place.

        Args:
            s: The scale factor.
        """
        self.x *= s
        self.y *= s

    def __itruediv__(mut self, s: Float64):
        """Divides each component by `s` in place (a true division, unlike
        `/`).

        Args:
            s: The divisor.
        """
        self.x /= s
        self.y /= s

    def __eq__(self, v: Self) -> Bool:
        """Whether both components are equal.

        Args:
            v: The vector to compare with.

        Returns:
            True if `self` and `v` are equal componentwise.
        """
        return self.x == v.x and self.y == v.y

    def __ne__(self, v: Self) -> Bool:
        """Whether some component differs.

        Args:
            v: The vector to compare with.

        Returns:
            True if `self` and `v` differ in some component.
        """
        return not (self == v)

    def write_to(self, mut writer: Some[Writer]):
        """Writes the vector as `<x, y>`.

        Args:
            writer: The writer to write to.
        """
        writer.write("<", self.x, ", ", self.y, ">")

    # ===------------------------------------------------------------------=== #
    # Norms and normalization
    # ===------------------------------------------------------------------=== #

    def norm(self) -> Float64:
        """The Euclidean length (the complex modulus).

        Returns:
            `|self|`.
        """
        return sqrt(self.x * self.x + self.y * self.y)

    def norm2(self) -> Float64:
        """The squared Euclidean length.

        Returns:
            `|self|^2`.
        """
        return self.x * self.x + self.y * self.y

    def normalize(self) -> Self:
        """The unit vector in the same direction.

        The zero vector gives NaNs.

        Returns:
            `self / |self|`.
        """
        var r = 1.0 / sqrt(self.x * self.x + self.y * self.y)
        return self * r

    def normalize_cutoff(self, mag: Float64 = 0.0) -> Self:
        """The unit vector in the same direction, or `self` unchanged if its
        length is at most `mag`.

        Args:
            mag: Vectors no longer than this are returned unscaled.

        Returns:
            `self / |self|`, or `self`.
        """
        var len = sqrt(self.x * self.x + self.y * self.y)
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
    # Rotations and complex functions
    # ===------------------------------------------------------------------=== #

    def rotate_cw(self, theta: Float64) -> Self:
        """The vector rotated clockwise by `theta`.

        Args:
            theta: The angle in radians.

        Returns:
            The rotated vector.
        """
        var cos_th = cos(theta)
        var sin_th = sin(theta)
        return Self(
            cos_th * self.x + sin_th * self.y,
            -sin_th * self.x + cos_th * self.y,
        )

    def rotate_ccw(self, theta: Float64) -> Self:
        """The vector rotated counterclockwise by `theta`.

        Args:
            theta: The angle in radians.

        Returns:
            The rotated vector.
        """
        return self.rotate_cw(-theta)

    def rotate90(self) -> Self:
        """The vector rotated a quarter turn counterclockwise (times `i`).

        Returns:
            `(-y, x)`.
        """
        return Self(-self.y, self.x)

    def pow(self, p: Float64) -> Self:
        """The complex power `self^p`, on the principal branch.

        Args:
            p: The real exponent.

        Returns:
            `self^p`.
        """
        # libstdc++'s std::pow(complex<double>, double)
        if self.y == 0.0 and self.x > 0.0:
            return Self(libm_pow(self.x, p), 0.0)
        var t_re = log(hypot(self.x, self.y))
        var t_im = atan2(self.y, self.x)
        var rho = libm_exp(p * t_re)
        var theta = p * t_im
        return Self(rho * cos(theta), rho * sin(theta))

    def pow(self, p: Self) -> Self:
        """The complex power `self^p`, on the principal branch.

        Args:
            p: The complex exponent.

        Returns:
            `exp(p log(self))`.
        """
        # glibc's cpow: cexp(p * clog(self))
        var l_re = log(hypot(self.x, self.y))
        var l_im = atan2(self.y, self.x)
        var e_re = p.x * l_re - p.y * l_im
        var e_im = p.x * l_im + p.y * l_re
        var rho = libm_exp(e_re)
        return Self(rho * cos(e_im), rho * sin(e_im))

    def conj(self) -> Self:
        """The complex conjugate.

        Returns:
            `(x, -y)`.
        """
        return Self(self.x, -self.y)

    def inv(self) -> Self:
        """The complex reciprocal.

        Returns:
            `1 / self` as a complex number.
        """
        return Self(1.0, 0.0) / self

    def arg(self) -> Float64:
        """The angle counterclockwise from the x axis (the complex
        argument).

        Returns:
            The angle in `[-pi, pi]`.
        """
        return atan2(self.y, self.x)

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
        return self.x * v.x + self.y * v.y

    def cross(self, v: Self) -> Float64:
        """The 2D cross product, the z component of the 3D one.

        Args:
            v: The other vector.

        Returns:
            `x v.y - y v.x`.
        """
        return self.x * v.y - self.y * v.x

    def cross3(self, v: Self) -> Vector3:
        """The 3D cross product of the two vectors in the xy plane.

        Args:
            v: The other vector.

        Returns:
            `(0, 0, x v.y - y v.x)`.
        """
        return Vector3(0.0, 0.0, self.x * v.y - self.y * v.x)

    def sum(self) -> Float64:
        """The sum of the components.

        Returns:
            `x + y`.
        """
        return self.x + self.y

    def angle(self, v: Self) -> Float64:
        """The unsigned angle between `self` and `v`.

        Args:
            v: The other vector.

        Returns:
            The angle in `[0, pi]`.
        """
        return acos(fmax(-1.0, fmin(1.0, self.unit().dot(v.unit()))))

    def oriented_angle(self, v: Self) -> Float64:
        """The signed angle from `self` to `v`, counterclockwise positive.

        Args:
            v: The vector the angle is measured to.

        Returns:
            The angle in `[-pi, pi]`.
        """
        var u_hat = self.unit()
        var v_hat = v.unit()
        return atan2(u_hat.cross(v_hat), u_hat.dot(v_hat))

    # ===------------------------------------------------------------------=== #
    # Other operations
    # ===------------------------------------------------------------------=== #

    def is_finite(self) -> Bool:
        """Whether both components are finite.

        Returns:
            True if no component is infinite or NaN.
        """
        return isfinite(self.x) and isfinite(self.y)

    def is_defined(self) -> Bool:
        """Whether no component is NaN.

        Returns:
            True if no component is NaN.
        """
        return not isnan(self.x) and not isnan(self.y)

    def clamp(self, low: Self, high: Self) -> Self:
        """The vector clamped componentwise to `[low, high]`.

        Args:
            low: The lower bounds.
            high: The upper bounds.

        Returns:
            Each component limited to its interval.
        """
        return Self(clamp(self.x, low.x, high.x), clamp(self.y, low.y, high.y))

    def componentwise_min(self, v: Self) -> Self:
        """The componentwise minimum, ignoring NaNs (C's `fmin`).

        Args:
            v: The other vector.

        Returns:
            The smaller of each pair of components.
        """
        return Self(fmin(self.x, v.x), fmin(self.y, v.y))

    def componentwise_max(self, v: Self) -> Self:
        """The componentwise maximum, ignoring NaNs (C's `fmax`).

        Args:
            v: The other vector.

        Returns:
            The larger of each pair of components.
        """
        return Self(fmax(self.x, v.x), fmax(self.y, v.y))
