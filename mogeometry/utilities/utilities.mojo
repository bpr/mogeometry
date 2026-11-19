# ===----------------------------------------------------------------------=== #
# mogeometry: geometry processing on surface meshes and point clouds in Mojo.
# Ported from geometry-central's include/geometrycentral/utilities/utilities.h
# (commit 019669d, MIT License, Copyright (c) 2017-2019 Nicholas Sharp and the
# geometry-central contributors); see NOTICE.
# ===----------------------------------------------------------------------=== #
"""
Scalar Helpers (`mogeometry.utilities.utilities`)
=================================================
The constants and scalar arithmetic of geometry-central's `utilities.h`
that the rest of the port uses: `PI`, `clamp`, `regularize_angle`, C's
`fmin` and `fmax` (which geometry-central calls through `std::`), and
`libm_pow` and `libm_exp`, the C library's `pow` and `exp`, which are
more accurate than Mojo's.

Not ported: the memory, type-name, and permutation helpers (C++
plumbing), `INVALID_IND` (handles will say which index is invalid), and
the random-number helpers, which arrive with their first user.
"""

from std.ffi import external_call
from std.math import floor, isnan

comptime PI: Float64 = 3.1415926535897932384
"""Pi, as geometry-central writes it (the nearest double)."""


def clamp(val: Float64, low: Float64, high: Float64) -> Float64:
    """`val` limited to the interval `[low, high]`.

    A NaN `val` passes through unchanged.

    Args:
        val: The value to clamp.
        low: The lower bound.
        high: The upper bound.

    Returns:
        `high` if `val > high`, `low` if `val < low`, else `val`.
    """
    if val > high:
        return high
    if val < low:
        return low
    return val


def regularize_angle(theta: Float64) -> Float64:
    """`theta` mapped into `[0, 2 pi)` by adding a multiple of `2 pi`.

    Args:
        theta: An angle in radians.

    Returns:
        The equivalent angle in `[0, 2 pi)`.
    """
    return theta - 2 * PI * floor(theta / (2 * PI))


def libm_pow(x: Float64, y: Float64) -> Float64:
    """`x` to the power `y`, by the C library's `pow`.

    Mojo's `**` with a fractional exponent is off by up to ~5e-12
    relative (see AGENTS.md, Known toolchain bugs); this is correctly
    rounded to within an ulp, as geometry-central's `std::pow` is.

    Args:
        x: The base.
        y: The exponent.

    Returns:
        `x ** y`.
    """
    return external_call["pow", Float64](x, y)


def libm_exp(x: Float64) -> Float64:
    """`e` to the power `x`, by the C library's `exp`.

    Mojo's `std.math.exp` on `Float64` is off by up to ~35,000 ulps
    (~8e-12 relative; see AGENTS.md, Known toolchain bugs); this is
    within an ulp, as geometry-central's `std::exp` is.

    Args:
        x: The exponent.

    Returns:
        `e ** x`.
    """
    return external_call["exp", Float64](x)


def fmin(a: Float64, b: Float64) -> Float64:
    """The smaller of `a` and `b`, ignoring a NaN, as C's `fmin`.

    Args:
        a: A number.
        b: A number.

    Returns:
        The smaller argument; the other one if either is NaN.
    """
    if isnan(a):
        return b
    if isnan(b):
        return a
    return a if a < b else b


def fmax(a: Float64, b: Float64) -> Float64:
    """The larger of `a` and `b`, ignoring a NaN, as C's `fmax`.

    Args:
        a: A number.
        b: A number.

    Returns:
        The larger argument; the other one if either is NaN.
    """
    if isnan(a):
        return b
    if isnan(b):
        return a
    return a if a > b else b
