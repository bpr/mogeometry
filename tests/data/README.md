# Reference data

Expected results from geometry-central itself, one directory per area,
each with the C++ driver that wrote it:

| Directory    | Driver                     | Tested by                    |
| ------------ | -------------------------- | ---------------------------- |
| `utilities/` | `utilities/reference.cpp`  | `tests/test_utilities.mojo`  |

Rebuild with `pixi run -e ref ref-data` (`build_reference.sh`). It
compiles each driver, with any geometry-central sources it needs,
against the checkout at `$GEOMETRY_CENTRAL_PATH` (default
`../../c++/geometry-central`). That checkout should be at the commit in
`NOTICE`; the script prints the commit it used. The `ref` pixi
environment provides the compiler and Eigen 3.3.9, the version
geometry-central's CMake downloads. Drivers are built with
`-ffp-contract=off`, so every operation rounds as the C++ is written.

## Format

`reference.txt` has one case per line:

```
name in1 in2 ... -> out1 out2 ...
```

Every number is a double written as its IEEE 754 bit pattern, as a
signed 64-bit decimal integer, so it reads back exactly (Mojo: `atol`,
then `bitcast[DType.float64]`). Booleans are the doubles 0 and 1. The
driver's `l.i(...)` and `l.o(...)` calls list each name's inputs and
outputs in order.
Inputs are random, from a fixed seed, and filtered to be well-conditioned
(see the driver), so the tests can use tight tolerances.
