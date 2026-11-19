# The build works: mogeometry imports, and MojoDEC's sparse solvers link.

from std.testing import TestSuite, assert_almost_equal, assert_equal

import mogeometry
from mojodec import SparseCholesky
from mojodec.core.sparse import CSRMatrix
from mojodec.core.utils import real_vector


def test_mojodec_solves() raises:
    # the 1D Laplacian [2 -1 0; -1 2 -1; 0 -1 2] x = (1, 0, 1): x = (1, 1, 1)
    var rows: List[Int] = [0, 0, 1, 1, 1, 2, 2]
    var cols: List[Int] = [0, 1, 0, 1, 2, 1, 2]
    var vals: List[Float64] = [2, -1, -1, 2, -1, -1, 2]
    var A = CSRMatrix.from_coo(3, 3, rows, cols, vals)
    var b = real_vector(3)
    b.unsafe_set(0, 1)
    b.unsafe_set(2, 1)
    var x = SparseCholesky(A).solve(b)
    assert_equal(len(x), 3)
    for i in range(3):
        assert_almost_equal(x.unsafe_get(i), 1, atol=1e-15, rtol=0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
