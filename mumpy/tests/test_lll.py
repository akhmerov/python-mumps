import numpy as np
from mumpy import lll
from kwant._common import ensure_rng

def test_lll():
    rng = ensure_rng(1)
    for i in range(50):
        x = rng.randint(4) + 1
        mat = rng.randn(x, x + rng.randint(2))
        c = 1.34 + .5 * rng.random_sample()
        reduced_mat, coefs = lll.lll(mat)
        assert lll.is_c_reduced(reduced_mat, c)
        assert np.allclose(np.dot(mat.T, coefs), reduced_mat.T)


def test_cvp():
    rng = ensure_rng(0)
    for i in range(10):
        mat = rng.randn(4, 4)
        mat = lll.lll(mat)[0]
        for j in range(4):
            point = 50 * rng.randn(4)
            assert np.array_equal(lll.cvp(point, mat, 10)[:3],
                                  lll.cvp(point, mat, 3))
