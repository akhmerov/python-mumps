# Copyright 2011-2016 Anton Akhmerov, Christoph Groth, and Michael Wimmer and
# Copyright 2017 Bas Nijholt.
#
# This file is part of mumpy. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# mumpy authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at https://github.com/basnijholt/mumpy

import numpy as np
from mumpy import lll


def ensure_rng(rng=None):
    """Turn rng into a random number generator instance.

    If rng is None, return the RandomState instance used by np.random.
    If rng is an integer, return a new RandomState instance seeded with rng.
    If rng is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if rng is None:
        return np.random.mtrand._rand
    if isinstance(rng, numbers.Integral):
        return np.random.RandomState(rng)
    if all(hasattr(rng, attr) for attr in ('random_sample', 'randn',
                                           'randint', 'choice')):
        return rng
    raise ValueError("Expecting a seed or an object that offers the "
                     "numpy.random.RandomState interface.")


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
