# Copyright 2011-2016 Anton Akhmerov, Christoph Groth, and Michael Wimmer and
# Copyright 2017 Bas Nijholt.
#
# This file is part of mumpy. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# mumpy authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at https://github.com/basnijholt/mumpy

from mumpy import MUMPSContext, schur_complement
from kwant.lattice import honeycomb
from kwant.builder import Builder, HoppingKind
import pytest
import numpy as np
import scipy.sparse as sp
from ._test_utils import _Random, assert_array_almost_equal


def test_lu_with_dense():
    def _test_lu_with_dense(dtype):
        rand = _Random()
        a = rand.randmat(5, 5, dtype)
        bmat = rand.randmat(5, 5, dtype)
        bvec = rand.randvec(5, dtype)

        ctx = MUMPSContext()
        ctx.factor(sp.coo_matrix(a))

        xvec = ctx.solve(bvec)
        xmat = ctx.solve(bmat)

        assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)
        assert_array_almost_equal(dtype, np.dot(a, xvec), bvec)

        # now "sparse" right hand side

        xvec = ctx.solve(sp.csc_matrix(bvec.reshape(5,1)))
        xmat = ctx.solve(sp.csc_matrix(bmat))

        assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)
        assert_array_almost_equal(dtype, np.dot(a, xvec),
                                  bvec.reshape(5,1))

    _test_lu_with_dense(np.complex128)


def test_schur_complement_with_dense():
    def _test_schur_complement_with_dense(dtype):
        rand = _Random()
        a = rand.randmat(10, 10, dtype)
        s = schur_complement(sp.coo_matrix(a), list(range(3)))
        assert_array_almost_equal(dtype, np.linalg.inv(s),
                                  np.linalg.inv(a)[:3, :3])

    _test_schur_complement_with_dense(np.complex128)


def test_error_minus_9(r=10):
    """Test if MUMPSError -9 is properly caught by increasing memory"""

    graphene = honeycomb()
    a, b = graphene.sublattices

    def circle(pos):
        x, y = pos
        return x**2 + y**2 < r**2

    syst = Builder()
    syst[graphene.shape(circle, (0,0))] = -0.0001
    for kind in [((0, 0), b, a), ((0, 1), b, a), ((-1, 1), b, a)]:
        syst[HoppingKind(*kind)] = - 1

    ham = syst.finalized().hamiltonian_submatrix(sparse=True)

    # No need to check result, it's enough if no exception is raised
    MUMPSContext().factor(ham)


def test_factor_warning():
    """Test that a warning is raised if factor is asked without analysis."""
    a = sp.identity(10, dtype=complex)
    with pytest.warns(RuntimeWarning):
        MUMPSContext().factor(a, reuse_analysis=True)
