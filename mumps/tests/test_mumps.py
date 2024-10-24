# Copyright 2011-2016 Kwant authors.
# Copyright 2018 Python-MUMPS Authors.
#
# This file is part of Python-MUMPS. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# Python-MUMPS authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at
# https://gitlab.kwant-project.org/kwant/python-mumps.

import pytest
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

from mumps import Context, MUMPSError, schur_complement, nullspace
from ._test_utils import _Random

# Decimal places of precision per datatype. These limits have been determined
# heuristically by inspecting the upper error bound reported by MUMPS for
# the random matrices used here.
precisions = {
    np.float32: 1,  # yes, really only 1 decimal place of precision!
    np.float64: 10,
    np.complex64: 1,
    np.complex128: 10,
}

dtypes = list(precisions.keys())


def assert_array_almost_equal(dtype, a, b):
    np.testing.assert_almost_equal(a, b, decimal=precisions[dtype])


def assert_almost_equal(dtype, a, b):
    np.testing.assert_almost_equal(a, b, decimal=precisions[dtype])


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [2, 10, 100], ids=str)
def test_lu_with_dense(dtype, mat_size):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    bmat = rand.randmat(mat_size, mat_size, dtype)
    bvec = rand.randvec(mat_size, dtype)

    ctx = Context()
    ctx.factor(sp.coo_matrix(a))

    xvec = ctx.solve(bvec)
    xmat = ctx.solve(bmat)

    assert_array_almost_equal(dtype, np.dot(a, xvec), bvec)
    assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)

    # now "sparse" right hand side
    xvec = ctx.solve(sp.csc_matrix(bvec.reshape(mat_size, 1)))
    xmat = ctx.solve(sp.csc_matrix(bmat))

    assert_array_almost_equal(dtype, np.dot(a, xvec), bvec.reshape(mat_size, 1))
    assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [5, 10], ids=str)
def test_schur_complement_with_dense(dtype, mat_size):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    s = schur_complement(sp.coo_matrix(a), list(range(3)))
    assert_array_almost_equal(dtype, np.linalg.inv(s), np.linalg.inv(a)[:3, :3])


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [5, 10], ids=str)
@pytest.mark.parametrize("symmetric_matrix", [True, False], ids=str)
def test_schur_complement_solution(dtype, mat_size, symmetric_matrix):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    if symmetric_matrix:
        a = a + a.T

    bvec = rand.randvec(mat_size, dtype)

    ctx = Context()
    ctx.set_matrix(a, symmetric=symmetric_matrix)
    ctx.schur(range(3))

    xvec = ctx.solve_schur(bvec)

    assert_array_almost_equal(dtype, a @ xvec, bvec)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_factor_error(dtype):
    """Test that an error is raised if factor is asked to reuse missing analysis."""
    a = sp.identity(10, dtype=dtype)
    with pytest.raises(ValueError):
        Context().factor(a, reuse_analysis=True)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_error_minus_19(dtype):
    """Test if MUMPSError -19 is properly caught by increasing memory"""
    ctx = Context()

    a = sp.eye(5000, dtype=dtype)

    # Create the context so we can modify it for factorization
    ctx.analyze(a)

    # We ensure that this first call creates a -19 error by
    # allocating only 1 MB for factorization
    ctx.mumps_instance.icntl[23] = 1  # Memory upper bound set to 1MB
    ctx.mumps_instance.icntl[14] = 1  # Initial memory relaxation to 1%
    ctx.mumps_instance.job = 2
    ctx.mumps_instance.call()
    # ensure that we really don't have enough memory
    assert ctx.mumps_instance.infog[1] == -19

    # This call should not raise any errors as
    # it would successfully allocate memory
    Context().factor(a)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [20, 50], ids=str)
@pytest.mark.parametrize("symmetric_matrix", [True, False], ids=str)
def test_nullspace(dtype, mat_size, symmetric_matrix):
    """Tests the nullspace wrapper by creating a rank deficient matrix
    with known deficiency

    This test is to check whether the wrapper works. The algorithm
    to obtain the nullspace basis within MUMPS is still rather
    unclear with the use of these null pivot threshold. Suppose
    that we have to live with this uncertainty, we at least create
    a test case where we know for sure the right nullspace
    dimensions and check if the subroutine has successfully
    calculated the nullspace basis.
    """

    if not symmetric_matrix:
        test = sp.coo_matrix([[1, 2], [2, 4]])
        try:
            nullspace(test, False)
        except MUMPSError as error:
            if error.error == -37:
                pytest.skip("Installed MUMPS does not support unsymmetric nullspace")
            else:
                raise

    # We ensure the dimension of the nullspace row vectors by creating a random
    # column vector, a of size (mat_size x 1)
    rand = _Random()
    a = sp.coo_matrix(rand.randmat(mat_size, 1, dtype))

    # We then create the symmetric matrix by making this operation:
    # A = a * a^T, where A has size (mat_size x mat_size). Also because the
    # matrix is created by an outer product of two column vectors, the rank of
    # A is always 1. Hence, we know the dimension of the nullspace row vectors
    # is always mat_size - 1. Also, since the two column vectors are the same,
    # A is then a symmetric matrix.
    if symmetric_matrix:
        b = a.copy()
    else:
        b = sp.coo_matrix(rand.randmat(mat_size, 1, dtype))

    A = a * b.transpose()
    A_original = a * b.transpose()

    # Remove the lower triangular elements in order to comply with the
    # symmetric solver of MUMPS
    if symmetric_matrix:
        for i in range(a.shape[0]):
            for j in range(i):
                A[i, j] = 0
        A.eliminate_zeros()

    # Start with a small pivot_threshold value
    pivot_threshold = 1e-5
    maximum_threshold = 1.0

    # We repeatedly add the threshold until we get the full nullspace basis
    while True:
        X = nullspace(A, symmetric_matrix, pivot_threshold)
        pivot_threshold += 2 * 1e-5
        if (X.shape[1] == A.shape[0] - 1) or (pivot_threshold > maximum_threshold):
            break

    # We know that the row vector dimension of the returned nullspace
    # has to be mat_size-1
    assert X.shape[1] == A_original.shape[0] - 1

    # Check that indeed the nullspace does indeed give 0 when multiplied
    # with the input matrix
    assert_almost_equal(dtype, la.norm(A_original * X), 0.0)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_one_by_one(dtype):
    """Test a 1x1 matrix.

    This is a regression test for
    https://gitlab.kwant-project.org/kwant/python-mumps/-/issues/18
    """
    ctx = Context()
    ctx.factor(sp.eye(1, dtype=dtype))
    assert_almost_equal(dtype, ctx.solve(np.array([1], dtype=dtype)), 1)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_zero_size_rhs(dtype):
    """Test that a 0xn rhs can be solved."""
    a = np.random.randn(10, 10).astype(dtype)
    ctx = Context()
    ctx.factor(a)
    rhs = np.zeros((10, 0), dtype=dtype)
    assert_almost_equal(dtype, ctx.solve(rhs), rhs)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_symmetric_matrix(dtype):
    """Test that a symmetric matrix can be solved."""
    n = 10
    a = np.random.randn(n, n).astype(dtype)
    if np.iscomplexobj(a):
        a += (1j * np.random.randn(n, n)).astype(dtype)
    a += a.T
    ctx = Context()
    ctx.set_matrix(a, symmetric=True)
    ctx.factor()
    rhs = np.random.randn(n, 1).astype(dtype)
    assert_almost_equal(dtype, ctx.solve(rhs), la.solve(a, rhs))
