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

from mumps import Context, MUMPSError, schur_complement, nullspace, complex_to_real
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
@pytest.mark.parametrize("blr", [False, True], ids=str)
def test_lu_with_dense(dtype, mat_size, blr):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    bmat = rand.randmat(mat_size, mat_size, dtype)
    bvec = rand.randvec(mat_size, dtype)

    ctx = Context()
    if blr:
        ctx.activate_blr(option=1)
    ctx.factor(sp.coo_matrix(a))

    xvec = ctx.solve(bvec)
    xmat = ctx.solve(bmat)
    if ctx.myid == 0:
        assert_array_almost_equal(dtype, np.dot(a, xvec), bvec)
        assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)

    # now "sparse" right hand side
    xvec = ctx.solve(sp.csc_matrix(bvec.reshape(mat_size, 1)))
    xmat = ctx.solve(sp.csc_matrix(bmat))
    if ctx.myid == 0:
        assert_array_almost_equal(dtype, np.dot(a, xvec), bvec.reshape(mat_size, 1))
        assert_array_almost_equal(dtype, np.dot(a, xmat), bmat)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_sparse_array_rhs_vector(dtype):
    rand = _Random()
    mat_size = 8
    a = rand.randmat(mat_size, mat_size, dtype)
    a += (mat_size + 1) * np.eye(mat_size, dtype=dtype)
    bvec = rand.randvec(mat_size, dtype)

    coords = np.arange(mat_size, dtype=int)
    b_sparse = sp.coo_array((bvec, (coords,)), shape=(mat_size,))
    assert b_sparse.ndim == 1

    ctx = Context()
    ctx.factor(sp.coo_array(a))
    x = ctx.solve(b_sparse)

    assert x.ndim == 1
    assert_array_almost_equal(dtype, a @ x, bvec)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_sparse_rhs_multiple_columns(dtype):
    rand = _Random()
    mat_size = 6
    a = rand.randmat(mat_size, mat_size, dtype)
    a += (mat_size + 1) * np.eye(mat_size, dtype=dtype)

    rows = np.array([0, 2, 5, 1, 4, 3], dtype=int)
    cols = np.array([0, 0, 0, 1, 1, 2], dtype=int)
    vals = rand.randvec(rows.size, dtype)
    b_sparse = sp.coo_matrix((vals, (rows, cols)), shape=(mat_size, 3)).tocsc()

    ctx = Context()
    ctx.factor(sp.coo_matrix(a))
    x = ctx.solve(b_sparse)

    assert_array_almost_equal(dtype, a @ x, b_sparse.toarray())


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [5, 10], ids=str)
@pytest.mark.mpi_skip
def test_schur_complement_with_dense(dtype, mat_size):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    s = schur_complement(sp.coo_matrix(a), list(range(3)))
    assert_array_almost_equal(dtype, np.linalg.inv(s), np.linalg.inv(a)[:3, :3])


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [5, 10], ids=str)
@pytest.mark.parametrize("symmetric_matrix", [True, False], ids=str)
@pytest.mark.mpi_skip
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
def test_schur_sparse_rhs_multiple_columns(dtype):
    rand = _Random()
    mat_size = 6
    a = rand.randmat(mat_size, mat_size, dtype)
    a += (mat_size + 1) * np.eye(mat_size, dtype=dtype)

    rows = np.array([0, 2, 5, 1, 4, 3], dtype=int)
    cols = np.array([0, 0, 0, 1, 1, 2], dtype=int)
    vals = rand.randvec(rows.size, dtype)
    b_sparse = sp.coo_matrix((vals, (rows, cols)), shape=(mat_size, 3)).tocsc()

    ctx = Context()
    ctx.set_matrix(a)
    ctx.schur(range(2))
    x = ctx.solve_schur(b_sparse)

    assert_array_almost_equal(dtype, a @ x, b_sparse.toarray())


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_schur_sparse_array_rhs_vector(dtype):
    rand = _Random()
    mat_size = 6
    a = rand.randmat(mat_size, mat_size, dtype)
    a += (mat_size + 1) * np.eye(mat_size, dtype=dtype)
    bvec = rand.randvec(mat_size, dtype)

    coords = np.arange(mat_size, dtype=int)
    b_sparse = sp.coo_array((bvec, (coords,)), shape=(mat_size,))
    assert b_sparse.ndim == 1

    ctx = Context()
    ctx.schur(range(2), sp.coo_array(a))
    x = ctx.solve_schur(b_sparse)

    assert x.ndim == 1
    assert_array_almost_equal(dtype, a @ x, bvec)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_factor_error(dtype):
    """Test that an error is raised if factor is asked to reuse missing analysis."""
    a = sp.identity(10, dtype=dtype)
    with pytest.raises(ValueError):
        Context().factor(a, reuse_analysis=True)

    Context().set_matrix(a)
    with pytest.raises(ValueError):
        Context().factor(reuse_analysis=True)


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
    sol = ctx.solve(np.array([1], dtype=dtype))
    if ctx.myid == 0:
        assert_almost_equal(dtype, sol, 1)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_zero_size_rhs(dtype):
    """Test that a 0xn rhs can be solved."""
    a = np.random.randn(10, 10).astype(dtype)
    ctx = Context()
    ctx.factor(a)
    rhs = np.zeros((10, 0), dtype=dtype)
    sol = ctx.solve(rhs)
    if ctx.myid == 0:
        assert_almost_equal(dtype, sol, rhs)


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
    sol = ctx.solve(rhs)
    if ctx.myid == 0:
        assert_almost_equal(dtype, sol, la.solve(a, rhs))


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [2, 10, 100], ids=str)
def test_slogdet_with_dense(dtype, mat_size):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    ctx = Context()
    sign, logabsdet = ctx.slogdet(sp.csr_matrix(a))
    # relative comparison of large numbers
    det = la.det(a)
    if ctx.myid == 0:
        assert_almost_equal(dtype, sign, det / np.abs(det))
        assert_almost_equal(dtype, logabsdet, np.log(np.abs(det)))

    # test singular matrix
    b = np.zeros((mat_size + 1, mat_size + 1), dtype)
    b[:mat_size][:, :mat_size] = a
    ctx = Context()
    sign, logabsdet = ctx.slogdet(sp.csr_matrix(b))
    if ctx.myid == 0:
        assert_almost_equal(dtype, sign, 0)
        assert_almost_equal(dtype, logabsdet, -np.inf)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", [2, 10, 100], ids=str)
def test_signature_with_dense(dtype, mat_size):
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    a += a.T.conj()

    sign_ref = np.sum(np.sign(la.eigvalsh(a)))

    if dtype in [np.complex64, np.complex128]:
        a = complex_to_real(sp.csr_array(a))
        sign_ref *= 2

    ctx = Context()
    if sp.issparse(a):
        sign = ctx.signature(a)
    else:
        sign = ctx.signature(sp.csr_array(a))
    if ctx.myid == 0:
        assert sign == sign_ref


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.mpi
def test_mpi_run(dtype):
    from mpi4py import MPI

    ctx = Context(comm=MPI.COMM_WORLD)
    S = sp.coo_array(
        (
            np.array([1.0, 2.0, 3.0, 4.0]),
            (np.array([0, 1, 2, 3]), np.array([0, 1, 2, 3])),
        ),
        dtype=dtype,
    )
    b = np.array([1.0, 2.0, 3.0, 4.0], dtype=dtype)
    ctx.set_matrix(S)
    ctx.analyze()
    ctx.factor()
    sol = ctx.solve(b)
    if ctx.myid == 0:
        assert np.allclose(sol, np.ones(4, dtype=dtype))
    else:
        assert sol is None


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.mpi
def test_mpi_subcomm_run(dtype):
    from mpi4py import MPI

    world = MPI.COMM_WORLD
    if world.Get_size() < 2:
        pytest.skip("Subcommunicator test requires at least 2 MPI ranks.")

    color = 0 if world.rank != 0 else MPI.UNDEFINED
    subcomm = world.Split(color=color, key=world.rank)

    if subcomm != MPI.COMM_NULL:
        ctx = Context(comm=subcomm)
        S = sp.coo_array(
            (
                np.array([1.0, 2.0, 3.0, 4.0]),
                (np.array([0, 1, 2, 3]), np.array([0, 1, 2, 3])),
            ),
            dtype=dtype,
        )
        b = np.array([1.0, 2.0, 3.0, 4.0], dtype=dtype)
        ctx.set_matrix(S)
        ctx.analyze()
        ctx.factor()
        sol = ctx.solve(b)
        if ctx.myid == 0:
            assert np.allclose(sol, np.ones(4, dtype=dtype))
        else:
            assert sol is None

    world.Barrier()
