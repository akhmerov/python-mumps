# Copyright 2011-2016 Anton Akhmerov, Christoph Groth, and Michael Wimmer
# Copyright 2018 Mumpy Authors.
#
# This file is part of mumpy. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# mumpy authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at
# https://gitlab.kwant-project.org/kwant/mumpy.

from mumpy import MUMPSContext, schur_complement, MUMPSError
import pytest
import numpy as np
import scipy.sparse as sp
from ._test_utils import _Random, assert_array_almost_equal

"""
Here we define the datatypes and the expected decimal precision from the algorithms.
For the double precision datatypes, the IEEE specification only guarantees up to approximately
16 decimals digits, but with so many operations involving these numbers,
practically it is good for up to approximately 12, hence, we check for only 10 decimal places.
As for the single precision datatypes, since the specification allows for up to only
approximately 6 decimal digits, we only can ensure precision up to 1 decimal place.
"""
precisions = {
    np.float32: 1,
    np.float64: 10,
    np.complex64: 1,
    np.complex128: 10,
}

dtypes = list(precisions.keys())

matrix_sizes = np.arange(100,600,100)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", matrix_sizes, ids=str)
def test_lu_with_dense(dtype, mat_size):
    precision = precisions.get(dtype)
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    bmat = rand.randmat(mat_size, mat_size, dtype)
    bvec = rand.randvec(mat_size, dtype)

    ctx = MUMPSContext()
    ctx.factor(sp.coo_matrix(a))

    xvec = ctx.solve(bvec)
    xmat = ctx.solve(bmat)

    np.testing.assert_array_almost_equal(np.max(np.abs(np.imag(np.dot(a, xmat) - bmat))), 
                                         0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.real(np.dot(a, xmat) - bmat))), 
                                         0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.imag(np.dot(a, xvec) - bvec))),
                                         0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.real(np.dot(a, xvec) - bvec))),
                                         0., precision)

    # now "sparse" right hand side
    xvec = ctx.solve(sp.csc_matrix(bvec.reshape(mat_size,1)))
    xmat = ctx.solve(sp.csc_matrix(bmat))

    np.testing.assert_array_almost_equal(np.max(np.abs(np.imag(np.dot(a, xmat) - bmat))), 0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.real(np.dot(a, xmat) - bmat))), 0., precision)
    
    np.testing.assert_array_almost_equal(np.max(np.abs(np.imag(np.dot(a, xvec) - 
                                            bvec.reshape(mat_size,1)))), 0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.real(np.dot(a, xvec) - 
                                            bvec.reshape(mat_size,1)))), 0., precision)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
@pytest.mark.parametrize("mat_size", matrix_sizes, ids=str)
def test_schur_complement_with_dense(dtype, mat_size):
    precision = precisions.get(dtype)
    rand = _Random()
    a = rand.randmat(mat_size, mat_size, dtype)
    s = schur_complement(sp.coo_matrix(a), list(range(3)))
    np.testing.assert_array_almost_equal(np.max(np.abs(np.imag(np.linalg.inv(s) 
                                                - np.linalg.inv(a)[:3, :3]))),
                                            0., precision)
    np.testing.assert_array_almost_equal(np.max(np.abs(np.real(np.linalg.inv(s) 
                                                - np.linalg.inv(a)[:3, :3]))),
                                            0., precision)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_factor_warning(dtype):
    """Test that a warning is raised if factor is asked without analysis."""
    a = sp.identity(10, dtype=dtype)
    with pytest.warns(RuntimeWarning):
        MUMPSContext().factor(a, reuse_analysis=True)


@pytest.mark.parametrize("dtype", dtypes, ids=str)
def test_error_minus_9(dtype):
    """Test if MUMPSError -9 is properly caught by increasing memory"""
    ctx = MUMPSContext()

    rand = _Random()
    a = sp.coo_matrix(rand.randmat(1000, 1000, dtype))

    # Create the context so we can modify it for factorization
    ctx.analyze(a)

    # We ensure that this first call creates a -9 error by allocating only 1 MB for factorization
    ctx.mumps_instance.icntl[23] = 1 # This parameter allocates the maximum size of the working memory in MBytes per processor
    ctx.mumps_instance.icntl[14] = 1 # 
    ctx.mumps_instance.job = 2
    ctx.mumps_instance.call()
    assert(ctx.mumps_instance.infog[1] == -9)

    # This call should not raise any errors as it would successfully allocate memory
    MUMPSContext().factor(a)
