# Copyright 2011-2016 Kwant authors.
# Copyright 2018 Python-MUMPS Authors.
#
# This file is part of Python-MUMPS. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# Python-MUMPS authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at
# https://gitlab.kwant-project.org/kwant/python-mumps.

"""Interface to the MUMPS sparse solver library"""

__all__ = [
    "Context",
    "schur_complement",
    "nullspace",
    "AnalysisStatistics",
    "FactorizationStatistics",
    "MUMPSError",
    "orderings",
]

import time
import numpy as np
import scipy.sparse
import scipy.linalg as la

from mumps import _mumps
from mumps.fortran_helpers import prepare_for_fortran

orderings = {
    "amd": 0,
    "amf": 2,
    "scotch": 3,
    "pord": 4,
    "metis": 5,
    "qamd": 6,
    "auto": 7,
}

ordering_name = [
    "amd",
    "user-defined",
    "amf",
    "scotch",
    "pord",
    "metis",
    "qamd",
]


def possible_orderings():
    """Return the ordering options that are available in the current
    installation of MUMPS.

    Which ordering options are actually available depends how MUMPS was
    compiled. Note that passing an ordering that is not available in the
    current installation of MUMPS will not fail, instead MUMPS will fall back
    to a supported one.

    Returns
    -------
    orderings : list of strings
       A list of installed orderings that can be used in the `ordering` option
       of MUMPS.
    """

    if not possible_orderings.cached:
        # Try all orderings on a small test matrix, and check which one was
        # actually used.

        possible_orderings.cached = ["auto"]
        for ordering in [0, 2, 3, 4, 5, 6]:
            data = np.asfortranarray([1, 1], dtype=np.complex128)
            row = np.asfortranarray([1, 2], dtype=_mumps.int_dtype)
            col = np.asfortranarray([1, 2], dtype=_mumps.int_dtype)

            instance = _mumps.zmumps()
            instance.set_assembled_matrix(2, row, col, data)
            instance.icntl[7] = ordering
            instance.job = 1
            instance.call()

            if instance.infog[7] == ordering:
                possible_orderings.cached.append(ordering_name[ordering])

    return possible_orderings.cached


possible_orderings.cached = None


error_messages = {
    -5: "Not enough memory during analysis phase",
    -6: "Matrix is singular in structure",
    -7: "Not enough memory during analysis phase",
    -10: "Matrix is numerically singular",
    -11: "The authors of MUMPS would like to hear about this",
    -12: "The authors of MUMPS would like to hear about this",
    -13: "Not enough memory",
    -37: "Nullspace basis solver not available for unsymmetric matrices",
}


class MUMPSError(RuntimeError):
    def __init__(self, infog):
        self.error = infog[1]
        if self.error in error_messages:
            msg = "{}. (MUMPS error {})".format(error_messages[self.error], self.error)
        else:
            msg = "MUMPS failed with error {}.".format(self.error)

        RuntimeError.__init__(self, msg)


class AnalysisStatistics:
    def __init__(self, inst, time=None):
        self.est_mem_incore = inst.infog[17]
        self.est_mem_ooc = inst.infog[27]
        self.est_nonzeros = (
            inst.infog[20] if inst.infog[20] > 0 else -inst.infog[20] * 1000000
        )
        self.est_flops = inst.rinfog[1]
        self.ordering = ordering_name[inst.infog[7]]
        self.time = time

    def __str__(self):
        parts = [
            "estimated memory for in-core factorization:",
            str(self.est_mem_incore),
            "mbytes\n",
            "estimated memory for out-of-core factorization:",
            str(self.est_mem_ooc),
            "mbytes\n",
            "estimated number of nonzeros in factors:",
            str(self.est_nonzeros),
            "\n",
            "estimated number of flops:",
            str(self.est_flops),
            "\n",
            "ordering used:",
            self.ordering,
        ]
        if hasattr(self, "time"):
            parts.extend(["\n analysis time:", str(self.time), "secs"])
        return " ".join(parts)


class FactorizationStatistics:
    def __init__(self, inst, time=None, include_ordering=False):
        # information about pivoting
        self.offdiag_pivots = inst.infog[12] if inst.sym == 0 else 0
        self.delayed_pivots = inst.infog[13]
        self.tiny_pivots = inst.infog[25]

        # possibly include ordering (used in schur_complement)
        if include_ordering:
            self.ordering = ordering_name[inst.infog[7]]

        # information about runtime efficiency
        self.memory = inst.infog[22]
        self.nonzeros = (
            inst.infog[29] if inst.infog[29] > 0 else -inst.infog[29] * 1000000
        )
        self.flops = inst.rinfog[3]
        if time:
            self.time = time

    def __str__(self):
        parts = [
            "off-diagonal pivots:",
            str(self.offdiag_pivots),
            "\n",
            "delayed pivots:",
            str(self.delayed_pivots),
            "\n",
            "tiny pivots:",
            str(self.tiny_pivots),
            "\n",
        ]
        if hasattr(self, "ordering"):
            parts.extend(["ordering used:", self.ordering, "\n"])
        parts.extend(
            [
                "memory used during factorization:",
                str(self.memory),
                "mbytes\n",
                "nonzeros in factored matrix:",
                str(self.nonzeros),
                "\n",
                "floating point operations:",
                str(self.flops),
            ]
        )
        if hasattr(self, "time"):
            parts.extend(["\n factorization time:", str(self.time), "secs"])
        return " ".join(parts)


class Context:
    """Context contains the internal data structures needed by the
    MUMPS library and contains a user-friendly interface.

    WARNING: Only complex numbers supported.

    Examples
    --------

    Solving a small system of equations.

    >>> import scipy.sparse as sp
    >>> a = sp.coo_array([[1., 0], [0, 2.]], dtype=complex)
    >>> ctx = mumps.Context()
    >>> ctx.factor(a)
    >>> ctx.solve([1., 1.])
    array([ 1.0+0.j,  0.5+0.j])

    Instance variables
    ------------------

    analysis_stats : `AnalysisStatistics`
        contains MUMPS statistics after an analysis step (i.e.  after a call to
        `analyze` or `factor`)
    factor_stats : `FactorizationStatistics`
        contains MUMPS statistics after a factorization step (i.e.  after a
        call to `factor`)

    """

    def __init__(self, verbose=False):
        """Init the Context class

        Parameters
        ----------

        verbose : True or False
            control whether MUMPS prints lots of internal statistics
            and debug information to screen.
        """
        self.mumps_instance = None
        self.dtype = None
        self.verbose = verbose
        self.factored = False
        self.schur_complement = None
        self.schur_indices = None
        self.schur_rhs = None
        self.schur_x = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # Force MUMPS to deallocate memory
        self.job = -2
        self.call()
        self.mumps_instance = None
        return False

    def call(self):
        """Execute the MUMPS subroutine.

        Compared to directly calling the MUMPS subroutine, this method
        automatically checks the return code and raises an exception if
        an error occurred. Additionally it returns the time spent.

        Raises
        ------
        MUMPSError
            if the MUMPS subroutine returned an error code.

        Returns
        -------
        time : float
            time spent in the MUMPS subroutine.
        """
        t1 = time.perf_counter()
        self.mumps_instance.call()
        t2 = time.perf_counter()
        if self.mumps_instance.infog[1] < 0:
            raise MUMPSError(self.mumps_instance.infog)
        return t2 - t1

    def set_matrix(self, a, overwrite_a=False, symmetric=False):
        """Set the matrix to be used in the next analysis or factorization step.

        Parameters
        ----------
        a : sparse SciPy array
            input matrix. Internally, the matrix is converted to `coo` format
            (so passing this format is best for performance).
        overwrite_a : True or False
            whether the data in a may be overwritten, which can lead to a small
            performance gain. Default is False.
        symmetric: True or False
            whether the matrix is symmetric. Default is False. If True, only
            the upper triangular part of the matrix is used.

        Note:
        -----
        On complex matrices ``symmetric=True`` means symmetric and not Hermitian.
        """
        a = scipy.sparse.coo_array(a)
        if symmetric:
            a = scipy.sparse.triu(a)

        if a.ndim != 2 or a.shape[0] != a.shape[1]:
            raise ValueError("Input matrix must be square!")

        dtype, row, col, data = _make_assembled_from_coo(a, overwrite_a)
        sym = 2 if symmetric else 0
        if self.dtype != dtype:
            self.mumps_instance = getattr(_mumps, dtype + "mumps")(self.verbose, sym)
            self.dtype = dtype
        # Note: We store the matrix data to avoid garbage collection.
        # See https://gitlab.kwant-project.org/kwant/python-mumps/-/issues/13
        self.n = a.shape[0]
        self.row = row
        self.col = col
        self.data = data
        self.mumps_instance.set_assembled_matrix(a.shape[0], row, col, data)

    def analyze(self, a=None, ordering="auto", overwrite_a=False):
        """Perform analysis step of MUMPS.

        In the analysis step, MUMPS figures out a reordering for the matrix and
        estimates number of operations and memory needed for the factorization
        time. This step usually needs not be called separately (it is done
        automatically by `factor`), but it can be useful to test which ordering
        would give best performance in the actual factorization, as MUMPS
        estimates are available in `analysis_stats`.

        Parameters
        ----------

        a : sparse SciPy matrix
            input matrix. Internally, the matrix is converted to `coo` format
            (so passing this format is best for performance). If `a` is not
            given, the matrix passed to `set_matrix` is used.
        ordering : { 'auto', 'amd', 'amf', 'scotch', 'pord', 'metis', 'qamd' }
            ordering to use in the factorization. The availability of a
            particular ordering depends on the MUMPS installation.  Default is
            'auto'.
        overwrite_a : True or False
            whether the data in a may be overwritten, which can lead to a small
            performance gain. Default is False.
        """
        if ordering not in orderings.keys():
            raise ValueError("Unknown ordering '" + ordering + "'!")

        if a is not None:
            self.set_matrix(a, overwrite_a)
        self.mumps_instance.icntl[7] = orderings[ordering]
        self.mumps_instance.job = 1
        t = self.call()
        self.factored = False
        self.analysis_stats = AnalysisStatistics(self.mumps_instance, t)

    def factor(
        self,
        a=None,
        ordering="auto",
        ooc=False,
        pivot_tol=0.01,
        reuse_analysis=False,
        overwrite_a=False,
    ):
        """Perform the LU factorization of the matrix.

        This LU factorization can then later be used to solve a linear system
        with `solve`. Statistical data of the factorization is stored in
        `factor_stats`.

        Parameters
        ----------

        a : sparse SciPy matrix
            input matrix. Internally, the matrix is converted to `coo` format
            (so passing this format is best for performance). If `a` is not
            given, the matrix passed to `analyze` is used.
        ordering : { 'auto', 'amd', 'amf', 'scotch', 'pord', 'metis', 'qamd' }
            ordering to use in the factorization. The availability of a
            particular ordering depends on the MUMPS installation.  Default is
            'auto'.
        ooc : True or False
            whether to use the out-of-core functionality of MUMPS.
            (out-of-core means that data is written to disk to reduce memory
            usage.) Default is False.
        pivot_tol: number in the range [0, 1]
            pivoting threshold. Pivoting is typically limited in sparse
            solvers, as too much pivoting destroys sparsity. 1.0 means full
            pivoting, whereas 0.0 means no pivoting. Default is 0.01.
        reuse_analysis: True or False
            whether to reuse the analysis done in a previous call to `analyze`
            or `factor`. If the structure of the matrix stays the same, and the
            numerical values do not change much, the previous analysis can be
            reused, saving some time.  WARNING: There is no check whether the
            structure of your matrix is compatible with the previous
            analysis. Also, if the values are not similar enough, there might
            be loss of accuracy, without a warning. Default is False.
        overwrite_a : True or False
            whether the data in a may be overwritten, which can lead to a small
            performance gain. Default is False.
        """
        if reuse_analysis and self.mumps_instance is None:
            raise ValueError("Missing analysis although reuse_analysis=True.")
        if a is not None:
            self.set_matrix(a, overwrite_a)
        if not reuse_analysis:
            self.analyze(ordering=ordering, overwrite_a=overwrite_a)

        self.mumps_instance.icntl[22] = 1 if ooc else 0
        self.mumps_instance.job = 2
        self.mumps_instance.cntl[1] = pivot_tol

        while True:
            try:
                t = self.call()
            except MUMPSError:
                # error -8, -9 (not enough allocated memory) is treated
                # specially, by increasing the memory relaxation parameter
                if self.mumps_instance.infog[1] in (-8, -9):
                    # Double the memory relaxation parameter
                    self.mumps_instance.icntl[14] *= 2
                else:
                    raise
            else:
                break

        self.factored = True
        self.factor_stats = FactorizationStatistics(self.mumps_instance, t)

    def _solve_sparse(self, b):
        b = b.tocsc()
        x = np.empty((b.shape[0], b.shape[1]), order="F", dtype=self.data.dtype)

        dtype, col_ptr, row_ind, data = _make_sparse_rhs_from_csc(b, self.data.dtype)

        if b.shape[0] != self.n:
            raise ValueError("Right hand side has wrong size")

        if self.dtype != dtype:
            raise ValueError(
                "Data type of right hand side is not "
                "compatible with the dtype of the "
                "linear system"
            )

        self.mumps_instance.set_sparse_rhs(col_ptr, row_ind, data)
        self.mumps_instance.set_dense_rhs(x)
        self.mumps_instance.job = 3
        self.mumps_instance.icntl[20] = 1
        self.call()

        return x

    def _solve_dense(self, b, overwrite_b=False):
        dtype, b = prepare_for_fortran(
            overwrite_b, b, np.zeros(1, dtype=self.data.dtype)
        )[:2]

        if b.shape[0] != self.n:
            raise ValueError("Right hand side has wrong size")

        if self.dtype != dtype:
            raise ValueError(
                "Data type of right hand side is not "
                "compatible with the dtype of the "
                "linear system"
            )

        self.mumps_instance.set_dense_rhs(b)
        self.mumps_instance.job = 3
        self.call()

        return b

    def solve(self, b, overwrite_b=False):
        """Solve a linear system after the LU factorization has previously
        been performed by `factor`.

        Supports both dense and sparse right hand sides.

        Parameters
        ----------

        b : dense (NumPy) matrix or vector or sparse (SciPy) matrix
            the right hand side to solve. Accepts both dense and sparse input;
            if the input is sparse 'csc' format is used internally (so passing
            a 'csc' matrix gives best performance).
        overwrite_b : True or False
            whether the data in b may be overwritten, which can lead to a small
            performance gain. Default is False.

        Returns
        -------

        x : NumPy array
            the solution to the linear system as a dense matrix (a vector is
            returned if b was a vector, otherwise a matrix is returned).
        """
        if b.ndim == 2 and b.shape[1] == 0:
            # Empty right hand side
            # We can return the copy directly because there is no data to be mutated.
            return b

        if not self.factored:
            raise RuntimeError("Factorization must be done before solving!")

        if scipy.sparse.isspmatrix(b):
            return self._solve_sparse(b)
        else:
            return self._solve_dense(b, overwrite_b)

    def schur(
        self,
        indices,
        a=None,
        ordering="auto",
        ooc=False,
        pivot_tol=0.01,
        overwrite_a=False,
        discard_factors=False,
    ):
        """Compute the Schur complement block of matrix a using MUMPS.

        Parameters:
        indices : 1d array
            indices (row and column) of the desired Schur complement block.  (The
            Schur complement block is square, so that the indices are both row and
            column indices.)
        a : sparse matrix
            input matrix. Internally, the matrix is converted to `coo` format (so
            passing this format is best for performance)
        ordering : { 'auto', 'amd', 'amf', 'scotch', 'pord', 'metis', 'qamd' }
            ordering to use in the factorization. The availability of a particular
            ordering depends on the MUMPS installation.  Default is 'auto'.
        ooc : True or False
            whether to use the out-of-core functionality of MUMPS.  (out-of-core
            means that data is written to disk to reduce memory usage.) Default is
            False.
        pivot_tol: number in the range [0, 1]
            pivoting threshold. Pivoting is typically limited in sparse solvers, as
            too much pivoting destroys sparsity. 1.0 means full pivoting, whereas
            0.0 means no pivoting. Default is 0.01.
        overwrite_a : True or False
            whether the data in a may be overwritten, which can lead to a small
            performance gain. Default is False.
        discard_factors: True or False
            whether to discard all matrix factors during factorization phase.
            Default is False.

        Returns
        -------

        schur_compl: NumPy array
            Schur complement block
        """
        if a is not None:
            self.set_matrix(a, overwrite_a)

        indices = np.asanyarray(indices)
        if indices.ndim != 1:
            raise ValueError("Schur indices must be specified in a 1d array!")
        self.schur_indices = indices = _makemumps_index_array(indices)
        self.schur_complement = np.empty(
            (indices.size, indices.size), order="C", dtype=self.data.dtype
        )

        self.mumps_instance.icntl[19] = 1
        self.mumps_instance.set_schur(self.schur_complement, indices)

        self.mumps_instance.icntl[31] = 1 if discard_factors else 0

        self.analyze(ordering=ordering)
        self.factor(reuse_analysis=True, ooc=ooc, pivot_tol=pivot_tol)

        return self.schur_complement

    def schur_condense(self, b, overwrite_b=False):
        """Return the condensed right hand side vector or matrix.

        This requires that the factorization has previously been performed by ``schur``.
        Supports both dense and sparse right hand sides.

        Parameters
        ----------

        b : dense (NumPy) matrix or vector or sparse (SciPy) matrix
            the right hand side to solve. Accepts both dense and sparse input;
            if the input is sparse 'csc' format is used internally (so passing
            a 'csc' matrix gives best performance).
        overwrite_b : True or False
            whether the data in b may be overwritten, which can lead to a small
            performance gain. Default is False.

        Returns
        -------

        schur_rhs : NumPy array
            the solution to the linear system as a dense matrix (a vector is
            returned if b was a vector, otherwise a matrix is returned).
        """

        if self.schur_complement is None:
            raise RuntimeError(
                "Factorization must be done by calling 'schur()' before solving!"
            )

        if b.shape[0] != self.n:
            raise ValueError("Right hand side has wrong size")

        if len(b.shape) > 1 and b.shape[0] > 0:
            raise ValueError("Right hand side must be a vector, not a matrix.")

        dtype = self.data.dtype

        self.schur_rhs = np.empty((self.schur_indices.size,), dtype=dtype)
        self.mumps_instance.set_schur_rhs(self.schur_rhs)

        if scipy.sparse.isspmatrix(b):
            b = b.tocsc()
            self.schur_x = np.empty((b.shape[0], b.shape[1]), order="F", dtype=dtype)
            b_dtype, col_ptr, row_ind, data = _make_sparse_rhs_from_csc(b, dtype)

            self.mumps_instance.set_sparse_rhs(col_ptr, row_ind, data)
            self.mumps_instance.set_dense_rhs(self.schur_x)
            self.mumps_instance.icntl[20] = 1

        else:
            b_dtype, b = prepare_for_fortran(overwrite_b, b, np.zeros(1, dtype=dtype))[
                :2
            ]
            self.mumps_instance.set_dense_rhs(b)
            self.schur_x = b

        if self.dtype != b_dtype:
            raise ValueError(
                "Data type of right hand side is not compatible with the dtype of the "
                "linear system"
            )

        self.mumps_instance.icntl[26] = 1  # Reduction/condensation phase
        self.mumps_instance.job = 3
        self.call()

        return self.schur_rhs

    def schur_expand(self, x2):
        """Perform the expansion step and return the complete solution.

        Parameters
        ----------
        x2 : NumPy array (vector or a matrix)
            the partial solution of the Schur system.

        Returns
        -------

        x : NumPy array
            the solution to the linear system as a dense matrix (a vector is
            returned if b was a vector, otherwise a matrix is returned).
        """
        if self.schur_rhs is None:
            raise RuntimeError(
                "Condensation must be done by calling 'schur_condense()' before expansion!"
            )

        self.schur_rhs[:] = x2
        self.mumps_instance.icntl[26] = 2  # Expansion phase
        self.mumps_instance.job = 3
        self.call()

        return self.schur_x

    def solve_schur(self, b, overwrite_b=False):
        """Solve a linear system using Schur complement method.

        This requires that the factorization has previously been performed by ``schur``.
        Supports both dense and sparse right hand sides.

        Parameters
        ----------

        b : dense (NumPy) matrix or vector or sparse (SciPy) matrix
            the right hand side to solve. Accepts both dense and sparse input;
            if the input is sparse 'csc' format is used internally (so passing
            a 'csc' matrix gives best performance).
        overwrite_b : True or False
            whether the data in b may be overwritten, which can lead to a small
            performance gain. Default is False.

        Returns
        -------

        x : NumPy array
            the solution to the linear system as a dense matrix (a vector is
            returned if b was a vector, otherwise a matrix is returned).

        Example
        -------

        Solving a system of equations.

        >>> import mumps
        >>> import numpy as np
        >>> import scipy.sparse as sp
        >>> row = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3])
        >>> col = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 1])
        >>> val = np.array([1, 2, 2, 1, 1, 3, -1, 2, 1, 1, 3, 1], dtype='d')
        >>> b = np.array([15, 12, 3, 5], dtype='d')
        >>> schur_indices = np.array([2, 3])
        >>> mtx = sp.coo_matrix((val, (row, col)), shape=(4, 4))
        >>> ctx = mumps.Context()
        >>> ctx.schur(schur_indices, mtx)
        >>> ctx.solve_schur(b)
        array([1., 2., 3., 4.])
        """
        schur_rhs = self.schur_condense(b, overwrite_b=overwrite_b)

        # solve dense system
        if self.mumps_instance.sym:
            # Schur matrix is lower triangular for symmetric matrices
            x2 = la.solve(self.schur_complement, schur_rhs, lower=True, assume_a="sym")
        else:
            x2 = la.solve(self.schur_complement, schur_rhs)

        x = self.schur_expand(x2)

        return x


def schur_complement(
    a,
    indices,
    ordering="auto",
    ooc=False,
    pivot_tol=0.01,
    calc_stats=False,
    overwrite_a=False,
):
    """Compute the Schur complement block of matrix a using MUMPS.

    Parameters:
    a : sparse matrix
        input matrix. Internally, the matrix is converted to `coo` format (so
        passing this format is best for performance)
    indices : 1d array
        indices (row and column) of the desired Schur complement block.  (The
        Schur complement block is square, so that the indices are both row and
        column indices.)
    ordering : { 'auto', 'amd', 'amf', 'scotch', 'pord', 'metis', 'qamd' }
        ordering to use in the factorization. The availability of a particular
        ordering depends on the MUMPS installation.  Default is 'auto'.
    ooc : True or False
        whether to use the out-of-core functionality of MUMPS.  (out-of-core
        means that data is written to disk to reduce memory usage.) Default is
        False.
    pivot_tol: number in the range [0, 1]
        pivoting threshold. Pivoting is typically limited in sparse solvers, as
        too much pivoting destroys sparsity. 1.0 means full pivoting, whereas
        0.0 means no pivoting. Default is 0.01.
    calc_stats: True or False
        whether to return the analysis and factorization statistics collected
        by MUMPS. Default is False.
    overwrite_a : True or False
        whether the data in a may be overwritten, which can lead to a small
        performance gain. Default is False.

    Returns
    -------

    s : NumPy array
        Schur complement block
    factor_stats: `FactorizationStatistics`
        statistics of the factorization as collected by MUMPS.  Only returned
        if ``calc_stats==True``.
    """
    with Context() as ctx:
        schur_compl = ctx.schur(
            indices,
            a,
            ordering=ordering,
            ooc=ooc,
            overwrite_a=overwrite_a,
            pivot_tol=pivot_tol,
            discard_factors=True,
        )

    if calc_stats:
        return [schur_compl, ctx.analysis_stats, ctx.factor_stats]
    else:
        return schur_compl


def nullspace(a, symmetric=False, pivot_threshold=0.0):
    """Compute the right nullspace using MUMPS.

    Parameters:
    a : scipy sparse matrix
        input matrix. Internally, the matrix is copied and converted
        to `coo` format.
    symmetric: bool, default: False
        If 'True', then 'a' must contain the diagonal and upper triangle
        of the (symmetric) matrix for which we want to find the nullspace.
        A ValueError is raised if this is not the case.
    pivot_threshold : double (-inf, inf)
        used to determine if a pivot of the preprocessed matrix is null.
        Controls the cntl(3) parameter. See the documentation for
        more information of the preprocessed matrix.

    Returns
    -------
    nullspace : 2D numpy array
                The columns vectors form part of the orthogonal basis for the
                right nullspace of the input matrix. The dimension of
                the row vectors are equal to the number of null pivots
                detected. It may be the complete null space basis.

    Notes
    -----
    Some versions of MUMPS do not support finding the nullspace of
    non-symmetric matrices. In this case a MUMPSError will be raised.
    """
    with Context() as ctx:
        ctx.set_matrix(a, symmetric=symmetric)

        ordering = "auto"
        ctx.mumps_instance.icntl[7] = orderings[ordering]
        ctx.mumps_instance.icntl[24] = 1
        ctx.mumps_instance.cntl[3] = pivot_threshold
        ctx.mumps_instance.job = 4

        # Find null pivots
        ctx.call()

        # Get the size of the null space
        n_null = ctx.mumps_instance.infog[28]
        if n_null == 0:
            return np.empty((a.shape[1], 0))

        # Initialize matrix for null space basis
        nullspace = np.zeros((a.shape[1], n_null), dtype=ctx.data.dtype, order="F")
        # Set RHS
        ctx.mumps_instance.set_dense_rhs(nullspace)
        ctx.mumps_instance.job = 3
        # Return all null space basis vectors, overwriting RHS
        ctx.mumps_instance.icntl[25] = -1
        ctx.call()

    # Orthonormalize
    nullspace, _ = la.qr(nullspace, mode="economic")

    return nullspace


# Some internal helper functions
def _make_assembled_from_coo(a, overwrite_a):
    """
    Takes the scipy sparse coo matrix and converts
    to the correct fortran arrays.

    Parameters:
    a : scipy sparse coo matrix
        input matrix.
    overwrite_a : bool
        controls the prepare_for_fortran function behaviour
        in fortran_helpers.py.

    Returns
    -------
    dtype : char 'SDCZ'
        a character indicating what fortran data type the
        data array contains.
        s - single precision real
        d - double precision real
        c - single precision complex
        z - double precision complex
    row, col : numpy integer fortran array
        the row and column indices of the non-zero
        elements in the input matrix.
    data : numpy fortran array (same type as dtype)
        array containing the non-zero value
        of the elements in the input matrix.
    """
    dtype, data = prepare_for_fortran(overwrite_a, a.data)

    row = np.asfortranarray(a.row.astype(_mumps.int_dtype))
    col = np.asfortranarray(a.col.astype(_mumps.int_dtype))

    # MUMPS uses Fortran indices.
    row += 1
    col += 1

    return dtype, row, col, data


def _make_sparse_rhs_from_csc(b, dtype):
    dtype, data = prepare_for_fortran(True, b.data, np.zeros(1, dtype=dtype))[:2]

    col_ptr = np.asfortranarray(b.indptr.astype(_mumps.int_dtype))
    row_ind = np.asfortranarray(b.indices.astype(_mumps.int_dtype))

    # MUMPS uses Fortran indices.
    col_ptr += 1
    row_ind += 1

    return dtype, col_ptr, row_ind, data


def _makemumps_index_array(a):
    a = np.asfortranarray(a.astype(_mumps.int_dtype))
    a += 1  # Fortran indices

    return a
