# This file is a part of python-mumps, copyright Python-mumps authors, BSD 2-Clause
# License. It also copies parts of the MUMPS source code, specifically the parameter
# descriptions, which is available under the license below.

# The license below and the provenance comments "Begin MUMPS snippet" are generated
# automatically by util/update_mumps_descriptions.py, and may not be edited manually.

# The script was last run using
# MUMPS VERSION: 5.8.1

# MUMPS LICENSE
#
#   Copyright 1991-2025 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
#   Mumps Technologies, University of Bordeaux.
#
#   This version of MUMPS is provided to you free of charge. It is
#   released under the CeCILL-C license
#   (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
#   https://cecill.info/licences/Licence_CeCILL-C_V1-en.html),
#   except for variants of AMD ordering and [sdcz]MUMPS_TRUNCATED_RRQR
#   derived from the LAPACK package distributed under BSD 3-clause
#   license (see headers of ana_orderings.F and [sdcz]lr_core.F),
#   and except for the external and optional ordering PORD provided
#   in a separate directory PORD (see PORD/README for License information).
#
#   You can acknowledge (using references [1] and [2]) the contribution
#   of this package in any scientific publication dependent upon the use
#   of the package. Please use reasonable endeavours to notify the authors
#   of the package of this publication.
#
#    [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
#    A fully asynchronous multifrontal solver using distributed dynamic
#    scheduling, SIAM Journal on Matrix Analysis and Applications,
#    Vol 23, No 1, pp 15-41 (2001).
#
#    [2] P. R. Amestoy, A. Buttari, J.-Y. L'Excellent and T. Mary,
#    Performance and scalability of the block low-rank multifrontal
#    factorization on multicore architectures,
#    ACM Transactions on Mathematical Software,
#    Vol 45, Issue 1, pp 2:1-2:26 (2019)
#
#   As a counterpart to the access to the source code and rights to copy,
#   modify and redistribute granted by the license, users are provided only
#   with a limited warranty  and the software's author,  the holder of the
#   economic rights,  and the successive licensors  have only  limited
#   liability.
#
#   In this respect, the user's attention is drawn to the risks associated
#   with loading,  using,  modifying and/or developing or reproducing the
#   software by the user in light of its specific status of free software,
#   that may mean  that it is complicated to manipulate,  and  that  also
#   therefore means  that it is reserved for developers  and  experienced
#   professionals having in-depth computer knowledge. Users are therefore
#   encouraged to load and test the software's suitability as regards their
#   requirements in conditions enabling the security of their systems
#   and/or data to be ensured and, more generally, to use and operate it
#   in the same conditions as regards security.
#
#   The fact that you are presently reading this means that you have had
#   knowledge of the CeCILL-C license and that you accept its terms.
#

# END MUMPS LICENSE

"""Enums for MUMPS control/info arrays.

This module provides typed wrappers for the MUMPS parameter arrays.

Scope: - Arrays: CNTL, ICNTL, INFO, RINFO, RINFOG.
"""

from __future__ import annotations

from ._param_base import (
    ParamArray,
    RawArray,
    _OneBasedArray,
    param,
)


# Canonical lengths (from ZMUMPSID signature and headers)
LEN_CNTL = 15
LEN_ICNTL = 60
LEN_INFO = 80
LEN_RINFO = 40
LEN_RINFOG = 40

# -------------
# ICNTL members
# -------------


class ICNTL(ParamArray):
    """Integer control parameters ICNTL(1..60).

    Exposes selected parameters with tokens; raw 1-based indexing is also available.
    """

    def __init__(self, size: int = LEN_ICNTL):
        super().__init__(size=size)

    # === Begin MUMPS snippet: ICNTL(1) page 72 from userguide_5.8.1.txt:4015-4019 ===
    # ICNTL(1) is the output stream for error messages.
    #       Possible values :
    #     ≤ 0: these messages will be suppressed.
    #     > 0 : is the output stream.
    #       Default value: 6 (standard output stream)
    # === End MUMPS snippet ===

    @param(index=1, page=72)
    class error_stream:
        """
        Output stream for error messages.

        Default: `standard_output`.

        Setting to `suppressed` disables error message output; otherwise a
        selected stream receives error messages.
        """

        suppressed = 0, "suppress error messages"
        standard_error = 5, "standard error stream"
        standard_output = 6, "standard output stream"

    # === Begin MUMPS snippet: ICNTL(2) page 72 from userguide_5.8.1.txt:4021-4027 ===
    # ICNTL(2) is the output stream for diagnostic printing and statistics local to each MPI process.
    #       Possible values :
    #     ≤ 0: these messages will be suppressed.
    #     > 0 : is the output stream.
    #       Default value: 0
    #       Remarks: If ICNTL(2) > 0 and ICNTL(4) ≥ 2, then information on advancement (flops done)
    #       is also printed.
    # === End MUMPS snippet ===

    @param(index=2, page=72)
    class local_diagnostic_stream:
        """
        Per-process diagnostic and statistics output stream.

        Default: `suppressed`.

        When not `suppressed` and `message_level` is at least `errors_and_warnings`,
        additional advancement information (e.g. flops done) will also be printed.
        """

        suppressed = 0, "suppress per-process diagnostics"
        standard_error = 5, "standard error stream"
        standard_output = 6, "standard per-process output stream"

    # === Begin MUMPS snippet: ICNTL(3) page 72 from userguide_5.8.1.txt:4029-4033 ===
    # ICNTL(3) is the output stream for global information, collected on the host.
    #       Possible values :
    #     ≤ 0: these messages will be suppressed.
    #     > 0 : is the output stream.
    #       Default value: 6 (standard output stream)
    # === End MUMPS snippet ===

    @param(index=3, page=72)
    class host_diagnostic_stream:
        """
        Host-level diagnostic and statistics output stream.

        Default: `standard_output`.
        """

        suppressed = 0, "suppress host diagnostics"
        standard_error = 5, "standard error stream"
        standard_output = 6, "standard host output stream"

    # === Begin MUMPS snippet: ICNTL(4) page 72 from userguide_5.8.1.txt:4035-4042 ===
    # ICNTL(4) is the level of printing for error, warning, and diagnostic messages.
    #       Possible values :
    #     ≤ 0: No messages output.
    #       1 : Only error messages printed.
    #       2 : Errors, warnings, and main statistics printed.
    #       3 : Errors and warnings and terse diagnostics (only first ten entries of arrays) printed.
    #     ≥ 4 : Errors, warnings and information on input, output parameters printed.
    #       Default value: 2 (errors, warnings and main statistics printed)
    # === End MUMPS snippet ===

    @param(index=4, page=72)
    class message_level:
        """
        Printing/verbosity level

        Default: Errors and warnings are printed (`errors_and_warnings`)
        """

        no_messages = 0, "no messages"
        errors_only = 1, "only error messages printed"
        errors_and_warnings = 2, "errors and warnings printed"
        terse_diagnostics = 3, "errors, warnings, and terse diagnostics"
        verbose_all = 4, "all information on input and output parameters printed"

    # === Begin MUMPS snippet: ICNTL(5) page 72 from userguide_5.8.1.txt:4044-4069 ===
    # ICNTL(5) controls the matrix input format (see Subsection 5.4.2).
    #       Phase: accessed by the host and only during the analysis phase
    #       Possible variables/arrays involved: N, NNZ (or NZ for backward compatibility), IRN, JCN,
    #       NNZ loc (or NZ loc for backward compatibility), IRN loc, JCN loc, A loc, NELT, ELTPTR,
    #       ELTVAR, and A ELT
    #       Possible values :
    #        0 : assembled format. The matrix must be input in the structure components N, NNZ (or NZ),
    #            IRN, JCN, and A if the matrix is centralized on the host (see Subsection 5.4.2.1) or in the
    #            structure components N, NNZ loc (or NZ loc), IRN loc, JCN loc, A loc if the matrix is
    #            distributed (see Subsection 5.4.2.2).
    #        1 : elemental format. The matrix must be input in the structure components N, NELT, ELTPTR,
    #            ELTVAR, and A ELT (see Subsection 5.4.2.3).
    #      Any other values will be treated as 0.
    #      Default value: 0 (assembled format)
    #      Related parameters: ICNTL(18)
    #      Incompatibility: If the matrix is in elemental format (ICNTL(5)=1), the BLR feature
    #      (ICNTL(35)≥ 1) is currently not available, see error -800.
    #      Remarks: NNZ and NNZ loc are 64-bit integers (NZ and NZ loc are 32-bit integers kept for
    #      backward compatibility and will be obsolete in future releases).
    #      Parallel analysis (ICNTL(28) =2) is only available for matrices in assembled format and, thus, an
    #      error will be raised for elemental matrices (ICNTL(5)=1).
    #      Elemental matrices can be input only centralized on the host (ICNTL(18)=0).
    # === End MUMPS snippet ===

    @param(index=5, page=72)
    class input_format:
        """
        Matrix input format.

        Default: `assembled`.

        For `assembled`, the matrix may be centralized on the host or distributed;
        for `element`, input uses the elemental data structures and must be
        centralized on the host (`ICNTL(18)=0`). Parallel analysis (`ICNTL(28)`) is
        available only with `assembled`. The BLR feature (`ICNTL(35)`) is not
        available with `element`. Values outside the defined members are treated as
        `assembled`. Related: `ICNTL(18)`.
        """

        assembled = 0, "assembled matrix"
        element = 1, "elemental matrix"

    # === Begin MUMPS snippet: ICNTL(6) page 73 from userguide_5.8.1.txt:4071-4099 ===
    # ICNTL(6) permutes the matrix to a zero-free diagonal and/or scale the matrix (see Subsection 3.2 and
    # Subsection 5.5.2).
    #      Phase: accessed by the host and only during sequential analysis (ICNTL(28)=1)
    #      Possible variables/arrays involved: optionally UNS PERM, mumps par%A, COLSCA and ROWSCA
    #      Possible values :
    #        0 : No column permutation is computed.
    #        1 : The permuted matrix has as many entries on its diagonal as possible. The values on the
    #            diagonal are of arbitrary size.
    #        2 : The permutation is such that the smallest value on the diagonal of the permuted matrix is
    #            maximized. The numerical values of the original matrix, (mumps par%A), must be provided
    #            by the user during the analysis phase.
    #        3 : Variant of option 2 with different performance. The numerical values of the original matrix
    #            (mumps par%A) must be provided by the user during the analysis phase.
    #        4 : The sum of the diagonal entries of the permuted matrix is maximized. The numerical values of
    #            the original matrix (mumps par%A) must be provided by the user during the analysis phase.
    #        5 : The product of the diagonal entries of the permuted matrix is maximized. Scaling vectors
    #            are also computed and stored in COLSCA and ROWSCA, if ICNTL(8) is set to -2 or 77.
    #            With these scaling vectors, the nonzero diagonal entries in the permuted matrix are one in
    #            absolute value and all the off-diagonal entries less than or equal to one in absolute value.
    #            For unsymmetric matrices, COLSCA and ROWSCA are meaningful on the permuted matrix
    #            A Qc (see Equation (5)). For symmetric matrices, COLSCA and ROWSCA are meaningful on
    #            the original matrix A. The numerical values of the original matrix, mumps par%A, must be
    #            provided by the user during the analysis phase.
    #        6 : Similar to 5 but with a more costly (time and memory footprint) algorithm. The numerical
    #            values of the original matrix, mumps par%A, must be provided by the user during the analysis
    #            phase.
    #        7 : Based on the structural symmetry of the input matrix and on the availability of the numerical
    #            values, the value of ICNTL(6) is automatically chosen by the software.
    #      Other values are treated as 0. On output from the analysis phase, INFOG(23) holds the value of
    # === End MUMPS snippet ===

    @param(index=6, page=73)
    class colperm_strategy:
        """
        Column permutation and optional scaling strategy used during analysis.

        Used only in sequential analysis (`ICNTL(28)=1`). Some strategies require
        access to numerical values during analysis. When scaling is requested
        (`ICNTL(8)` analysis-driven or automatic), certain strategies may compute
        and store scaling vectors in `COLSCA` and `ROWSCA`. Values outside the
        defined members are treated as `no_permutation`.
        """

        no_permutation = 0, "no column permutation"
        maximize_diagonal_nonzeros = 1, "maximize diagonal nonzeros"
        maximize_min_diagonal = (
            2,
            "maximize smallest diagonal entry (needs values during analysis)",
        )
        maximin_variant = (
            3,
            "variant of previous strategy (needs values during analysis)",
        )
        maximize_diagonal_sum = (
            4,
            "maximize sum of diagonal entries (needs values during analysis)",
        )
        maximize_diagonal_product = (
            5,
            "maximize product of diagonal entries; may compute scaling vectors",
        )
        maximize_diagonal_product_expensive = (
            6,
            "like previous but more costly; may compute scaling vectors",
        )
        automatic = 7, "automatic choice based on structure/values"

    # === Begin MUMPS snippet: ICNTL(7) page 74 from userguide_5.8.1.txt:4127-4185 ===
    #      ICNTL(7) computes a symmetric permutation (ordering) to determine the pivot order to be used for the
    #      factorization in case of sequential analysis (ICNTL(28)=1). See Subsection 3.2 and Subsection 5.6.
    #            Phase: accessed by the host and only during the sequential analysis phase (ICNTL(28) = 1).
    #            Possible variables/arrays involved: PERM IN, SYM PERM
    #            Possible values :
    #              0 : Approximate Minimum Degree (AMD) [6] is used,
    #              1 : The pivot order should be set by the user in PERM IN, on the host processor. In that case,
    #                  PERM IN must be allocated on the host by the user and PERM IN(i), (i=1, ... N) must hold the
    #                  position of variable i in the pivot order. In other words, row/column i in the original matrix
    #                  corresponds to row/column PERM IN(i) in the reordered matrix.
    #              2 : Approximate Minimum Fill (AMF) is used,
    #              3 : SCOTCH9 [42] package is used if previously installed by the user otherwise treated as 7.
    #              4 : PORD10 [46] is used if previously installed by the user otherwise treated as 7.
    #              5 : the Metis11 [34] package is used if previously installed by the user otherwise treated as 7.
    #                  It is possible to modify some components of the internal options array of Metis (see
    #                  Metis manual) in order to fine-tune and modify various aspects of the internal algorithms
    #                  used by Metis. This can be done by setting some elements (see the file metis.h in the
    #                  Metis installation to check the position of each option in the array) of the MUMPS array
    #                  mumps par%METIS OPTIONS after the MUMPS initialization phase (JOB= –1) and before
    #                  the analysis phase. Note that the METIS OPTIONS array of the MUMPS structure is of size 40,
    #                  which is large enough for both Metis 4.x and Metis 5.x versions. It is passed by MUMPS as the
    #                  argument "options" to the METIS ordering routine METIS NodeAND (METIS NodeWND is
    #                  sometimes also called in case MUMPS was installed with Metis 4.x) during the analysis phase.
    #              6 : Approximate Minimum Degree with automatic quasi-dense row detection (QAMD) is used.
    #              7 : Automatic choice by the software during analysis phase. This choice will depend on the
    #                  ordering packages made available, on the matrix (type and size), and on the number of
    #                  processors.
    #            Other values are treated as 7.
    #            Default value: 7 (automatic choice)
    #            Incompatibility: ICNTL(7) is meaningless if the parallel analysis is chosen (ICNTL(28)=2).
    #            Related parameters: ICNTL(28)
    # 9 See http://gforge.inria.fr/projects/scotch/ to obtain a copy.
    # 10 Distributed within MUMPS by permission of J. Schulze (University of Paderborn).
    # 11 See http://glaros.dtc.umn.edu/gkhome/metis/metis/overview to obtain a copy.
    #      Remarks: Even when the ordering is provided by the user, the analysis must be performed before
    #      numerical factorization.
    #      For assembled matrices (centralized or distributed) (ICNTL(5)=0) all the options are available.
    #      For elemental matrices (ICNTL(5)=1), only options 0, 1, 5 and 7 are available, with option 7
    #      leading to an automatic choice between AMD and Metis (options 0 or 5); other values are treated
    #      as 7.
    #      If the user asks for a Schur complement matrix (ICNTL(19)= 1, 2, 3) and
    #         – the matrix is assembled (ICNTL(5)=0) then only options 0, 1, 5 and 7 are currently available.
    #           Other options are treated as 7.
    #         – the matrix is elemental (ICNTL(5)=1) only options 0, 1 and 7 are currently available. Other
    #           options are treated as 7 which will (currently) be treated as 0 (AMD).
    #         – in both cases (assembled or elemental matrix) if the pivot order is given by the user
    #           (ICNTL(7)=1) then the following property should hold: PERM IN(LISTVAR SCHUR(i)) =
    #           N-SIZE SCHUR+i, for i=1,SIZE SCHUR.
    #      For matrices with relatively dense rows, we highly recommend option 6 which may significantly
    #      reduce the time for analysis.
    #      On output, the pointer array SYM PERM provides access, on the host processor, to the symmetric
    #      permutation that is effectively computed during the analysis phase by the MUMPS package, and
    #      INFOG(7) to the ordering option that was effectively chosen. In fact, the option corresponding to
    # === End MUMPS snippet ===

    @param(index=7, page=74)
    class ordering:
        """Ordering strategy for pivoting.

        Default is automatic choice during analysis. If the user wishes
        to force the pivot order provided in the IS array, set this parameter
        to 1 before analysis. The effective ordering actually used is stored
        in KEEP(256). Several orderings require external packages (SCOTCH,
        PORD, Metis) and will only be available if those packages are installed.
        For values other than the ones enumerated below, a suitable pivot order
        will be chosen automatically.
        """

        amd = 0, "approximate minimum degree (AMD; included in MUMPS)"
        user = 1, "ordering provided by the user (use pivot order in IS)"
        amf = 2, "approximate minimum fill (AMF; included in MUMPS)"
        scotch = 3, "SCOTCH ordering (requires external SCOTCH package)"
        pord = 4, "PORD ordering (Juergen Schulze's PORD; external)"
        metis = 5, "Metis ordering (requires external Metis package)"
        amd_auto_dense = (
            6,
            "AMD with automatic quasi-dense row detection (use when AMD ordering time is large)",
        )
        automatic = 7, "automatic choice during analysis"

    # === Begin MUMPS snippet: ICNTL(8) page 75 from userguide_5.8.1.txt:4192-4238 ===
    # ICNTL(8) describes the scaling strategy (see Subsection 5.5).
    #      Phase: accessed by the host during analysis phase (that need be sequential ICNTL(28)=1) or on
    #      entry to numerical factorization phase
    #      Possible variables/arrays involved: COLSCA, ROWSCA
    #      Possible values :
    #       -2: Scaling done during analysis (see [23, 24] for the unsymmetric case and [25] for the symmetric
    #           case), under certain conditions. The original matrix has to be centralized (ICNTL(18)=0)
    #           and the user has to provide its numerical values (mumps par%A) on entry to the analysis,
    #           otherwise the scaling will not be computed. Also, ICNTL(6) must be activated for the scaling
    #           to be computed.
    #           The effective value of the scaling applied is reported in INFOG(33). We recommend that the
    #           user checks for INFOG(33) in order to know if the scaling was applied. When not applied,
    #           other scalings can be performed during the factorization.
    #       -1: Scaling provided by the user. Scaling arrays must be provided in COLSCA and ROWSCA on
    #           entry to the numerical factorization phase by the user, who is then responsible for allocating
    #           and freeing them. If the input matrix is symmetric (SYM= 1 or 2), then the user should ensure
    #           that the array ROWSCA is equal to (or points to the same location as) the array COLSCA.
    #       0 : No scaling applied/computed.
    #       1 : Diagonal scaling computed during the numerical factorization phase,
    #       3 : Column scaling computed during the numerical factorization phase,
    #       4 : Row and column scaling based on infinite row/column norms, computed during the numerical
    #           factorization phase,
    #       7 : Simultaneous row and column iterative scaling (based on [45, 15, 36, 35]) computed during
    #           the numerical factorization phase.
    #       8 : Similar to 7 but more rigorous and expensive to compute; computed during the numerical
    #           factorization phase.
    #       77 : Automatic choice of the value of ICNTL(8) done during analysis.
    #       Other values are treated as 77.
    #       Default value: 77 (automatic choice done by the package)
    #       Related parameters: ICNTL(6), ICNTL(12)
    #       Remarks: If ICNTL(8) = 77, then an automatic choice of the scaling option may be performed,
    #       either during the analysis or the factorization. The effective value used for ICNTL(8) is returned in
    #       INFOG(33). If the scaling arrays are computed during the analysis, then they are ready to be used
    #       by the factorization phase. Note that scalings can be efficiently computed during analysis when
    #       requested (see ICNTL(6) and ICNTL(12)).
    #       If the input matrix is real and symmetric with SYM= 1 then automatic choice is no scaling. However,
    #       the user may want to scale the matrix when BLR feature is activated (see ICNTL(35)).
    #       Incompatibility: If the input matrix is symmetric (SYM= 1 or 2), then only options –2, –1, 0, 1, 7,
    #       8 and 77 are allowed and other options are treated as 0. If the input matrix is in elemental format
    #       (ICNTL(5) = 1), then only options –1 and 0 are allowed and other options are treated as 0. If
    #       the input matrix is assembled and distributed (ICNTL(18)=1,2,3 and ICNTL(5) = 0), then only
    #       options 7, 8 and 77 are allowed, otherwise no scaling is applied.
    #       If block format is exploited (ICNTL(15)̸= 0) then scaling is not applied.
    # === End MUMPS snippet ===

    @param(index=8, page=75)
    class scaling:
        """
        Scaling strategy.

        Default: `automatic` (effective choice reported in `INFOG(33)`). Interacts
        with `ICNTL(6)` and `ICNTL(12)`. Some options are restricted depending on
        symmetry (`SYM`), input format (`ICNTL(5)`), and distribution (`ICNTL(18)`).
        When block format is used (`ICNTL(15) != 0`), scaling is not applied. If
        computed during analysis, the scaling arrays are ready to be used by the
        factorization phase. Values outside the defined members are treated as
        `automatic`.
        """

        analysis = (
            -2,
            "compute during analysis (centralized matrix; needs numerical values; ICNTL(6) active)",
        )
        user_provided = -1, "provided by user in COLSCA/ROWSCA"
        none = 0, "no scaling"
        diagonal = 1, "diagonal scaling during factorization"
        column = 3, "column scaling during factorization"
        row_col_norm = 4, "row and column norm-based scaling during factorization"
        iterative = 7, "simultaneous row/column iterative scaling during factorization"
        iterative_rigorous = 8, "more rigorous and expensive iterative scaling"
        automatic = 77, "automatic choice"

    # === Begin MUMPS snippet: ICNTL(9) page 76 from userguide_5.8.1.txt:4240-4248 ===
    # ICNTL(9) computes the solution using A or AT
    #       Phase: accessed by the host during the solve phase.
    #       Possible values :
    #        1 : AX = B is solved.
    #     ̸= 1 : AT X = B is solved.
    #       Default value: 1
    #       Related parameters: ICNTL(10), ICNTL(11), ICNTL(21), ICNTL(32)
    #       Remarks: when a forward elimination is performed during the factorization (see ICNTL(32))
    #       only ICNTL(9)=1 is allowed.
    # === End MUMPS snippet ===

    @param(index=9, page=76)
    class transpose_solve:
        """
        Solve form: normal vs transpose.

        Default: `normal`. When forward elimination during factorization is
        enabled (`ICNTL(32)`), only `normal` is allowed. Related: `ICNTL(10)`,
        `ICNTL(11)`, `ICNTL(21)`, `ICNTL(32)`.
        """

        normal = 1, "solve A x = b"
        transpose = 0, "solve A^T x = b"

    # === Begin MUMPS snippet: ICNTL(10) page 76 from userguide_5.8.1.txt:4250-4283 ===
    # ICNTL(10) applies the iterative refinement to the computed solution (see Subsection 5.8).
    #       Phase: accessed by the host during the solve phase.
    #       Possible variables/arrays involved: NRHS
    #       Possible values :
    #     < 0 : Fixed number of steps of iterative refinement. No stopping criterion is used.
    #       0 : No iterative refinement.
    #     > 0 : Maximum number of steps of iterative refinement. A stopping criterion is used, therefore a
    #           test for convergence is done at each step of the iterative refinement algorithm.
    #       Default value: 0 (no iterative refinement)
    #       Related parameters: CNTL(2)
    #       Incompatibility: If ICNTL(21)=1 (solution kept distributed) or if ICNTL(32)=1 (forward
    #       elimination during factorization), or if NRHS>1 (multiple right hand sides), or if ICNTL(20)=10
    #       or 11 (distributed right hand sides), then iterative refinement is disabled and ICNTL(10) is treated
    #       as 0.
    #       Remarks: Note that if ICNTL(10)< 0, |ICNTL(10)| steps of iterative refinement are performed,
    #       without any test of convergence (see Algorithm 3). This means that the iterative refinement may
    #       diverge, that is the solution instead of being improved may be worse from an accuracy point of view.
    #       But it has been shown [19] that with only two to three steps of iterative refinement the solution can
    #       often be significantly improved. So if the convergence test should not be done we recommend to
    #       set ICNTL(10) to -2 or -3.
    #      Note also that it is not necessary to activate the error analysis option (ICNTL(11)= 1,2) to be
    #      able to run the iterative refinement with stopping criterion (ICNTL(10) > 0). However, since
    #      the backward errors ω1 and ω2 have been computed, they are still returned in RINFOG(7) and
    #      RINFOG(8), respectively.
    #      It must be noticed that iterative refinement with stopping criterion (ICNTL(10) > 0) will stop
    #      when
    #        1. either the requested accuracy is reached (ω1 + ω2 < CNTL(2))
    #        2. or when the convergence rate is too slow (ω1 + ω2 does not decrease by at least a factor of 2)
    #        3. or when exactly ICNTL(10) steps have been performed.
    #      In the first two cases the number of iterative refinement steps (INFOG(15)) may be lower than
    # === End MUMPS snippet ===

    @param(index=10, page=76)
    class iterative_refinement:
        """
        Number of iterative refinement steps.

        Default: no iterative refinement. Negative values perform a fixed number of
        steps without a stopping criterion; positive values set a maximum number of
        steps with a stopping criterion (accuracy target and convergence rate), and
        the process may terminate earlier. Disabled when the solution is kept
        distributed, when forward elimination is used during factorization, when
        multiple right-hand sides are present, or with specific distributed RHS
        modes. Related: `CNTL(2)`. Backward error measures are returned regardless
        of whether full error analysis is requested (`ICNTL(11)`).
        """

    # === Begin MUMPS snippet: ICNTL(11) page 77 from userguide_5.8.1.txt:4286-4314 ===
    # ICNTL(11) computes statistics related to an error analysis of the linear system solved (Ax = b or
    # AT x = b (see ICNTL(9))). See Subsection 5.9.
    #      Phase: accessed by the host and only during the solve phase.
    #      Possible variables/arrays involved: NRHS
    #      Possible values :
    #       0 : no error analysis is performed (no statistics).
    #       1 : compute all the statistics (very expensive).
    #       2 : compute main statistics (norms, residuals, componentwise backward errors), but not the most
    #           expensive ones like (condition number and forward error estimates).
    #      Values different from 0, 1, and 2 are treated as 0.
    #      Default value: 0 (no statistics).
    #      Incompatibility: If ICNTL(21)=1 (solution kept distributed) or if ICNTL(32)=1 (forward
    #      elimination during factorization), or if NRHS>1 (multiple right hand sides), or if ICNTL(20)=10
    #      or 11 (distributed right hand sides), or if ICNTL(25)=-1 (computation of the null space basis),
    #      then error analysis is not performed and ICNTL(11) is treated as 0.
    #      Related parameters: ICNTL(9)
    #      Remarks: The computed statistics are returned in various informational parameters:
    #         – If ICNTL(11)= 2, then the infinite norm of the input matrix (∥A∥∞ or ∥AT ∥∞ in
    #           RINFOG(4)), the infinite norm of the computed solution (∥x̄∥∞ in RINFOG(5)), and the
    #                             ∥Ax̄−b∥∞
    #           scaled residual ∥A∥ ∞ ∥x̄∥∞
    #                                       in RINFOG(6), a componentwise backward error estimate in
    #           RINFOG(7) and RINFOG(8) are computed.
    #         – If ICNTL(11)= 1, then in addition to the above statistics also an estimate for the error in
    #           the solution in RINFOG(9), and condition numbers for the linear system in RINFOG(10) and
    #           RINFOG(11) are also returned.
    #      If performance is critical, ICNTL(11) should be set to 0. If both performance is critical and statistics
    #      are requested, then ICNTL(11) should be set to 2. If ICNTL(11)=1, the error analysis is very costly
    #      (typically significantly more costly than the solve phase itself).
    # === End MUMPS snippet ===

    @param(index=11, page=77)
    class error_analysis:
        """
        Error analysis statistics during solve.

        Default: `none`. When enabled, statistics are computed and returned in
        informational arrays. This feature is disabled when the solution is kept
        distributed, when forward elimination during factorization is enabled,
        when multiple right-hand sides are present, when distributed right-hand
        sides are used, or when the null space basis is requested. Related:
        `transpose_solve`.
        """

        none = 0, "no error analysis"
        all_stats = 1, "compute all statistics (expensive)"
        main_stats = 2, "compute main statistics (not the most expensive ones)"

    # === Begin MUMPS snippet: ICNTL(12) page 77 from userguide_5.8.1.txt:4316-4334 ===
    # ICNTL(12) defines an ordering strategy for symmetric matrices (SYM = 2) (see [25] for more details)
    # and is used, in conjunction with ICNTL(6), to add constraints to the ordering algorithm (ICNTL(7)
    # option).
    #      Phase: accessed by the host and only during the analysis phase.
    #      Possible values :
    #       0 : automatic choice
    #       1 : usual ordering (nothing done)
    #       2 : ordering on the compressed graph associated with the matrix.
    #       3 : constrained ordering, only available with AMF (ICNTL(7)=2).
    #      Other values are treated as 1.
    #      Default value: 0 (automatic choice).
    #      Incompatibility: If the matrix is unsymmetric (SYM=0) or symmetric definite positive matrices
    #      (SYM=1), or the matrix is in elemental format (ICNTL(5)=1), or the matrix is initially distributed
    #      (ICNTL(18)=1,2,3) or the ordering is provided by the user (ICNTL(7)=1), or the Schur option
    #      (ICNTL(19) ̸= 0) is required, or the analysis is performed by blocks (ICNTL(15) ̸= 0),
    # === End MUMPS snippet ===

    @param(index=12, page=77)
    class ldlt_strategy:
        """
        Ordering strategy for symmetric matrices used with LDLT variants.

        Default: `automatic`. Used during analysis and combined with
        `colperm_strategy` to add constraints to the ordering. Certain settings
        are not meaningful or allowed depending on symmetry, input format,
        initial distribution, user-provided ordering, Schur option, or block
        compression settings.
        """

        automatic = 0, "automatic choice"
        usual = 1, "usual ordering (no additional constraint)"
        compressed_graph = 2, "ordering on compressed graph"
        constrained_amf = 3, "constrained ordering (only with AMF)"

    # === Begin MUMPS snippet: ICNTL(13) page 78 from userguide_5.8.1.txt:4349-4374 ===
    # ICNTL(13) controls the parallelism of the root node (enabling or not the use of ScaLAPACK) and also
    # its splitting.
    #      Phase: accessed by the host during the analysis phase.
    #      Possible values :
    #    < -1 : treated as 0.
    #      -1 : force splitting of the root node in all cases (even sequentially)
    #       0 : parallel factorization of the root node based on ScaLAPACK. If the size of the root frontal node
    #           (last Schur complement to be factored) is larger than an internal threshold, then ScaLAPACK
    #           will be used for factorizing it. Otherwise, the root node will be processed by a single MPI
    #           process.
    #    > 0 : ScaLAPACK is not used (recommended value is 1 to partly recover parallelism of the root
    #           node). It forces a sequential factorization of the root node (ScaLAPACK will not be used). To
    #           recover parallelism lost by the fact of not using ScaLAPACK, splitting of the root node can
    #           be activated: if the number of working processes is strictly larger than ICNTL(13) (always
    #           the case with ICNTL(13)=1) then splitting of the root node is performed to enable node level
    #           parallelism.
    #      Default value: 0 (parallel factorization on the root node)
    #      Remarks: Processing the root sequentially (ICNTL(13) > 0) can be useful when the user is
    #      interested in the inertia of the matrix (see INFO(12) and INFOG(12)), or when the user wants
    #      to detect null pivots (see Subsection 5.12) or to activate BLR compression (Subsection 5.19) on the
    #      root node.
    #      Although ICNTL(13) controls the efficiency of the factorization and solve phases, preprocessing
    #      work is performed during analysis and this option must be set on entry to the analysis phase.
    #      With SYM=1, if ScaLAPACK is allowed (ICNTL(13)≤ 0) then Cholesky factorization will be
    #      performed on the root node and thus negative pivots will raise an error (code -40 is returned in
    #      INFOG(1)).
    # === End MUMPS snippet ===

    @param(index=13, page=78)
    class scalapack_root_split:
        """
        Control ScaLAPACK use and splitting at the root node.

        Default: ScaLAPACK is allowed; large root fronts are factorized in
        parallel, otherwise the root is processed by a single process. You can
        force splitting of the root (even in sequential contexts), or disable
        ScaLAPACK entirely and process the root sequentially while optionally
        enabling splitting to regain some parallelism when enough processes are
        available. Sequential root processing is useful for inertia inspection,
        detecting null pivots, or enabling BLR compression on the root. Set this
        parameter before analysis.

        Value meanings:
        - -1: Force splitting of the root in all cases (including sequential runs).
        -  0: Allow ScaLAPACK at the root; if the root is large enough it will be
            factorized in parallel, otherwise a single process handles it.
        - >0: Disable ScaLAPACK (force sequential root). If the number of working
            processes is strictly larger than this value (commonly 1), splitting
            of the root is performed to regain node-level parallelism.
        - < -1: Treated as 0.
        """

    # === Begin MUMPS snippet: ICNTL(14) page 83 from userguide_5.8.1.txt:4601-4619 ===
    # ICNTL(23). The relaxation ICNTL(14) is first applied to the internal integer workarray IS and to
    # communication and I/O buffers; the remaining available space is then shared between the main (and
    # often most critical) real/complex internal workarray S holding the factors, the stack of contribution
    # blocks and dynamic workarrays that are used either to expand the S array or to store low-rank
    # dynamic structures.
    # Lower bounds for ICNTL(23), in case ICNTL(23) is provided only on the host:
    #    – In case of full-rank factors only (ICNTL(35)=0 or 3), a lower bound for ICNTL(23) (if ICNTL(14),
    #      has not been modified since the analysis) is given by INFOG(16) if the factorization is in-core
    #      (ICNTL(22)=0), and by INFOG(26) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank factors (ICNTL(35)=1 or 2) only (ICNTL(37)=0), a lower bound for ICNTL(23)
    #      (if ICNTL(14), has not been modified since the analysis and ICNTL(38) is a good approximation
    #      of the average compression rate of the factors) is given by INFOG(36) if the factorization is in-core
    #      (ICNTL(22)=0), and by INFOG(38) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank contribution blocks (CB) only (ICNTL(35)=0,3 and ICNTL(37)=1), a lower bound
    #      for ICNTL(23) (if ICNTL(14), has not been modified since the analysis and ICNTL(39) is a good
    #      approximation of the average compression rate of the CB) is given by INFOG(44) if the factorization is
    #      in-core (ICNTL(22)=0), and by INFOG(46) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank factors and contribution blocks (ICNTL(35)=1,2 and ICNTL(37)=1), a lower
    #      bound for ICNTL(23) (if ICNTL(14), has not been modified since the analysis, and ICNTL(38) and
    # === End MUMPS snippet ===

    @param(index=14, page=83)
    class mem_relaxation:
        """
        Memory relaxation used when allocating internal work arrays.

        Applied first to the internal integer work array (IS) and to communication/I/O
        buffers. The remaining space is then shared between: (a) the main
        real/complex work array that holds factors, (b) the stack of contribution
        blocks, and (c) dynamic work arrays used to expand that space or store
        low-rank dynamic structures.

        Works in tandem with `max_memory_mb` (`ICNTL(23)`), which caps the
        available working memory. See the manual for lower-bound guidance on
        `ICNTL(23)` based on INFO/INFOG values.
        """

    # === Begin MUMPS snippet: ICNTL(15) page 79 from userguide_5.8.1.txt:4388-4412 ===
    # ICNTL(15) exploits compression of the input matrix resulting from a block format, see Subsection 5.7.
    #      Phase: accessed by the host process during the analysis phase.
    #      Possible variables/arrays involved: NBLK, BLKPTR, BLKVAR
    #      Possible values :
    #        0: no compression
    #       -k: all blocks are of fixed size k> 0. N (the order of the matrix A) must be a multiple of k. NBLK
    #           and BLKPTR should not be provided by the user and will be computed internally. Concerning
    #           BLKVAR, please refer to the Remarks below.
    #        1: block format provided by the user. NBLK need be provided on the host by the user and
    #           holds the number of blocks. BLKPTR(1:NBLK+1) must be provided by the user on the host.
    #           Concerning BLKVAR, please refer to the Remarks below.
    #      Any other values will be treated as 0.
    #
    #      Default value: 0
    #      Remarks: If BLKVAR is not provided by the user then BLKVAR is internally treated as the identity
    #      (BLKVAR(i)=i, (i=1, ..., N)). It corresponds to contiguous variables in blocks.
    #         – If ICNTL(15)=1 then BLKVAR(BLKPTR(iblk):BLKPTR(iblk+1)-1), (iblk=1, NBLK) holds
    #           the variables associated to block iblk.
    #         – If ICNTL(15) < 0 then BLKPTR need not be provided by the user and NBLK = N/k where
    #           N must be a multiple of k.
    #      In case the pivot order is provided on entry by the user at the analysis phase (ICNTL(7)= 1) then
    #      PERM IN should be compatible with the compression. This means that PERM IN, of size N, should
    #      result from an expansion of a pivot order on the compressed matrix, i.e., variables in a block should
    #      be consecutive in the pivot order.
    #      Incompatibility: With element entry format ICNTL(5)= 1, with Schur complement
    # === End MUMPS snippet ===

    @param(index=15, page=79)
    class graph_compression:
        """
        Enable use of block (compressed) input format when the matrix has a block structure.

        Default: no compression. Two modes exist:
        - Fixed-size blocks: a negative value denotes blocks of uniform size k>0
            (N must be a multiple of k). The number of blocks and block pointers are
            computed internally; variable indices per block default to contiguous
            ranges if not provided.
        - User-provided block format: a positive value enables user-specified
            blocks. Provide the number of blocks and block pointers on the host;
            variable indices per block are given in the corresponding array.

        When a pivot order is provided by the user at analysis, it must be
        compatible with the block compression (variables of each block remain
        consecutive in the pivot order). Not available with elemental input
        format (`input_format` set to `element`).
        """

    # === Begin MUMPS snippet: ICNTL(16) page 79 from userguide_5.8.1.txt:4416-4427 ===
    # ICNTL(16) controls the setting of the number of OpenMP threads, see Subsection 5.21 by MUMPS when
    # the setting of multithreading is not possible outside MUMPS (see Subsection 3.13).
    #      Phase: accessed by the host at the beginning of all phases
    #      Possible values :
    #       0 : (recommended value) nothing is done, MUMPS uses the number of OpenMP threads configured
    #           by the calling application.
    #     > 0 : on entry MUMPS sets the number of OpenMP threads to the value of ICNTL(16) and on exit
    #           reset the number of OpenMP threads to its value before the MUMPS call.
    #      Other values are treated as 0.
    #      Default value: 0 (nothing is done)
    #      Remarks: Please note that the use of ICNTL(16) should be limited to the case when for some
    #      reasons, setting the number of OpenMP threads is not possible outside MUMPS.
    # === End MUMPS snippet ===

    @param(index=16, page=79)
    class omp_threads:
        """
        Number of OpenMP threads to use.

        Default: do not change threading; the application/environment controls
        the thread count. Setting a positive value instructs the solver to set
        the OpenMP thread count on entry and restore the previous value on exit.
        Use only when thread control cannot be managed outside the solver.
        """

    # === Begin MUMPS snippet: ICNTL(18) page 79 from userguide_5.8.1.txt:4429-4451 ===
    # ICNTL(18) defines the strategy for the distributed input matrix (only for assembled matrix, see
    # Subsection 5.4.2).
    #      Phase: accessed by the host during the analysis phase.
    #      Possible values :
    #       0 : the input matrix is centralized on the host (see Subsection 5.4.2.1).
    #       1 : the user provides the structure of the matrix on the host at analysis, MUMPS returns a mapping
    #           and the user should then provide the matrix entries distributed according to the mapping on
    #           entry to the numerical factorization phase (see Subsection 5.4.2.2).
    #       2 : the user provides the structure of the matrix on the host at analysis, and the distributed
    #           matrix entries on all slave processors at factorization. Any distribution is allowed (see
    #           Subsection 5.4.2.2).
    #       3 : user directly provides the distributed matrix, pattern and entries, input both for analysis and
    #           factorization (see Subsection 5.4.2.2).
    #      Other values are treated as 0.
    #      Default value: 0 (input matrix centralized on the host)
    #      Related parameters: ICNTL(5)
    #      Remarks: In case of distributed matrix, we recommend options 2 or 3. Among them, we
    #      recommend option 3 which is easier to use. Option 1 is kept for backward compatibility but is
    #      deprecated and we plan to suppress it in a future release.
    # === End MUMPS snippet ===

    @param(index=18, page=79)
    class distributed_input_strategy:
        """
        Strategy for distributed input of assembled matrices.

        Default: `centralized_on_host`.

        Related: `input_format`.

        Recommendations: prefer `host_structure_distributed_entries` or
        `fully_distributed`; `host_maps_then_redistribute` is kept for backward
        compatibility and may be removed in a future release.
        """

        centralized_on_host = 0, "input matrix centralized on host"
        host_maps_then_redistribute = (
            1,
            "host provides structure; mapping returned; redistribute at factorization",
        )
        host_structure_distributed_entries = (
            2,
            "host structure at analysis; distributed entries at factorization",
        )
        fully_distributed = (
            3,
            "pattern and entries provided distributed for analysis and factorization",
        )

    # === Begin MUMPS snippet: ICNTL(19) page 80 from userguide_5.8.1.txt:4453-4487 ===
    # ICNTL(19) computes the Schur complement matrix (see Subsection 5.18).
    #      Phase: accessed by the host during the analysis phase.
    #      Possible variables/arrays involved: SIZE SCHUR, LISTVAR SCHUR, NPROW, NPCOL, MBLOCK,
    #      NBLOCK, SCHUR, SCHUR MLOC, SCHUR NLOC, and SCHUR LLD
    #      Possible values :
    #        0 : complete factorization. No Schur complement is returned.
    #        1 : the Schur complement matrix will be returned centralized by rows on the host after the
    #            factorization phase. On the host before the analysis phase, the user must set the integer variable
    #            SIZE SCHUR to the size of the Schur matrix, the integer pointer array LISTVAR SCHUR to
    #            the list of indices of the Schur matrix.
    #   2 or 3 : the Schur complement matrix will be returned distributed by columns: the Schur will
    #            be returned on the slave processors in the form of a 2D block cyclic distributed matrix
    #            (ScaLAPACK style) after factorization. Workspace should be allocated by the user before
    #            the factorization phase in order for MUMPS to store the Schur complement (see SCHUR,
    #            SCHUR MLOC, SCHUR NLOC, and SCHUR LLD in Subsection 5.18). On the host before the
    #            analysis phase, the user must set the integer variable SIZE SCHUR to the size of the Schur
    #            matrix, the integer pointer array LISTVAR SCHUR to the list of indices of the Schur matrix.
    #            The integer variables NPROW, NPCOL, MBLOCK, NBLOCK may also be defined (default values
    #            will otherwise be provided).
    #      Values not equal to 1, 2 or 3 are treated as 0.
    #      Default value: 0 (complete factorization)
    #      Incompatibility: Since the Schur complement is a partial factorization of the global matrix (with
    #      partial ordering of the variables provided by the user), the following options of MUMPS are
    #      incompatible with the Schur option: rank-revealing factorization (ICNTL(56)=1), maximum
    #      transversal, scaling, iterative refinement, error analysis and parallel analysis.
    #      Related parameters: ICNTL(7), ICNTL(26)
    #      Remarks: If the ordering is given (ICNTL(7)= 1) then the following property should hold:
    #      PERM IN(LISTVAR SCHUR(i)) = N-SIZE SCHUR+i, for i=1,SIZE SCHUR.
    #      Note that, in order to have a centralized Schur complement matrix by columns (see
    #      Subsection 5.18.3), it is possible (and recommended) to use a particular case of the distributed
    #      Schur complement (ICNTL(19)=2 or 3), where the Schur complement is assigned to only one
    #      processor (NPCOL × NPROW = 1).
    # === End MUMPS snippet ===

    @param(index=19, page=80)
    class schur_option:
        """
        Schur complement extraction.

        Default: `complete_factorization` (no Schur complement). When enabled,
        a Schur complement is returned either centralized on the host or
        distributed in a 2D block-cyclic format. Incompatible with rank-revealing
        factorization, maximum transversal, scaling, iterative refinement, error
        analysis, and parallel analysis. Related: `ordering`.
        """

        complete_factorization = 0, "no Schur complement returned"
        centralized_by_rows = 1, "return centralized Schur on host (by rows)"
        distributed = 2, "return distributed Schur (2D block-cyclic)"

    # === Begin MUMPS snippet: ICNTL(20) page 81 from userguide_5.8.1.txt:4489-4528 ===
    # ICNTL(20) determines the format (dense, sparse, or distributed) of the right-hand sides
    #       Phase: accessed by the host during the solve phase and before a JOB= 9 call.
    #       Possible variables/arrays involved: RHS, NRHS, LRHS, IRHS SPARSE, RHS SPARSE,
    #       IRHS PTR, NZ RHS, Nloc RHS, LRHS loc, IRHS loc, RHS loc, INFO(23).
    #       Possible values :
    #        0 : the right-hand side is in dense format in the structure component RHS, NRHS, LRHS (see
    #            Subsection 5.17.1)
    #    1,2,3 : the right-hand side is in sparse format in the structure components IRHS SPARSE,
    #            RHS SPARSE, IRHS PTR and NZ RHS.
    #             1 : The decision of exploiting sparsity of the right-hand side to accelerate the solution phase
    #                 is done automatically.
    #             2 : Sparsity of the right-hand side is NOT exploited to improve solution phase.
    #             3 : Sparsity of the right-hand side is exploited during solution phase.
    #   10, 11 : distributed right-hand side.
    #            The right-hand side is provided distributed in the structure components Nloc RHS,
    #            LRHS loc, IRHS loc, RHS loc (see Subsection 5.17.3).
    #            When provided before a JOB= 9 call, values 10 and 11 indicate which distribution MUMPS
    #            should build and return it to the user in IRHS loc. In this case, the user should provide a
    #            workarray IRHS loc on each MPI process of size at least INFO(23), where INFO(23) is
    #            returned after the factorization phase.
    #            10 : fill IRHS loc to minimize internal communications of right-hand side data during the
    #                 solve phase.
    #            11 : fill IRHS loc to match the distribution of the solution in case of the distribution of the
    #                 solution is imposed by MUMPS (ICNTL(21)=1).
    #            In any case, values 10 and 11 have the same meaning entering the solve phase: use the
    #            distributed right-hand side in IRHS loc, RHS loc without taking into account how the
    #            distribution has been computed.
    #       Values different from 0, 1, 2, 3, 10, 11 are treated as 0. For a sparse right-hand side, the
    #       recommended value is 1.
    #       Default value: 0 (dense right-hand sides)
    #       Incompatibility: When NRHS > 1 (multiple right-hand side), the functionalities related to iterative
    #       refinement ( ICNTL(10)) and error analysis (ICNTL(11)) are currently disabled.
    #       With sparse right-hand sides (ICNTL(20)=1,2,3), the forward elimination during the factorization
    #       (ICNTL(32)=1) is not currently available.
    #       Remarks: For details on how to set the input parameters, see Subsection 5.17.1, Subsection 5.17.2
    #       and Subsection 5.17.3. Please note that duplicate entries in the sparse or distributed right-hand sides
    #       are summed. A JOB= 9 call can only be done after a successful factorization phase and its results
    #       depends on the transpose option ICNTL(9), which should not be modified between a JOB= 9 and
    #       JOB= 3 call. The distributed right-hand side feature enables the user to provide a sparse structured
    #       RHS (i.e., a RHS with some empty rows that will be considered equal to zero).
    # === End MUMPS snippet ===

    @param(index=20, page=81)
    class rhs_format:
        """
        Right-hand side format for the solve phase.

        Default: `dense`. Supports dense, sparse (with automatic detection, no
        exploitation, or explicit exploitation of sparsity), and distributed RHS
        modes that minimize communications or match the imposed solution
        distribution. Related: `distributed_solution`.

        When multiple right-hand sides are present, iterative refinement and error
        analysis are disabled. With sparse RHS modes, forward elimination during
        factorization is not available.
        """

        dense = 0, "dense RHS in RHS/NRHS/LRHS"
        sparse_auto = 1, "sparse RHS; automatic exploitation"
        sparse_no_exploit = 2, "sparse RHS; do not exploit sparsity"
        sparse_exploit = 3, "sparse RHS; explicitly exploit sparsity"
        distributed_minimize_comm = 10, "distributed RHS; minimize communications"
        distributed_match_solution = (
            11,
            "distributed RHS; match solution distribution if imposed",
        )

    # === Begin MUMPS snippet: ICNTL(21) page 81 from userguide_5.8.1.txt:4530-4550 ===
    # ICNTL(21) determines the distribution (centralized or distributed) of the solution vectors.
    #       Phase: accessed by the host during the solve phase.
    #       Possible variables/arrays involved: RHS, ISOL loc, SOL loc, LSOL loc, INFO(23).
    #       Possible values :
    #        0 : the solution vector is assembled and stored in the structure component RHS (gather phase),
    #            that must have been allocated earlier by the user (see Subsection 5.17.5).
    #        1 : the solution vector is kept distributed on each slave processor in the structure components
    #            ISOL loc and SOL loc. The distribution of the solution is chosen by MUMPS and is returned
    #            in ISOL loc. The arrays ISOL loc and SOL loc must have been allocated by the user and
    #            must be of size at least INFO(23), where INFO(23) has been returned by MUMPS at the
    #            end of the factorization phase (see Subsection 5.17.6).
    #      Values different from 0 and 1 are currently treated as 0.
    #      Default value: 0 (assembled centralized format)
    #      Incompatibility: If the solution is kept distributed, error analysis and iterative refinement (controlled
    #      by ICNTL(10) and ICNTL(11)) are not applied.
    # === End MUMPS snippet ===

    @param(index=21, page=81)
    class distributed_solution:
        """
        Solution distribution after solve.

        Default: `assembled` (gathered on the host). When `distributed`, the
        solution is kept distributed across processes. In distributed mode, error
        analysis and iterative refinement are not applied. Related arrays must be
        allocated according to the distribution chosen by the solver.
        """

        assembled = 0, "assemble solution on host in RHS"
        distributed = 1, "keep solution distributed (ISOL_loc/SOL_loc)"

    # === Begin MUMPS snippet: ICNTL(22) page 82 from userguide_5.8.1.txt:4552-4583 ===
    # ICNTL(22) controls the in-core/out-of-core (OOC) factorization and solve.
    #      Phase: accessed by the host during the factorization phase.
    #      Possible variables/arrays involved: OOC TMPDIR and OOC PREFIX
    #      Possible values :
    #        0: In-core factorization and solution phases (default standard version).
    #        1: Out-of-core factorization and solve phases. The complete matrix of factors is written to disk
    #           (see Subsection 3.15).
    #      Other values are treated as 0 (in-core factorization).
    #      Default value: 0 (in-core factorization)
    #      Related parameters: ICNTL(35) since factors in low-rank form are not written to disk
    #      Remarks: The variables OOC TMPDIR and OOC PREFIX are used to indicate the directory and
    #      the prefix, respectively, where to store the files containing the factors. They must be set after the
    #      initialization phase (JOB= –1) and before the factorization phase (JOB= 2,4,5 or 6). Otherwise,
    #      MUMPS will use the /tmp directory and arbitrary file names. Note MUMPS accesses to the variables
    #      OOC TMPDIR and OOC PREFIX only during the factorization phase. Several files under the same
    #      directory and with the same prefix are created to store the factors. Their names contain a unique
    #      hash and MUMPS is in charge of keeping trace of them.
    #      The files containing the factors will be deleted if a new factorization starts or when a deallocation
    #      phase (JOB= –2 or JOB= –4) is called, except if the save/restore feature has been used and the files
    #      containing the factors are associated to a saved instance. See Section Subsection 5.20.4).
    #      Note that, in case of abnormal termination of an application calling MUMPS (for example, a
    #      termination of the calling process with a segmentation fault, or, more generally, a termination of
    #      the calling process without a call to MUMPS with JOB= –2), the files containing the factors are not
    #      deleted. It is then the user’s responsibility to delete them, as shown in bold in the example below,
    #      where the application calling MUMPS is launched from a bash script and environment variables are
    #      used to define the OOC environment:
    #      #!/bin/bash
    #      export MUMPS OOC TMPDIR="/local/mumps data/"
    #      export MUMPS OOC PREFIX="job myapp "
    #      mpirun -np 128 ./myapplication
    #      # Suppress MUMPS OOC files in case of bad application termination
    #      rm -f ${MUMPS OOC TMPDIR}/${MUMPS OOC PREFIX}*
    # === End MUMPS snippet ===

    @param(index=22, page=82)
    class ooc_setting:
        """
        In-core vs out-of-core (OOC) factorization/solve.

        Default: `in_core`. In OOC mode, the factor matrices are written to disk
        using the configured temporary directory and file prefix, which must be
        set between initialization and factorization. Related: BLR setting since
        low-rank factors are not written to disk.
        """

        in_core = 0, "keep factors in memory (standard mode)"
        out_of_core = 1, "write factors to disk (OOC)"

    # === Begin MUMPS snippet: ICNTL(23) page 83 from userguide_5.8.1.txt:4601-4619 ===
    # ICNTL(23). The relaxation ICNTL(14) is first applied to the internal integer workarray IS and to
    # communication and I/O buffers; the remaining available space is then shared between the main (and
    # often most critical) real/complex internal workarray S holding the factors, the stack of contribution
    # blocks and dynamic workarrays that are used either to expand the S array or to store low-rank
    # dynamic structures.
    # Lower bounds for ICNTL(23), in case ICNTL(23) is provided only on the host:
    #    – In case of full-rank factors only (ICNTL(35)=0 or 3), a lower bound for ICNTL(23) (if ICNTL(14),
    #      has not been modified since the analysis) is given by INFOG(16) if the factorization is in-core
    #      (ICNTL(22)=0), and by INFOG(26) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank factors (ICNTL(35)=1 or 2) only (ICNTL(37)=0), a lower bound for ICNTL(23)
    #      (if ICNTL(14), has not been modified since the analysis and ICNTL(38) is a good approximation
    #      of the average compression rate of the factors) is given by INFOG(36) if the factorization is in-core
    #      (ICNTL(22)=0), and by INFOG(38) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank contribution blocks (CB) only (ICNTL(35)=0,3 and ICNTL(37)=1), a lower bound
    #      for ICNTL(23) (if ICNTL(14), has not been modified since the analysis and ICNTL(39) is a good
    #      approximation of the average compression rate of the CB) is given by INFOG(44) if the factorization is
    #      in-core (ICNTL(22)=0), and by INFOG(46) if the factorization is out-of-core (ICNTL(22)=1).
    #    – In case of low-rank factors and contribution blocks (ICNTL(35)=1,2 and ICNTL(37)=1), a lower
    #      bound for ICNTL(23) (if ICNTL(14), has not been modified since the analysis, and ICNTL(38) and
    # === End MUMPS snippet ===

    @param(index=23, page=83)
    class max_memory_mb:
        """
        Maximum allowed working memory per process in megabytes.

        If provided only on the host, the same limit is used on other processes.
        Lower bounds can be estimated from INFO/INFOG depending on the factor
        format and in-core/out-of-core mode. When unset, the solver allocates a
        workspace based on analysis estimates and relaxation (`mem_relaxation`).
        """

    # === Begin MUMPS snippet: ICNTL(24) page 83 from userguide_5.8.1.txt:4644-4664 ===
    # ICNTL(24) controls the detection of “null pivot rows”.
    #      Phase: accessed by the host during the factorization phase
    #      Possible variables/arrays involved: PIVNUL LIST
    #      Possible values :
    #        0: Nothing done. A null pivot row will result in error INFO(1)=-10.
    #        1: Null pivot row detection.
    #       Other values are treated as 0.
    #       Default value: 0 (no null pivot row detection)
    #       Related parameters: CNTL(3), CNTL(5), ICNTL(13), ICNTL(25), ICNTL(56)
    #       Remarks: CNTL(3) is used to compute the threshold to decide if a pivot row is “null”.
    #       It can be used alone or combined with the rank-revealing option (see ICNTL(56)). In this case
    #       CNTL(3) it is used to determine if a pivot should be postponed.
    #       Note that when ScaLAPACK is applied on the root node (see ICNTL(13) = 0), then exact null
    #       pivots on the root will stop the factorization (INFO(1)=-10) while if tiny pivots are present on
    #       the root node the ScaLAPACK routine will factorize the root matrix. Computing the root node
    #       factorization sequentially (this can be forced by setting ICNTL(13) to 1) will help with the correct
    #       detection of null pivots but may degrade performance.
    # === End MUMPS snippet ===

    @param(index=24, page=83)
    class null_pivot_detection:
        """
        Detection of null pivot rows during factorization.

        Default: `off`. When enabled, pivots deemed null trigger specific
        handling. Interacts with thresholds and fixation settings. Consider the
        impact on root processing when ScaLAPACK is used.
        """

        off = 0, "no null pivot detection (null pivots raise an error)"
        on = 1, "enable null pivot detection"

    # === Begin MUMPS snippet: ICNTL(25) page 84 from userguide_5.8.1.txt:4666-4690 ===
    # ICNTL(25) allows the computation of a solution of a deficient matrix and also of a null space basis.
    #       Phase: accessed by the host during the solution phase
    #       Possible variables/arrays involved: RHS, PIVNUL LIST, ISOL loc and SOL loc
    #       Possible values :
    #         0: A normal solution step is performed. If the matrix was found singular during factorization
    #            then one of the possible solutions is returned.
    #         i: with 1 ≤ i ≤ INFOG(28). The i-th vector of the null space basis is computed.
    #        -1: The complete null space basis is computed.
    #       Default value: 0 (normal solution step)
    #       Incompatibility: Iterative refinement, error analysis, and the option to solve the transpose system
    #       (ICNTL(9) ̸= 1) are ignored when the solution step is used to return vectors from the null space
    #       (ICNTL(25) ̸= 0).
    #       Related parameters: ICNTL(56), ICNTL(21), ICNTL(24)
    #       Remarks: Null space basis computation can be activated when a zero-pivot detection option
    #       (ICNTL(24) ̸= 0) and/or SVD decomposition on the root node (ICNTL(56) = 1) was requested
    #       during the factorization and the matrix was found to be deficient (INFOG(28) > 0).
    #       Note that when vectors from the null space are requested (ICNTL(25) ̸= 0), both centralized
    #       (ICNTL(21)=0) and distributed (ICNTL(21)=1) solutions options can be used. If the solution
    #       is centralized (ICNTL(21)=0), then the null space vectors are returned to the user in the array
    #       RHS, allocated by the user on the host. If the solution is distributed (ICNTL(21)=1), then the null
    #       space vectors are returned in the array SOL loc, which must be allocated by the user on all working
    #       processes (see Subsection 5.17.6). In both cases, the number of columns of RHS or SOL loc must
    #       be equal to the number of vectors requested, so that NRHS must be equal to:
    #         – 1 if 1 ≤ ICNTL(25) ≤ INFOG(28)
    #         – INFOG(28) if ICNTL(25)=-1.
    # === End MUMPS snippet ===

    @param(index=25, page=84)
    class null_space_solution:
        "null space solution control during solve (0 normal; i compute i-th vector; -1 compute all)"

    # === Begin MUMPS snippet: ICNTL(26) page 87 from userguide_5.8.1.txt:4861-4892 ===
    # ICNTL(26)=1
    # Possible values :
    #   0: standard factorization not involving right-hand sides.
    #   1: forward elimination (Equation (3)) of the right-hand side vectors is performed during
    #      factorization (JOB= 2). The solve phase (JOB= 3) will then only involve backward
    #      substitution (Equation (4)).
    # Other values are treated as 0.
    # Default value: 0 (standard factorization)
    # Related parameters: ICNTL(31),ICNTL(26)
    # Incompatibility: This option is incompatible with sparse right-hand sides (ICNTL(20)=1,2,3),
    # with the solution of the transposed system (ICNTL(9) ̸= 1), with the computation of entries of
    # the inverse (ICNTL(30)=1), and with BLR factorizations (ICNTL(35)=1,2,3). In such cases,
    # error -43 is raised.
    # Furthermore, iterative refinement (ICNTL(10)) and error analysis (ICNTL(11)) are disabled.
    # Finally, the current implementation imposes that all right-hand sides are processed in one pass
    # during the backward step. Therefore, the blocking size (ICNTL(27)) is ignored.
    # Remarks: The right-hand sides must be dense to use this functionality: RHS, NRHS, and LRHS
    # should be provided as described in Subsection 5.17.1. They should be provided at the beginning of
    # the factorization phase (JOB= 2) rather than at the beginning of the solve phase (JOB= 3).
    # For unsymmetric matrices, if the forward elimination is performed during factorization
    # (ICNTL(32) = 1), the L factor (see ICNTL(31)) may be discarded to save space. In fact,
    # the L factor will then always be discarded (even when ICNTL(31)=0) in the case of a full-rank
    # factorization (ICNTL(35)=0) or BLR factorization with full-rank solve (ICNTL(35)=3). In the
    # case of a BLR factorization with ICNTL(35)=1 or 2, only the L factor corresponding to full-rank
    # frontal matrices are discarded in the current version.
    # We advise to use this option only for a reasonably small number of dense right-hand side vectors
    # because of the additional associated storage required when this option is activated and the number
    # of right-hand sides is large compared to ICNTL(27).
    # === End MUMPS snippet ===

    @param(index=26, page=87)
    class schur_solve_mode:
        "with Schur ON: 0 normal; 1 forward to build reduced RHS; 2 inject reduced solution"

    # === Begin MUMPS snippet: ICNTL(27) page 85 from userguide_5.8.1.txt:4724-4737 ===
    # ICNTL(27) controls the blocking size for multiple right-hand sides.
    #      Phase: accessed by the host during the solution phase
    #      Possible variables/arrays involved: id%NRHS
    #      Possible values :
    #     < 0 : an automatic setting is performed by the solver:
    #           (i) the blocksize = min(id%NRHS, −2×ICNTL(27)) if the factors are on disk
    #           (ICNTL(22)=1);
    #           (ii) the blocksize = min(id%NRHS,−ICNTL(27)) if the factors are in-core (ICNTL(22)=0)
    #       0 : no blocking, it is treated as 1.
    #     > 0 : blocksize = min(id%NRHS,ICNTL(27))
    #      Default value: -32
    #      Remarks: It influences both the memory usage (see INFOG(30) and INFOG(31)) and the
    #      solution time. Larger values of ICNTL(27) lead to larger memory requirements and a better
    #      performance (except if the larger memory requirements induce swapping effects). Tuning
    # === End MUMPS snippet ===

    @param(index=27, page=85)
    class rhs_blocking:
        """
        Blocking size for multiple right-hand sides during solve.

        Default: automatic with value -32. Larger block sizes generally improve
        performance at the cost of memory.

        Value meanings (let b = rhs_blocking, m = n_rhs):
        - b < 0: Automatic.
            - In-core (ICNTL(22)=0): blocksize = min(m, -b)
            - Out-of-core (ICNTL(22)=1): blocksize = min(m, -2*b)
        - b = 0: Treated as 1 (no blocking implies blocksize=1).
        - b > 0: blocksize = min(m, b)
        """

    # === Begin MUMPS snippet: ICNTL(28) page 92 from userguide_5.8.1.txt:5124-5139 ===
    # ICNTL(28)= 0, 1
    # Possible values :
    #   1: Symbolic factorization based on quotient graph, mixing right looking and left looking updates
    #   2: Column count based symbolic factorization based on [29]
    #   Other values are treated as 2.
    #   Default value: 2
    #   Related parameters: ICNTL(7), ICNTL(28)
    #   Remarks: When symbolic factorization is not performed within the ordering (case of ordering
    #   given, ICNTL(7)=1 or centralized Metis ordering, ICNTL(7)=5) then symbolic factorization
    #   will be automatically performed. When SCOTCH is used, ICNTL(7)=3, a fast block
    #   symbolic factorization (exploiting graph separator information) provided within SCOTCH library,
    #   libesmumps.a, is used.
    # === End MUMPS snippet ===

    @param(index=28, page=92)
    class analysis_mode:
        "analysis and ordering mode"

        auto = 0, "automatic (defaults to sequential)"
        sequential = 1, "sequential analysis (ordering via ICNTL(7))"
        parallel = 2, "parallel analysis (ordering via ICNTL(29))"

    # === Begin MUMPS snippet: ICNTL(29) page 85 from userguide_5.8.1.txt:4747-4764 ===
    #      ICNTL(29) parameter is meaningless (choice of the parallel ordering tool).
    #   2: parallel computation. A parallel ordering and parallel symbolic factorization is requested by
    #      the user. For that, one of the parallel ordering tools (or all) must be available, and the matrix
    #      should not be too small. The ordering method is set by ICNTL(29) and the ICNTL(7)
    #      parameter is meaningless.
    # Any other values will be treated as 0.
    # Default value: 0 (automatic choice)
    # Incompatibility: The parallel analysis is not available when the Schur complement feature is
    # requested (ICNTL(19)=1,2 or 3), when a maximum transversal is requested on the input matrix
    # (i.e., ICNTL(6)=1, 2, 3, 4, 5 or 6) or when the input matrix is an unassembled matrices
    # (ICNTL(5)=1). When the number of processes available for parallel analysis is equal to 1,
    # or when the initial matrix is extremely small, a sequential analysis is indeed performed, even if
    # === End MUMPS snippet ===

    @param(index=29, page=85)
    class parallel_ordering:
        "ordering tool in parallel analysis (0 auto -> PT-SCOTCH; 1 PT-SCOTCH; 2 ParMetis)"

    # === Begin MUMPS snippet: ICNTL(30) page 86 from userguide_5.8.1.txt:4790-4820 ===
    # ICNTL(30) computes a user-specified set of entries in the inverse A−1 of the original matrix (see
    # Subsection 5.17.4).
    #      Phase: accessed during the solution phase.
    #      Possible variables/arrays involved: NZ RHS, NRHS, RHS SPARSE, IRHS SPARSE, IRHS PTR
    #      Possible values :
    #        0: no entries in A−1 are computed.
    #        1: computes entries in A−1 .
    #      Other values are treated as 0.
    #      Default value: 0 (no entries in A−1 are computed)
    #      Incompatibility: Error analysis and iterative refinement will not be performed, even if the
    #      corresponding options are set (ICNTL(10) and ICNTL(11)). Because the entries of A−1 are
    #      returned in RHS SPARSE on the host, this functionality is incompatible with the distributed solution
    #      option (ICNTL(21)). Furthermore, computing entries of A−1 is not possible in the case of partial
    #      factorizations with a Schur complement (ICNTL(19)). Option to compute solution using A or
    #      AT (ICNTL(9)) is meaningless and thus ignored.
    #      Related parameters: ICNTL(27)
    #      Remarks: When a set of entries of A−1 is requested, the associated set of columns will be
    #      computed in blocks of size ICNTL(27). Larger ICNTL(27) values will most likely decrease the
    #      amount of factor accesses, enable more parallelism and thus reduce the solution time [48, 44, 13].
    #      The user must specify on input to a call of the solve phase in the arrays IRHS PTR and
    #      IRHS SPARSE the target entries. The array RHS SPARSE should be allocated but not initialized.
    #      Note that since selected entries of the inverse of the matrix are requested, NRHS must be set to N. On
    #      output the arrays IRHS PTR, IRHS SPARSE and RHS SPARSE will hold the requested entries. If
    #      duplicate target entries are provided then duplicate solutions will be returned.
    #      When entries of A−1 are requested (ICNTL(30) = 1), mumps par%RHS needs not be allocated.
    # === End MUMPS snippet ===

    @param(index=30, page=86)
    class inverse_entries:
        "compute selected entries of A^-1 during solve when nonzero (uses sparse RHS interface; overrides ICNTL(9))"

    # === Begin MUMPS snippet: ICNTL(31) page 87 from userguide_5.8.1.txt:4822-4855 ===
    # ICNTL(31) indicates which factors may be discarded during the factorization.
    #      Phase: accessed by the host during the analysis phase.
    #      Possible values :
    #        0 : the factors are kept during the factorization, phase except in the case of unsymmetric matrices
    #            when the forward elimination is performed during factorization (ICNTL(32) = 1). In this
    #            case, since it will not be used during the solve phase, the L factor is discarded.
    #         1: all factors are discarded during the factorization phase. The user is not interested in solving
    #            the linear system (Equations (3) or (4)) and will not call MUMPS solution phase (JOB= 3).
    #            This option is meaningful when only a Schur complement is needed (see ICNTL(19)), or
    #            when only statistics from the factorization, such as (for example) definiteness, value of the
    #            determinant, number of entries in factors after numerical pivoting, number of negative or null
    #            pivots are required. In this case, the memory allocated for the factorization will rely on the
    #            out-of-core estimates (and factors will not be written to disk).
    #         2: this setting is meaningful only for unsymmetric matrices and has no impact on symmetric
    #            matrices: only the U factor is kept after factorization so that exclusively a backward
    #            substitution is possible during the solve phase (JOB= 3). This can be useful when:
    #            −the user is only interested in the computation of a null space basis (see ICNTL(25))
    #            during the solve phase, or
    #            −the forward elimination is performed during the factorization (ICNTL(32)=1). Note that
    #            for unsymmetric matrices, if the forward elimination is performed during the factorization
    #            (ICNTL(32) = 1) then the L factor is always discarded during factorization. In this case
    #            (ICNTL(32) = 1), both ICNTL(31) = 0 and ICNTL(31) = 2 have the same behaviour.
    #
    #      Other values are treated as 0.
    #      Default value: 0 (the factors are kept during the factorization phase in order to be able to solve the
    #      system).
    #      Incompatibility: ICNTL(31) = 2 is not meaningful for symmetric matrices.
    #      Related parameters: ICNTL(32), forward elimination during factorization, ICNTL(33),
    #      computation of the determinant, ICNTL(25) computation of a null space basis, ICNTL(22)
    #      out-of-core factors.
    #      Remarks: For unsymmetric matrices and ICNTL(32)=2, MUMPS currently discards the L
    #      factors corresponding to full-rank frontal matrices but not of low-rank frontal matrices (except
    #      if ICNTL(35)=3). In a future version, discarding all the L factors in case of BLR factorization
    #      and ICNTL(32)=2 may lead to a further memory reduction.
    # === End MUMPS snippet ===

    @param(index=31, page=87)
    class param_31:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(32) page 87 from userguide_5.8.1.txt:4857-4860 ===
    # ICNTL(32) performs the forward elimination of the right-hand sides (Equation (3)) during the
    # factorization (JOB= 2). (see Subsection 5.16).
    #      Phase: accessed by the host during the analysis phase.
    #      Possible variables/arrays involved: RHS, NRHS, LRHS, and possibly REDRHS, LREDRHS when
    # === End MUMPS snippet ===

    @param(index=32, page=87)
    class param_32:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(33) page 88 from userguide_5.8.1.txt:4894-4926 ===
    # ICNTL(33) computes the determinant of the input matrix.
    #      Phase: accessed by the host during the factorization phase.
    #      Possible values :
    #       0 : the determinant of the input matrix is not computed.
    #     ̸= 0: computes the determinant of the input matrix. The determinant is obtained by computing
    #           (a + ib) × 2c where a =RINFOG(12), b =RINFOG(13) and c = INFOG(34). In real
    #           arithmetic b=RINFOG(13) is equal to 0.
    #      Default value: 0 (determinant is not computed)
    #      Related parameters: ICNTL(31)
    #      Remarks: In case a Schur complement was requested (see ICNTL(19)), elements of the Schur
    #      complement are excluded from the computation of the determinant so that the determinant is that
    #      one of matrix A1,1 (using notations of Subsection 3.18).
    #      Although we recommend to compute the determinant on non-singular matrices, null pivot rows
    #      (ICNTL(24)) and static pivots (CNTL(4)) are excluded from the determinant so that a non-zero
    #      determinant is still returned on singular or near-singular matrices. This determinant is then not
    #      unique and will depend on which equations were excluded.
    #      Furthermore, we recommend to switch off scaling (ICNTL(8)) in such cases. If not (ICNTL(8)
    #      ̸= 0), we describe in the following the current behaviour of the package:
    #         – if static pivoting (CNTL(4)) is activated: all entries of the scaling arrays ROWSCA and
    #           COLSCA are currently taken into account in the computation of the determinant.
    #         – if the null pivot row detection (ICNTL(24)) is activated, then entries of ROWSCA and
    #           COLSCA corresponding to pivots in PIVNUL LIST are excluded from the determinant so
    #           that
    #              * for symmetric matrices (SYM=1 or 2), the returned determinant correctly corresponds to
    #                the matrix excluding rows and columns of PIVNUL LIST.
    #              * for unsymmetric matrices (SYM=0), scaling may perturb the value of the determinant in
    #                case off-diagonal pivoting has occurred (INFOG(12)̸=0).
    #      Note that if the user is interested in computing only the determinant, we recommend to discard the
    #      factors during factorization ICNTL(31).
    # === End MUMPS snippet ===

    @param(index=33, page=88)
    class determinant:
        "determinant computation before factorization (0 off; otherwise on; excludes null/static pivots)"

    # === Begin MUMPS snippet: ICNTL(34) page 89 from userguide_5.8.1.txt:4928-4939 ===
    # ICNTL(34) controls the conservation of the OOC files during JOB= –3 (See Subsection 5.20).
    #      Phase: accessed by the host during the save/restore files deletion phase (JOB= –3) in case of out-
    #      of-core (ICNTL(22)=1).
    #      Possible values :
    #        0: the out-of-core files are marked out for deletion
    #        1: the out-of-core files should not be deleted because another saved instance references them.
    #      Other values are treated as 0.
    #      Default value: 0 (out-of-core files associated to a saved instance are marked out for deletion at the
    #      end of the out-of-core file lifetime)
    #      Remarks: MUMPS will delete only the out-of-core files that are referenced in the saved data
    #      identified by the value of SAVE DIR and SAVE PREFIX. Extra out-of-core files with the same
    #      OOC TMPDIR and OOC PREFIX are not deleted.
    # === End MUMPS snippet ===

    @param(index=34, page=89)
    class save_restore_cleanup:
        "when JOB=-3 and save/restore used: 0 keep restore files; 1 clean files"

    # === Begin MUMPS snippet: ICNTL(35) page 89 from userguide_5.8.1.txt:4941-4972 ===
    # ICNTL(35) controls the activation of the BLR feature (see Subsection 5.19).
    #      Phase: accessed by the host during the analysis and during the factorization phases
    #      Possible values :
    #        0 : Standard analysis and factorization (BLR feature is not activated).
    #        1 : BLR feature is activated and automatic choice of BLR option is performed by the software.
    #        2 : BLR feature is activated during both the factorization and solution phases, which allows for
    #            memory gains by storing the factors in low-rank.
    #        3 : BLR feature is activated during the factorization phase but not the solution phase, which is still
    #            performed in full-rank. As a consequence, the full-rank factors must be kept and no memory
    #            gains can be obtained. In an OOC context, (ICNTL(22)=1) this option enables the user to
    #            write all factors to disk which is not the case with ICNTL(35)=2 since factors in low-rank
    #            form are not written to disk.
    #      Other values are treated as 0.
    #      Default value: 0 (standard multifrontal factorization).
    #      Related parameters: CNTL(7) (BLR approximations accuracy), ICNTL(36) (BLR factorization
    #      variant), ICNTL(37) (compression of the contribution blocks), ICNTL(38) (estimation of the
    #      compression rate of the factors) and ICNTL(39) (estimation of the compression rate of the
    #      contribution blocks).
    #      Incompatibility: Note that the activation of the BLR feature is currently incompatible with
    #      elemental matrices (ICNTL(5) = 1) (see error -800, subject to change in the future), and when
    #      the forward elimination during the factorization is requested (ICNTL(32) = 1), see error -43.
    #      Remarks: If ICNTL(35)=1, then the automatic choice of BLR option is to activate BLR
    #      feature during both factorization and solution phases (ICNTL(35)=2). In order to activate the
    #      BLR factorization, ICNTL(35) must be equal to 1, 2 or 3 before the analysis, where some
    #      preprocessing on the graph of the matrix is needed to prepare the low-rank factorization. The value
    #      of ICNTL(35) can then be set to any of the above values on entry to the factorization (e.g., taking
    #      into account the values returned by the analysis). On the other hand, if ICNTL(35)=0 at analysis,
    #      only ICNTL(35)=0 is allowed for the factorization (full-rank factorization). When activating
    #      BLR, it is recommended to set ICNTL(35) to 1 or 2 rather than 3 to benefit from memory gains.
    # === End MUMPS snippet ===

    @param(index=35, page=89)
    class param_35:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(36) page 90 from userguide_5.8.1.txt:4974-4986 ===
    # ICNTL(36) controls the choice of BLR factorization variant (see Subsection 5.19).
    #      Phase: accessed by the host during the factorization phase when ICNTL(35)=1, 2 or 3
    #      Possible values :
    #        0 : Standard UFSC variant with low-rank updates accumulation (LUA)
    #        1 : UCFS variant with low-rank updates accumulation (LUA). This variant consists in performing
    #            the compression earlier in order to further reduce the number of operations. Although it
    #            may have a numerical impact, the current implementation is still compatible with numerical
    #            pivoting.
    #      Other values are treated as 0.
    #      Default value: 0 (UFSC variant).
    #      Related parameters: ICNTL(35) and CNTL(1)
    #      Remarks: If numerical pivoting is not required and thus CNTL(1) can be set to 0.0, further
    #      performance gains can be expected with the UCFS version.
    # === End MUMPS snippet ===

    @param(index=36, page=90)
    class param_36:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(37) page 90 from userguide_5.8.1.txt:4988-4998 ===
    # ICNTL(37) controls the BLR compression of the contribution blocks (see Subsection 5.19).
    #      Phase: accessed by the host during the factorization phase when ICNTL(35)=1, 2 or 3
    #      Possible values :
    #        0 : contribution blocks are not compressed
    #        1 : contribution blocks are compressed, reducing the memory consumption at the cost of some
    #            additional operations
    #      Other values are treated as 0.
    #      Default value: 0 (contribution blocks not compressed).
    #      Related parameters: ICNTL(35), CNTL(7), ICNTL(40)
    #      Remarks: This feature should be activated if memory consumption is a primary concern; note
    #      that it is likely to increase the factorization time. We recommend to combine this feature with
    # === End MUMPS snippet ===

    @param(index=37, page=90)
    class param_37:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(38) page 90 from userguide_5.8.1.txt:5002-5013 ===
    # ICNTL(38) estimates compression rate of LU factors (see Subsection 5.19).
    #      Phase: accessed by the host during the analysis and the factorization phases when ICNTL(35)=1,
    #      2 or 3
    #      Possible values : between 0 and 1000 (1000 is no compression and 0 is full compression); other
    #      values are treated as 0; ICNTL(38)/10 is a percentage representing the typical compression of the
    #                                                         compressed factors
    #      factor matrices in BLR fronts: ICNTL(38)/10 = uncompressed factors × 100.
    #      Default value: 600 (when factors of BLR fronts are compressed, their size is 60.0% of their full-
    #      rank size).
    #      Related parameters: ICNTL(35), CNTL(7)
    #      Remarks: Influences statistics provided in INFO(29), INFO(30), INFO(31), INFOG(36),
    #      INFOG(37), INFOG(38), INFOG(39) , but also INFO(32-35) and INFOG(40-43)
    # === End MUMPS snippet ===

    @param(index=38, page=90)
    class param_38:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(39) page 83 from userguide_5.8.1.txt:4620-4642 ===
    #      ICNTL(39) are good approximations of the average compression rate of respectively the factors and the
    #      CB) is given by INFOG(40) if the factorization is in-core (ICNTL(22)=0), and by INFOG(42) if the
    #      factorization is out-of-core (ICNTL(22)=1).
    # Lower bounds for ICNTL(23), in case ICNTL(23) is provided locally to each MPI process:
    #    – Full-rank factors only (ICNTL(35)=0 or 3) ⇒ INFO(15) if the factorization is in-core
    #      (ICNTL(22)=0), INFO(17) if the factorization is out-of-core (ICNTL(22)=1).
    #    – Low-rank factors (ICNTL(35)=1 or 2) only (ICNTL(37)=0) ⇒ INFO(30) if the factorization is
    #      in-core (ICNTL(22)=0), INFO(31) if the factorization is out-of-core (ICNTL(22)=1).
    #    – Low-rank factors and contribution blocks (ICNTL(35)=1,2 and ICNTL(37)=1) ⇒ is given by
    #      INFO(34) if the factorization is in-core (ICNTL(22)=0), INFO(35) if the factorization is out-of-
    #      core (ICNTL(22)=1).
    # The above lower bounds include memory for the real/complex internal workarray S holding the
    # factors and stack of contribution blocks. In case WK USER is provided, the above quantities should
    # be diminished by the estimated memory for S/WK USER. This estimated memory can be obtained
    # from INFO(8), INFO(9), or INFO(20) (depending on MUMPS settings) by taking their absolute
    # value, if negative, or by dividing them by 106 , if positive. See also the paragraph Recommended
    # values of LWK USER below.
    # If ICNTL(23) is left to its default value 0 then MUMPS will allocate for the factorization phase
    # a workspace based on the estimates computed during the analysis if ICNTL(14) has not been
    # modified since analysis, or larger if ICNTL(14) was increased. Note that even with full-rank
    # factorization, these estimates are only accurate in the sequential version of MUMPS but they can
    # be inaccurate in the parallel case, especially for the out-of-core version. Therefore, in parallel, we
    # recommend to use ICNTL(23) and provide a value larger than the provided estimations.
    # === End MUMPS snippet ===

    @param(index=39, page=83)
    class param_39:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: ICNTL(48) page 91 from userguide_5.8.1.txt:5037-5058 ===
    # ICNTL(48) multithreading with tree parallelism (see Subsection 5.23).
    #      Phase: accessed by the host during the analysis phase (need be set after the initialization phase
    #      (JOB= –1)) and during the solution phase
    #      Possible values :
    #        0 : not activated
    #        1 : multithreaded tree parallelism activated
    #      Other values are treated as 0.
    #      Default value: 1
    #      Related parameters: if tree parallelism is activated, then the number of threads per MPI process set
    #      through the OMP NUM THREADS environment variable or ICNTL(16) control parameter should
    #      be at least 2. It will influence the choice of the L0 layer at the analysis phase and should not be
    #      modified during the subsequent numerical phases.
    #      Remarks:
    #         – Please note that, once ICNTL(48) = 1 is set prior to the analysis phase it will be used for
    #           both analysis and factorization phases. Tree parallelism can be switched off (ICNTL(48)=0)
    #           prior to the solution phase (JOB= 3).
    #         – The number of threads per MPI process should be identical and if not multithreaded tree
    #           parallelism will be switched off during analysis.
    #         – If tree parallelism is not effective during analysis then it cannot be activated during the
    #           factorization and solve phases.
    #         – The memory provided with WK USER will not be used under the L0 layer making WK USER
    #           not recommended when ICNTL(48) = 1.
    # === End MUMPS snippet ===

    @param(index=48, page=91)
    class l0_omp:
        "Controls L0_OMP feature (analysis/factorization/solve activation); out-of-range treated as 0"

    # === Begin MUMPS snippet: ICNTL(49) page 91 from userguide_5.8.1.txt:5060-5073 ===
    # ICNTL(49) compact workarray id%S at the end of factorization phase (see Subsection 5.22).
    #      Phase: accessed by the host during factorization phase
    #      Possible values :
    #        0 : nothing is done.
    #        1 : compact workarray id%S(MAXS) at the end of the factorization phase while satisfying the
    #            memory constraint that might have been provided with ICNTL(23) feature.
    #        2 : compact workarray id%S(MAXS) at the end of the factorization phase. The memory
    #            constraint that might have been provided with ICNTL(23) feature does not apply to this
    #            process.
    #      Other values are treated as 0.
    #      Default value: 0
    #      Incompatibility: With the use of LWK USER / WK USER feature.
    #      Remarks: ICNTL(49)=1,2 might require intermediate memory allocation to reallocate id%S of
    #      minimal size. If the memory allocation fails, then a warning is returned and nothing is done. If
    # === End MUMPS snippet ===

    @param(index=49, page=91)
    class compact_workarray:
        "compact workarray before solve: 0 none; 1 compact respecting ICNTL(23); 2 compact ignoring it"

    # === Begin MUMPS snippet: ICNTL(56) page 92 from userguide_5.8.1.txt:5087-5116 ===
    # ICNTL(56) detects pseudo-singularities during factorization and factorizes the root node with a rank-
    # revealing method. See Subsection 5.13 for details
    #      Phase: accessed by the host during the analysis and factorization phases
    #      Possible values :
    #        0 : standard factorization is performed.
    #        1 : Postponing and rank-revealing factorization on root node based on a singular value
    #            decomposition.
    #      Values different from 1 are treated as 0.
    #      Default value: 0 (standard factorization)
    #      Related parameters: ICNTL(13), ICNTL(24), ICNTL(25), CNTL(1), CNTL(3)
    #      Remarks: Note that ICNTL(56) must be set before analysis, during which a valid positive value
    #      prepares the data for later use of the rank-revealing functionality.
    #      The root is processed sequentially and ICNTL(13) setting is ignored.
    #      Please note that to improve the numerical behaviour of the factorization, the default value of
    #      CNTL(1) has been increased in case of activation of the rank-revealing feature. On numerically
    #      difficult problems the value of CNTL(1) may be further increased.
    #      CNTL(3) it is used to determine if a pivot should be postponed.
    #      Related data:
    #        mumps par%INFOG(28) (integer):
    #            INFOG(28) is set during factorization to the deficiency of the matrix. It counts null pivot
    #            rows if ICNTL(24) = 1 and null singular values (ICNTL(56) = 1).
    #        mumps par%PIVNUL LIST (integer array, dimension N):
    #            If INFOG(28) ̸= 0 then PIVNUL LIST(1:INFOG(28)) will hold, on the host, the row
    #            indices corresponding to both the null pivot rows (if ICNTL(24)= 1) and the null singular
    #            values (ICNTL(56) = 1).
    #        mumps par%SINGULAR VALUES(1:mumps par%NB SINGULAR VALUES) (real pointer
    #            array) which holds, on the host, all the singular values corresponding to SVD decomposition
    #            on the root node.
    #        If the matrix was found to be deficient (INFOG(28) > 0), the solution phase (JOB= 3) can
    #            then be used to either provide a “regular” solution or to compute the null-space basis (see
    # === End MUMPS snippet ===

    @param(index=56, page=92)
    class null_space_analysis:
        "null space prep/analysis (analysis: >0 prepare; factorization: 1 SVD, 2 QR if prepared)"

    # === Begin MUMPS snippet: ICNTL(58) page 92 from userguide_5.8.1.txt:5122-5123 ===
    # ICNTL(58) defines options for symbolic factorization
    #      Phase: accessed by the host during the analysis when centralized ordering is performed,
    # === End MUMPS snippet ===

    @param(index=58, page=92)
    class symb_factorization:
        "symbolic factorization strategy with METIS/given ordering (1 SYMBQAMD; 2 column count; default 2)"


# ------------
# CNTL members
# ------------


class CNTL(ParamArray):
    """Floating-point control parameters CNTL(1..15)."""

    def __init__(self, size: int = LEN_CNTL):
        # floats by default
        self._array = _OneBasedArray(size, default=0.0)
        self._array_alias = type(self)._array_alias or type(self).__name__.lower()

    # === Begin MUMPS snippet: CNTL(1) page 93 from userguide_5.8.1.txt:5145-5162 ===
    # CNTL(1) is the relative threshold for numerical pivoting. See Subsection 3.10
    #        Phase: accessed by the host during the factorization phase.
    #        Possible values :
    #    < 0.0: Automatic choice
    #    = 0.0: no numerical pivoting performed and the subroutine will fail if a zero pivot is encountered.
    #    > 0.0: numerical pivoting performed.
    #           For unsymmetric matrices values greater than 1.0 are treated as 1.0
    #           For symmetric matrices values greater than 0.5 are treated as 0.5
    #        Default value: -1.0 (automatic choice):
    #        0.1: in case of rank-revealing (ICNTL(56)= 1)
    #       0.01: for unsymmetric or general symmetric matrices
    #        0.0: for symmetric positive definite matrices
    #        Related parameters: CNTL(4)
    #        Remarks: It forms a trade-off between preserving sparsity and ensuring numerical stability during
    #        the factorization. In general, a larger value of CNTL(1) increases fill-in but leads to a more accurate
    #        factorization.
    #        Note that for diagonally dominant matrix, setting CNTL(1) to zero will decrease the factorization
    #        time while still providing a stable decomposition.
    # === End MUMPS snippet ===

    @param(index=1, page=93)
    class threshold:
        "relative pivoting threshold"

    # === Begin MUMPS snippet: CNTL(2) page 93 from userguide_5.8.1.txt:5164-5176 ===
    # CNTL(2) is the stopping criterion for iterative refinement
    #        Phase: accessed by the host during the solve phase.
    #        Possible values :
    #                                           √
    #    < 0.0: values < 0 are treated as        ϵ, where ϵ holds the machine precision and depends on the
    #           arithmetic version.
    #    ≥ 0.0: stopping criterion
    #                     √
    #      Default value: ϵ
    #        Related parameters: ICNTL(10), RINFOG(7), RINFOG(8)
    #        Remarks: Let ω1 and ω2 be the backward errors as defined in Subsection 3.3.2. Iterative refinement
    #        (Subsection 5.8) will stop when either the requested accuracy is reached (ω1 + ω2 < CNTL(2))
    #        or when the convergence rate is too slow (ω1 + ω2 does not decrease by at least a factor of 2).
    # === End MUMPS snippet ===

    @param(index=2, page=93)
    class stopping:
        "stopping criterion (iterative refinement tolerance)"

    # === Begin MUMPS snippet: CNTL(3) page 93 from userguide_5.8.1.txt:5178-5205 ===
    # CNTL(3) it is used to determine null pivot rows when the null pivot row detection option is enabled
    # (ICNTL(24) = 1) and/or singularities at the root node when the rank-revealing option is enabled
    # (ICNTL(56) = 1). It is also used to determine small pivots whose elimination will be postponed if
    # rank-revealing option is enabled (ICNTL(56) = 1).
    #       Phase: accessed by the host during the numerical factorization phase.
    #       Possible values : we define the threshold thres as follows
    #    > 0.0: thres = CN T L(3) × ∥Apre ∥
    #                                √
    #    = 0.0: thres = ϵ × ∥Apre ∥ × Nh
    #    < 0.0: thres = |CN T L(3)|
    #       where Apre is the preprocessed matrix to be factorized (see Equation (5)), Nh number of variables
    #       on the deepest branch of the elimination tree, ϵ is the machine precision and ∥.∥ is the infinite norm.
    #       Default value: 0.0
    #       Related parameters: ICNTL(24), ICNTL(56)
    #       Remarks:
    #         – If rank-revealing is enabled (ICNTL(56)=1) then
    #              * singular values smaller than thres are considered as null.
    #              * thres is also used to automatically define pivot rows that are small (with respect to the
    #                infinite norm of its row/column) and that should be postponed possibly up to the root
    #                node on which rank-revealing will be performed.
    #              * when null pivot row detection is also enabled (ICNTL(24)=1), the threshold value to
    #                determine null pivot rows also relies on thres.
    #         – When only null pivot row detection is enabled (ICNTL(24)=1). A pivot is considered to be
    #           null if the infinite norm of its row/column is smaller than thres.
    # === End MUMPS snippet ===

    @param(index=3, page=93)
    class rr_thresholds:
        "thresholds for RR and null pivot detection (see docs)"

    # === Begin MUMPS snippet: CNTL(4) page 94 from userguide_5.8.1.txt:5207-5225 ===
    # CNTL(4) determines the threshold for static pivoting. See Subsection 3.10
    #       Related parameters: CNTL(1), INFOG(25)
    #       Phase: accessed by the host, and must be set either before the factorization phase, or before the
    #       analysis phase.
    #       Possible values :
    #    < 0.0: static pivoting is not activated.
    #    > 0.0: static pivoting is activated and the pivots whose magnitude is smaller than CNTL(4) will be
    #           set to CNTL(4).
    #    = 0.0: static pivoting is activated and the threshold value to define a small pivot is determined
    #                                                                            √
    #           automatically. In the current version, this threshold is equal to ϵ × ∥Apre ∥, where Apre is
    #           the preprocessed matrix to be factored (see Equation (5)).
    #       Default value: -1.0 (no static pivoting)
    #       Related parameters: CNTL(1)
    #       Incompatibility: This option is incompatible with null pivot row detection (ICNTL(24)= 1) or
    #       with rank-revealing factorization (ICNTL(56) = 1) and will be ignored.
    #       Remarks: By static pivoting (as in [39]) we mean replacing small pivots whose elimination should
    #       be postponed because of partial threshold pivoting and would thus result in an increase of our
    #       estimations (memory and operations), by a small perturbation of the original matrix controlled by
    # === End MUMPS snippet ===

    @param(index=4, page=94)
    class static_pivoting:
        "static pivoting threshold (-1 off; 0 MACHEPS^1/2 ||A||; >0 threshold)"

    # === Begin MUMPS snippet: CNTL(5) page 94 from userguide_5.8.1.txt:5228-5246 ===
    # CNTL(5) defines the fixation for null pivots and is effective only when null pivot row detection is active
    # (ICNTL(24) = 1).
    #       Phase: accessed by the host during the numerical factorization phase.
    #       Possible values :
    #    ≤ 0.0: In the symmetric case (SYM = 2), the pivot column of the L factors is set to zero and the pivot
    #           entry in matrix D is set to one.
    #           In the unsymmetric case (SYM = 0), the fixation is automatically set to a large positive value
    #           and the pivot row of the U factors is set to zero.
    #    > 0.0: when a pivot piv is detected as null, in order to limit the impact of this pivot on the rest of
    #           the matrix, it is set to sign(piv) CNTL(5) ×∥Apre ∥, where Apre is the preprocessed matrix
    #           to be factored (see Equation (5)). We recommend setting CNTL(5) to a large floating-point
    #           value (e.g. 1020 ).
    #       Default value: 0.0
    #       Related parameters: ICNTL(24)
    # === End MUMPS snippet ===

    @param(index=5, page=94)
    class null_pivot_fixation:
        "null pivot fixation control (active if ICNTL(24)=1)"

    # === Begin MUMPS snippet: CNTL(7) page 95 from userguide_5.8.1.txt:5249-5267 ===
    # CNTL(7) defines the precision of the dropping parameter used during BLR compression (see
    # Subsection 5.19).
    #       Phase: accessed by the host during the factorization phase when ICNTL(35)=1, 2 or 3
    #       Possible values :
    #     0.0 : full precision approximation.
    #   > 0.0 : the dropping parameter is CNTL(7).
    #       Default value: 0.0 (full precision (i.e., no approximation)).
    #       Related parameters: ICNTL(35)
    #       Remarks: The value of CNTL(7) is used as a stopping criterion for the compression of BLR
    #       blocks which is achieved through a truncated Rank Revealing QR factorization. More precisely, to
    #       compute the low-rank form of a block, we perform a QR factorization with column pivoting which
    #       is stopped as soon as a diagonal coefficient of the R factor falls below the threshold, i.e., when
    #       ∥rkk ∥ < ε. This is implemented as a variant of the LAPACK [18] GEQP3 routine. Larger values
    #       of this parameter lead to more compression at the price of a lower accuracy. Note that ϵ is used as
    #       an absolute tolerance, i.e., not relative to the input matrix, or the frontal matrix or the block norms;
    #       for this reason we recommend to scale the matrix or or let the solver automatically preprocess (e.g.,
    #       scale) the input matrix.
    #       Note that, depending on the application, gains can be expected even with small values (close to
    #       machine precision) of CNTL(7).
    # === End MUMPS snippet ===

    @param(index=7, page=95)
    class blr_tolerance:
        "dropping parameter for BLR truncated RRQR (0.0 => no approximation)"


# ------------
# INFO members
# ------------


class INFO(RawArray):
    """Integer info array INFO(1..80). Read-only in most workflows."""

    def __init__(self, size: int = LEN_INFO):
        super().__init__(size, default=0)

    # === Begin MUMPS snippet: INFO(1) page 98 from userguide_5.8.1.txt:5412-5415 ===
    # INFO(1) is 0 if the call to MUMPS was successful, negative if an error occurred (see Section 8), or
    #    positive if a warning is returned. In particular, after successfully saving or restoring an instance
    #    (call to MUMPS with JOB= 7 or JOB= 8), INFO(1) will be 0 even if INFO(1) was different from
    #    0 at the moment of saving the MUMPS instance to disk.
    # === End MUMPS snippet ===

    @param(index=1, page=98)
    class param_1:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(2) page 98 from userguide_5.8.1.txt:5416-5417 ===
    # INFO(2) holds additional information about the error or the warning. If INFO(1) = -1, INFO(2) is
    #    the processor number (in communicator COMM) on which the error was detected.
    # === End MUMPS snippet ===

    @param(index=2, page=98)
    class param_2:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(3) page 98 from userguide_5.8.1.txt:5418-5423 ===
    # INFO(3) - after analysis: Estimated size of the real/complex space needed on the processor to
    #    store the factors, assuming the factors are stored in full-rank format (ICNTL(35)=0 or 3 during
    #    factorization). If INFO(3) is negative, then its absolute value corresponds to millions of
    #    real/complex entries used to store the factor matrices. Assuming that the factors will be stored
    #    in full-rank format during the factorization (ICNTL(35)=0 or 3), a rough estimation of the size of
    #    the disk space in bytes of the files written by the concerned processor can be obtained by multiplying
    # === End MUMPS snippet ===

    @param(index=3, page=98)
    class param_3:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(4) page 98 from userguide_5.8.1.txt:5436-5437 ===
    # INFO(4) - after analysis: Estimated integer space needed on the processor for factors (assuming a
    #    full-rank storage for the factors)
    # === End MUMPS snippet ===

    @param(index=4, page=98)
    class param_4:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(5) page 98 from userguide_5.8.1.txt:5438-5438 ===
    # INFO(5) - after analysis: Estimated maximum front size on the processor.
    # === End MUMPS snippet ===

    @param(index=5, page=98)
    class param_5:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(6) page 98 from userguide_5.8.1.txt:5439-5440 ===
    # INFO(6) - after analysis: Number of nodes in the complete tree. The same value is returned on all
    #    processors.
    # === End MUMPS snippet ===

    @param(index=6, page=98)
    class param_6:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(7) page 98 from userguide_5.8.1.txt:5441-5442 ===
    # INFO(7) - after analysis: Minimum estimated size of the main internal integer workarray IS to run the
    #    numerical factorization in-core .
    # === End MUMPS snippet ===

    @param(index=7, page=98)
    class param_7:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(8) page 98 from userguide_5.8.1.txt:5443-5452 ===
    # INFO(8) - after analysis: Minimum estimated size of the main internal real/complex workarray S to
    #    run the numerical factorization in-core when factors are stored full-rank (ICNTL(35)=0 or 3).
    #    If negative, then the absolute value corresponds to millions of real/complex entries needed in this
    #    workarray. It is also the estimated minimum size of LWK USER in that case, if the user provides
    #    WK USER.
    # === End MUMPS snippet ===

    @param(index=8, page=98)
    class param_8:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(9) page 99 from userguide_5.8.1.txt:5454-5458 ===
    # INFO(9) - after factorization: Size of the real/complex space used on the processor to store the factor
    #    matrices, possibly including low-rank factor matrices (ICNTL(35)=1 or 2). If negative, then the
    #    absolute value corresponds to millions of real/complex entries used to store the factor matrices.
    #    Finally, the total size of the factor matrices for all processors (sum of the INFO(9) values over all
    #    processors) is returned in INFOG(9).
    # === End MUMPS snippet ===

    @param(index=9, page=99)
    class param_9:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(10) page 99 from userguide_5.8.1.txt:5459-5460 ===
    # INFO(10) - after factorization: Size of the integer space used on the processor to store the factor
    #    matrices.
    # === End MUMPS snippet ===

    @param(index=10, page=99)
    class param_10:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(11) page 99 from userguide_5.8.1.txt:5461-5461 ===
    # INFO(11) - after factorization: Order of the largest frontal matrix processed on the processor.
    # === End MUMPS snippet ===

    @param(index=11, page=99)
    class param_11:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(12) page 99 from userguide_5.8.1.txt:5462-5473 ===
    # INFO(12) - after factorization: Number of off-diagonal pivots selected on the processor if SYM=0
    #    or number of negative pivots on the processor if SYM=1 or 2. If ICNTL(13)=0 (the default),
    #    this excludes pivots from the parallel root node treated by ScaLAPACK. (This means that the user
    #    should set ICNTL(13)=1 or use a single processor in order to get the exact number of off-diagonal
    #    or negative pivots rather than a lower bound.) If rank-revealing option is on (ICNTL(56)=1), the
    #    inertia might be incorrect if INFO(1) ≥ 16. Furthermore, when ICNTL(24) is set to 1 and
    #    SYM=1 or 2, or when ICNTL(56)=1, INFOG(12) excludes the null12 pivots, even if their sign
    #    is negative. In other words, a pivot cannot be both null and negative. Note that for complex
    #    symmetric matrices (SYM=1 or 2), INFO(12) will be 0. See also INFOG(12), which provides
    #    the total number of off-diagonal or negative pivots over all processors. For real symmetric matrices,
    #    see also INFO(40) and INFOG(50), which provide the local (resp. global) number of negative
    #    pivots among the null pivots detected when ICNTL(24) (and/or ICNTL(56))) is activated.
    # === End MUMPS snippet ===

    @param(index=12, page=99)
    class param_12:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(13) page 99 from userguide_5.8.1.txt:5474-5474 ===
    # INFO(13) - after factorization: The number of postponed elimination because of numerical issues.
    # === End MUMPS snippet ===

    @param(index=13, page=99)
    class param_13:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(14) page 99 from userguide_5.8.1.txt:5475-5475 ===
    # INFO(14) - after factorization: Number of memory compresses.
    # === End MUMPS snippet ===

    @param(index=14, page=99)
    class param_14:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(15) page 99 from userguide_5.8.1.txt:5476-5480 ===
    # INFO(15) - after analysis: estimated size in MegaBytes (millions of bytes) of all working space
    #    to perform full-rank numerical phases (factorization/solve) in-core (ICNTL(22)=0 for the
    #    factorization). The maximum and sum over all processors are returned respectively in INFOG(16)
    #    and INFOG(17). See also INFO(22) which provides the actual memory that was needed but
    #    only after factorization.
    # === End MUMPS snippet ===

    @param(index=15, page=99)
    class param_15:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(16) page 99 from userguide_5.8.1.txt:5481-5484 ===
    # INFO(16) - after factorization: total size (in millions of bytes) of all MUMPS internal data allocated
    #    during the numerical factorization. This excludes the memory for WK USER, in the case where
    #    WK USER is provided. The maximum and sum over all processors are returned respectively in
    #    INFOG(18) and INFOG(19).
    # === End MUMPS snippet ===

    @param(index=16, page=99)
    class param_16:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(17) page 99 from userguide_5.8.1.txt:5485-5487 ===
    # INFO(17) - after analysis: estimated size in MegaBytes (millions of bytes) of all working space to
    #    run the numerical phases out-of-core (ICNTL(22)̸=0) with the default strategy. The maximum
    #    and sum over all processors are returned respectively in INFOG(26) and INFOG(27). See also
    # === End MUMPS snippet ===

    @param(index=17, page=99)
    class param_17:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(18) page 99 from userguide_5.8.1.txt:5489-5490 ===
    # INFO(18) - after factorization: local number of null pivot rows detected locally when ICNTL(24)=1
    #    or ICNTL(56)= 1.
    # === End MUMPS snippet ===

    @param(index=18, page=99)
    class param_18:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(19) page 99 from userguide_5.8.1.txt:5491-5492 ===
    # INFO(19) - after analysis: Estimated size of the main internal integer workarray IS to run the
    #    numerical factorization out-of-core .
    # === End MUMPS snippet ===

    @param(index=19, page=99)
    class param_19:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(20) page 99 from userguide_5.8.1.txt:5493-5496 ===
    # INFO(20) - after analysis: Estimated size of the main internal real/complex workarray S to run the
    #    numerical factorization out-of-core . If negative, then the absolute value corresponds to millions of
    #    real/complex entries needed in this workarray. It is also the estimated minimum size of LWK USER
    #    in that case, if the user provides WK USER.
    # === End MUMPS snippet ===

    @param(index=20, page=99)
    class param_20:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(21) page 99 from userguide_5.8.1.txt:5497-5499 ===
    # INFO(21) - after factorization: Effective space used in the main real/complex workarray S– or in the
    #    workarray WK USER, in the case where WK USER is provided. If negative, then the absolute value
    #    corresponds to millions of real/complex entries needed in this workarray.
    # === End MUMPS snippet ===

    @param(index=21, page=99)
    class param_21:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(22) page 99 from userguide_5.8.1.txt:5500-5512 ===
    #       INFO(22) - after factorization: Size in millions of bytes of memory effectively used during
    #          factorization. This includes the part of the memory effectively used from the workarray WK USER,
    #          in the case where WK USER is provided. The maximum and sum over all processors are
    #          returned respectively in INFOG(21) and INFOG(22). The difference between estimated and
    # 12 i.e., whose magnitude is smaller than the tolerance defined by CNTL(3).
    #    effective memory may results from numerical pivoting difficulties, parallelism and BLR effective
    #    compression rates.
    # === End MUMPS snippet ===

    @param(index=22, page=99)
    class param_22:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(23) page 100 from userguide_5.8.1.txt:5513-5515 ===
    # INFO(23) - after factorization: total number of pivots eliminated on the processor. It may be used
    #    in the case of a distributed right-hand side (see ICNTL(20)) and/or of a distributed solution (see
    #    ICNTL(21)).
    # === End MUMPS snippet ===

    @param(index=23, page=100)
    class param_23:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(24) page 100 from userguide_5.8.1.txt:5516-5520 ===
    # INFO(24) - after analysis: estimated number of entries in the factor matrices on the processor. If
    #    negative, then the absolute value corresponds to millions of entries in the factors. Note that in the
    #    unsymmetric case, INFO(24)=INFO(3). In the symmetric case, however, INFO(24) < INFO(3).
    #    The total number of entries in the factor matrices for all processors (sum of the INFO(24) values
    #    over all processors) is returned in INFOG(20).
    # === End MUMPS snippet ===

    @param(index=24, page=100)
    class param_24:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(25) page 100 from userguide_5.8.1.txt:5521-5522 ===
    # INFO(25) - after factorization: number of tiny pivots (number of pivots modified by static pivoting)
    #    detected on the processor (see INFOG(25) for the the total number of tiny pivots).
    # === End MUMPS snippet ===

    @param(index=25, page=100)
    class param_25:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(26) page 100 from userguide_5.8.1.txt:5523-5525 ===
    # INFO(26) - after solution: effective size in MegaBytes (millions of bytes) of all working space to run
    #    the solution phase. (The maximum and sum over all processors are returned in INFOG(30) and
    #    INFOG(31), respectively).
    # === End MUMPS snippet ===

    @param(index=26, page=100)
    class param_26:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(27) page 100 from userguide_5.8.1.txt:5526-5528 ===
    # INFO(27) - after factorization: effective number of entries in factor matrices assuming full-rank
    #    factorization has been performed. If negative, then the absolute value corresponds to millions of
    #    entries in the factors. Note that in case full-rank storage of factors (ICNTL(35)=0 or 3), we have
    # === End MUMPS snippet ===

    @param(index=27, page=100)
    class param_27:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(28) page 100 from userguide_5.8.1.txt:5531-5534 ===
    # INFO(28) - after factorization: effective number of entries in factors on the processor taking into
    #    account BLR compression. If negative, then the absolute value corresponds to millions of entries in
    #    the factors. It is equal to INFO(27) when BLR functionality (see ICNTL(35)) is not activated
    #    or leads to no compression.
    # === End MUMPS snippet ===

    @param(index=28, page=100)
    class param_28:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(29) page 100 from userguide_5.8.1.txt:5535-5539 ===
    # INFO(29) - after analysis: minimum estimated size of the main internal real/complex workarray S
    #    to run the numerical factorization in-core when factors are stored low-rank (ICNTL(35)=1,2).
    #    If negative, then the absolute value corresponds to millions of real/complex entries needed in this
    #    workarray. It is also the estimated minimum size of LWK USER in that case, if the user provides
    #    WK USER.
    # === End MUMPS snippet ===

    @param(index=29, page=100)
    class param_29:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(30) page 100 from userguide_5.8.1.txt:5540-5549 ===
    # INFO(30) and INFO(31) - after analysis: estimated size in MegaBytes (millions of bytes) of all
    #    working space to perform low-rank numerical phases (factorization/solve) with low-rank factors
    #    (ICNTL(35)=1,2) and estimated compression rate given by ICNTL(38).
    #       • —– (30) in-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(36) and INFOG(37).
    #       • —– (31) out-of-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(38) and INFOG(39).
    #    See also INFO(22) which provides the actual memory that was needed but only after
    #    factorization. Numerical pivoting difficulties and the effective compression of the factors (in case
    #    ICNTL(35)=1,2) typically impact the difference between estimated and effective memory.
    # === End MUMPS snippet ===

    @param(index=30, page=100)
    class param_30:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(31) page 100 from userguide_5.8.1.txt:5540-5549 ===
    # INFO(30) and INFO(31) - after analysis: estimated size in MegaBytes (millions of bytes) of all
    #    working space to perform low-rank numerical phases (factorization/solve) with low-rank factors
    #    (ICNTL(35)=1,2) and estimated compression rate given by ICNTL(38).
    #       • —– (30) in-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(36) and INFOG(37).
    #       • —– (31) out-of-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(38) and INFOG(39).
    #    See also INFO(22) which provides the actual memory that was needed but only after
    #    factorization. Numerical pivoting difficulties and the effective compression of the factors (in case
    #    ICNTL(35)=1,2) typically impact the difference between estimated and effective memory.
    # === End MUMPS snippet ===

    @param(index=31, page=100)
    class param_31:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(32) page 100 from userguide_5.8.1.txt:5550-5554 ===
    # INFO(32) - after analysis: minimum estimated size of the main internal real/complex workarray S
    #    to run the numerical factorization in-core when factors and contribution blocks are stored low-
    #    rank (ICNTL(35)=1,2 and ICNTL(37)=1). If negative, then the absolute value corresponds to
    #    millions of real/complex entries needed in this workarray. It is also the estimated minimum size of
    #    LWK USER in that case, if the user provides WK USER.
    # === End MUMPS snippet ===

    @param(index=32, page=100)
    class param_32:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(33) page 100 from userguide_5.8.1.txt:5555-5559 ===
    # INFO(33) - after analysis: minimum estimated size of the main internal real/complex workarray S to
    #    run the numerical factorization out-of-core when factors and contribution blocks are stored low-
    #    rank (ICNTL(35)=1,2 and ICNTL(37)=1). If negative, then the absolute value corresponds to
    #    millions of real/complex entries needed in this workarray. It is also the estimated minimum size of
    #    LWK USER in that case, if the user provides WK USER.
    # === End MUMPS snippet ===

    @param(index=33, page=100)
    class param_33:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(34) page 100 from userguide_5.8.1.txt:5560-5572 ===
    # INFO(34) and INFO(35) - after analysis: estimated size in MegaBytes (millions of bytes) of
    #    all working space to perform low-rank numerical phases (factorization/solve) with low-rank
    #    factors and low-rank contribution blocks (ICNTL(35)=1,2 and ICNTL(37)=1) and estimated
    #    compression rates given by ICNTL(38) and ICNTL(39) relatively.
    #         • —– (34) in-core factorization and solve. The maximum and sum over all processors are
    #           returned respectively in INFOG(40) and INFOG(41).
    #         • —– (35) out-of-core factorization and solve The maximum and sum over all processors are
    #           returned respectively in INFOG(42) and INFOG(43).
    #       See also INFO(22) which provides the actual memory that was needed but only after factorization.
    # === End MUMPS snippet ===

    @param(index=34, page=100)
    class param_34:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(35) page 100 from userguide_5.8.1.txt:5560-5572 ===
    # INFO(34) and INFO(35) - after analysis: estimated size in MegaBytes (millions of bytes) of
    #    all working space to perform low-rank numerical phases (factorization/solve) with low-rank
    #    factors and low-rank contribution blocks (ICNTL(35)=1,2 and ICNTL(37)=1) and estimated
    #    compression rates given by ICNTL(38) and ICNTL(39) relatively.
    #         • —– (34) in-core factorization and solve. The maximum and sum over all processors are
    #           returned respectively in INFOG(40) and INFOG(41).
    #         • —– (35) out-of-core factorization and solve The maximum and sum over all processors are
    #           returned respectively in INFOG(42) and INFOG(43).
    #       See also INFO(22) which provides the actual memory that was needed but only after factorization.
    # === End MUMPS snippet ===

    @param(index=35, page=100)
    class param_35:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(36) page 101 from userguide_5.8.1.txt:5573-5577 ===
    # INFO(36) - after analysis: minimum estimated size of the main internal real/complex workarray
    #    S to run the numerical factorization out-of-core when contribution blocks are stored low-rank
    #    (ICNTL(35)=0,3 and ICNTL(37)=1). If negative, then the absolute value corresponds to
    #    millions of real/complex entries needed in this workarray. It is also the estimated minimum size
    #    of LWK USER in that case, if the user provides WK USER.
    # === End MUMPS snippet ===

    @param(index=36, page=101)
    class param_36:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(37) page 101 from userguide_5.8.1.txt:5578-5586 ===
    # INFO(37) and INFO(38) - after analysis: estimated size in MegaBytes (millions of bytes) of
    #    all working space to perform low-rank numerical phases (factorization/solve) with low-rank
    #    contribution blocks only (ICNTL(35)=0,3 and ICNTL(37)=1) and estimated compression rate
    #    given by ICNTL(39).
    #        • —– (37) in-core factorization and solve. The maximum and sum over all processors are
    #          returned respectively in INFOG(44) and INFOG(45).
    #        • —– (38) out-of-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(46) and INFOG(47).
    #      See also INFO(22) which provides the actual memory that was needed but only after factorization.
    # === End MUMPS snippet ===

    @param(index=37, page=101)
    class param_37:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(38) page 101 from userguide_5.8.1.txt:5578-5586 ===
    # INFO(37) and INFO(38) - after analysis: estimated size in MegaBytes (millions of bytes) of
    #    all working space to perform low-rank numerical phases (factorization/solve) with low-rank
    #    contribution blocks only (ICNTL(35)=0,3 and ICNTL(37)=1) and estimated compression rate
    #    given by ICNTL(39).
    #        • —– (37) in-core factorization and solve. The maximum and sum over all processors are
    #          returned respectively in INFOG(44) and INFOG(45).
    #        • —– (38) out-of-core factorization and solve The maximum and sum over all processors are
    #          returned respectively in INFOG(46) and INFOG(47).
    #      See also INFO(22) which provides the actual memory that was needed but only after factorization.
    # === End MUMPS snippet ===

    @param(index=38, page=101)
    class param_38:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(39) page 101 from userguide_5.8.1.txt:5587-5589 ===
    # INFO(39) - after factorization: effective size of the main internal real/complex workarray S (allocated
    #    internally or by the user when WK USER is provided) to run the numerical factorization. If negative,
    #    then the absolute value corresponds to millions of real/complex entries needed in this workarray.
    # === End MUMPS snippet ===

    @param(index=39, page=101)
    class param_39:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: INFO(40) page 101 from userguide_5.8.1.txt:5592-5596 ===
    # INFO(40) is the number of negative pivots among the null pivots/deficiency detected. Note that,
    # for singular matrices, INFO(40) may vary from one run to another due to floating-point rounding
    # effects. A pivot counted in INFO(12), the number of negative non-null pivots, will not be counted
    # in INFO(40). See also INFOG(28) which provides the number of null pivots/deficiency over all
    # processors and INFOG(50), which provides the number of negative null pivots over all processors.
    # === End MUMPS snippet ===

    @param(index=40, page=101)
    class param_40:
        "Placeholder parameter, not processed yet."


# -------------
# RINFO members
# -------------


class RINFO(RawArray):
    """Floating-point info array RINFO(1..40)."""

    def __init__(self, size: int = LEN_RINFO):
        super().__init__(size, default=0.0)

    # === Begin MUMPS snippet: RINFO(1) page 97 from userguide_5.8.1.txt:5374-5375 ===
    # RINFO(1) - after analysis: The estimated number of floating-point operations on the processor for the
    #    elimination process.
    # === End MUMPS snippet ===

    @param(index=1, page=97)
    class param_1:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(2) page 97 from userguide_5.8.1.txt:5376-5377 ===
    # RINFO(2) - after factorization: The number of floating-point operations on the processor for the
    #    assembly process.
    # === End MUMPS snippet ===

    @param(index=2, page=97)
    class param_2:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(3) page 97 from userguide_5.8.1.txt:5378-5380 ===
    # RINFO(3) - after factorization: The number of floating-point operations on the processor for the
    #    elimination process. In case the BLR feature is activated (ICNTL(35)=1, 2 or 3), RINFO(3)
    #    represents the theoretical number of operations for the standard full-rank factorization.
    # === End MUMPS snippet ===

    @param(index=3, page=97)
    class param_3:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(4) page 97 from userguide_5.8.1.txt:5381-5384 ===
    # RINFO(4) - after factorization: The effective number of floating-point operations on the processor
    #    for the elimination process. It is equal to RINFO(3) when the BLR feature is not activated
    #    (ICNTL(35)=0) and will typically be smaller than RINFO(3) when the BLR feature is activated
    #    and leads to compression.
    # === End MUMPS snippet ===

    @param(index=4, page=97)
    class param_4:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(5) page 97 from userguide_5.8.1.txt:5389-5399 ===
    # RINFO(5) is computed for a full-rank factorization (ICNTL(35)=0 also for the factorization).
    # If ICNTL(35)=1, 2 or 3 at analysis, then RINFO(5) is computed assuming a low-rank (in-
    # core) storage of the factors of the BLR fronts during the factorization (ICNTL(35)=1 or 2 during
    # factorization). In case ICNTL(35)=1, 2 or 3 at analysis and the factors are stored in full-rank
    # format (ICNTL(35)=0 or 3 for the factorization), we refer the user to INFO(3) in order to obtain
    # a rough estimate of the necessary disk space for the concerned processor.
    #  The effective size in MegaBytes of the files written by the current processor will be returned in
    # === End MUMPS snippet ===

    @param(index=5, page=97)
    class param_5:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(6) page 98 from userguide_5.8.1.txt:5402-5404 ===
    # RINFO(6) - after factorization: in the case of an out-of-core execution (ICNTL(22)=1), the size in
    #    MegaBytes of the disk space used by the files written by the concerned processor is provided. The
    #    total disk space (for all processors) is returned in RINFOG(16).
    # === End MUMPS snippet ===

    @param(index=6, page=98)
    class param_6:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(7) page 98 from userguide_5.8.1.txt:5405-5406 ===
    # RINFO(7) - after each job: The size (in MegaBytes) of the file used to save the data on the processor
    #    (See Subsection 5.20).
    # === End MUMPS snippet ===

    @param(index=7, page=98)
    class param_7:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFO(8) page 98 from userguide_5.8.1.txt:5407-5407 ===
    # RINFO(8) - after each job: The size (in MegaBytes) of the MUMPS structure.
    # === End MUMPS snippet ===

    @param(index=8, page=98)
    class param_8:
        "Placeholder parameter, not processed yet."


# --------------
# RINFOG members
# --------------


class RINFOG(RawArray):
    """Floating-point global info array RINFOG(1..40)."""

    def __init__(self, size: int = LEN_RINFOG):
        super().__init__(size, default=0.0)

    # === Begin MUMPS snippet: RINFOG(1) page 101 from userguide_5.8.1.txt:5604-5605 ===
    # RINFOG(1) - after analysis: the estimated number of floating-point operations (on all processors) for
    #    the elimination process.
    # === End MUMPS snippet ===

    @param(index=1, page=101)
    class param_1:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(2) page 101 from userguide_5.8.1.txt:5606-5607 ===
    # RINFOG(2) - after factorization: the total number of floating-point operations (on all processors) for
    #    the assembly process.
    # === End MUMPS snippet ===

    @param(index=2, page=101)
    class param_2:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(3) page 101 from userguide_5.8.1.txt:5608-5610 ===
    # RINFOG(3) - after factorization: the total number of floating-point operations (on all processors) for
    #    the elimination process. In case the BLR feature is activated (ICNTL(35)=1, 2 or 3), RINFOG(3)
    #    represents the theoretical number of operations for the standard full-rank factorization.
    # === End MUMPS snippet ===

    @param(index=3, page=101)
    class param_3:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(4) page 101 from userguide_5.8.1.txt:5611-5612 ===
    # RINFOG(4) to RINFOG(8) - after solve with error analysis: Only returned if ICNTL(11) = 1 or 2.
    #    See description in Subsection 5.9 .
    # === End MUMPS snippet ===

    @param(index=4, page=101)
    class param_4:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(5) from userguide_5.8.1.txt:2287-2290 ===
    # RINFOG(4)), the infinite norm of the computed solution (∥x̄∥∞ in RINFOG(5)), and the
    #                   ∥Ax̄−b∥∞
    # scaled residual ∥A∥ ∞ ∥x̄∥∞
    #                             in RINFOG(6), a componentwise backward error estimate in
    # === End MUMPS snippet ===

    # === Begin MUMPS snippet: RINFOG(7) from userguide_5.8.1.txt:2291-2293 ===
    #   RINFOG(7) and RINFOG(8) are computed.
    # – If ICNTL(11)= 1, then in addition to the above statistics also an estimate for the error in
    #   the solution in RINFOG(9), and condition numbers for the linear system in RINFOG(10) and
    # === End MUMPS snippet ===

    # === Begin MUMPS snippet: RINFOG(8) page 101 from userguide_5.8.1.txt:5611-5612 ===
    # RINFOG(4) to RINFOG(8) - after solve with error analysis: Only returned if ICNTL(11) = 1 or 2.
    #    See description in Subsection 5.9 .
    # === End MUMPS snippet ===

    @param(index=8, page=101)
    class param_8:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(9) page 101 from userguide_5.8.1.txt:5613-5614 ===
    # RINFOG(9) to RINFOG(11) - after solve with error analysis: Only returned if ICNTL(11) = 1.
    #    See description in Subsection 5.9 .
    # === End MUMPS snippet ===

    @param(index=9, page=101)
    class param_9:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(11) page 101 from userguide_5.8.1.txt:5613-5614 ===
    # RINFOG(9) to RINFOG(11) - after solve with error analysis: Only returned if ICNTL(11) = 1.
    #    See description in Subsection 5.9 .
    # === End MUMPS snippet ===

    @param(index=11, page=101)
    class param_11:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(12) page 101 from userguide_5.8.1.txt:5615-5621 ===
    # RINFOG(12) - after factorization: if the computation of the determinant was requested (see
    #    ICNTL(33)), RINFOG(12) contains the real part of the determinant. The determinant may
    #    contain an imaginary part in case of complex arithmetic (see RINFOG(13)). It is obtained by
    #    multiplying (RINFOG(12), RINFOG(13)) by 2 to the power INFOG(34).
    # === End MUMPS snippet ===

    @param(index=12, page=101)
    class param_12:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(13) page 102 from userguide_5.8.1.txt:5623-5625 ===
    # RINFOG(13) - after factorization: if the computation of the determinant was requested (see
    #    ICNTL(33)), RINFOG(13) contains the imaginary part of the determinant. The determinant
    #    is then obtained by multiplying (RINFOG(12), RINFOG(13)) by 2 to the power INFOG(34).
    # === End MUMPS snippet ===

    @param(index=13, page=102)
    class param_13:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(14) page 102 from userguide_5.8.1.txt:5626-5629 ===
    # RINFOG(14) - after factorization: the total effective number of floating-point operations (on all
    #    processors) for the elimination process. It is equal to RINFOG(3) when the BLR feature is
    #    not activated (ICNTL(35)=0) and will typically be smaller than RINFOG(3) when the BLR
    #    functionality is activated and leads to compression.
    # === End MUMPS snippet ===

    @param(index=14, page=102)
    class param_14:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(15) page 102 from userguide_5.8.1.txt:5634-5640 ===
    # RINFOG(15) is computed for a full-rank factorization (ICNTL(35)=0 also for the factorization).
    # If ICNTL(35)=1, 2 or 3 at analysis, then RINFOG(15) is computed assuming a low-rank (in-
    # core) storage of the factors of the BLR fronts during the factorization (ICNTL(35)=2 during
    # factorization). In case ICNTL(35)=1, 2 or 3 for the analysis and the factors will be stored in full-
    # rank format (ICNTL(35)=0 or 3 for the factorization), we refer the user to INFOG(3) in order
    # to obtain a rough estimate of the necessary disk space for all processors.
    # The effective size in Megabytes of the files written by all processors will be returned in
    # === End MUMPS snippet ===

    @param(index=15, page=102)
    class param_15:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(16) page 102 from userguide_5.8.1.txt:5642-5643 ===
    # RINFOG(16) - after factorization: in the case of an out-of-core execution (ICNTL(22)=1), the total
    #    size in MegaBytes of the disk space used by the files written by all processors is provided.
    # === End MUMPS snippet ===

    @param(index=16, page=102)
    class param_16:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(17) page 102 from userguide_5.8.1.txt:5644-5645 ===
    # RINFOG(17) - after each job: sum over all processors of the sizes (in MegaBytes) of the files used to
    #    save the instance (See Subsection 5.20).
    # === End MUMPS snippet ===

    @param(index=17, page=102)
    class param_17:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(18) page 102 from userguide_5.8.1.txt:5646-5647 ===
    # RINFOG(18) - after each job: sum over all processors of the sizes (in MegaBytes) of the MUMPS
    #    structures.
    # === End MUMPS snippet ===

    @param(index=18, page=102)
    class param_18:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(19) page 102 from userguide_5.8.1.txt:5648-5650 ===
    # RINFOG(19) - after factorization: smallest pivot in absolute value selected during factorization of
    #    the preprocessed matrix Apre (see Equation (5)) and considering ALSO small pivots selected as
    #    null-pivots (see ICNTL(24)) and pivots on which static pivoting (see CNTL(4)) is effective.
    # === End MUMPS snippet ===

    @param(index=19, page=102)
    class param_19:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(20) page 102 from userguide_5.8.1.txt:5651-5655 ===
    # RINFOG(20) - after factorization: smallest pivot in absolute value selected during factorization of the
    #    preprocessed matrix Apre (see Equation (5)) and NOT considering small pivots selected as null-
    #    pivots (see ICNTL(24)) and pivots on which static pivoting (see CNTL(4)) is effective. A huge
    #    value of RINFOG(20) indicates that all pivots were either null-pivots or pivots on which static
    #    pivoting was performed.
    # === End MUMPS snippet ===

    @param(index=20, page=102)
    class param_20:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(21) page 102 from userguide_5.8.1.txt:5656-5657 ===
    # RINFOG(21) (experimental, subject to change in the future) - after factorization: largest pivot in
    #    absolute value selected during factorization of the preprocessed matrix Apre (see Equation (5)).
    # === End MUMPS snippet ===

    @param(index=21, page=102)
    class param_21:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(22) page 102 from userguide_5.8.1.txt:5658-5659 ===
    # RINFOG(22) - after factorization: total number of floating-point operations offloaded to the
    #    accelerator(s) by all MPI processes. See also RINFO(9).
    # === End MUMPS snippet ===

    @param(index=22, page=102)
    class param_22:
        "Placeholder parameter, not processed yet."

    # === Begin MUMPS snippet: RINFOG(23) page 102 from userguide_5.8.1.txt:5660-5661 ===
    # RINFOG(23) - after factorization: average (over all MPI processes) time spent in operations offloaded
    #    to the accelerator(s), including communication. See also RINFO(10).
    # === End MUMPS snippet ===

    @param(index=23, page=102)
    class param_23:
        "Placeholder parameter, not processed yet."


__all__ = [
    "ICNTL",
    "CNTL",
    "INFO",
    "RINFO",
    "RINFOG",
]
