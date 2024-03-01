from enum import Enum


class IntEnum(int, Enum):
    def __new__(cls, value, description=None, values=None):
        """IntEnum with description and values which can be Enum."""
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        obj.values = values
        return obj


class Jobs(IntEnum):
    INITIALIZE = -1
    TERMINATE = -2
    REMOVE_SAVED_DATA = -3
    ANALYZE = 1
    FACTORIZE = 2
    SOLVE = 3
    ANALYZE_FACTORIZE = 4
    FACTORIZE_SOLVE = 5
    ANALYZE_FACTORIZE_SOLVE = 6
    SAVE = 7
    RESTORE = 8


class Stream(IntEnum):
    NULL = -1
    STDOUT = 6  # Fortran convention


class PrintLevels(IntEnum):
    NO_MESSAGES = 0, "No messages output"
    ERROR_MESSAGES = 1, "Only error messages printed"
    ERRORS_WARNINGS_STATISTICS = 2, "Errors, warnings, and main statistics printed"
    ERRORS_WARNINGS_DIAGNOSTICS = 3, "Errors and warnings and terse diagnostics printed"
    ERRORS_WARNINGS_INFORMATION = (
        4,
        "Errors, warnings and information on input, output parameters printed",
    )


class MatrixInputFormats(IntEnum):
    ASSEMBLED = 0, (
        "The matrix must be input in the structure components N, NNZ (or NZ), IRN, JCN,"
        " and A if the matrix is centralized on the host (see Subsection 5.2.2.1) or "
        "in the structure components N, NNZ loc (or NZ loc), IRN loc, JCN loc, A loc if"
        " the matrix is distributed (see Subsection 5.2.2.3)."
    )
    ELEMENTAL = 1, (
        "The matrix must be input in the structure components N, NELT, ELTPTR, ELTVAR, "
        "and A ELT (see Subsection 5.2.2.3)."
    )


class PermuteScale(IntEnum):
    NO_COLUMN_PERMUTATION = 0, "No column permutation is computed."
    MAXIMIZE_DIAGONAL_ENTRIES = 1, (
        "The permuted matrix has as many entries on its diagonal as possible. "
        "The values on the diagonal are of arbitrary size."
    )
    MAXIMIZE_MIN_DIAGONAL_ENTRY = 2, (
        "The permutation is such that the smallest value on the diagonal of the "
        "permuted matrix is maximized. The numerical values of the original matrix, "
        "(mumps par%A), must be provided by the user during the analysis phase."
    )
    VARIANT_MAXIMIZE_MIN_DIAGONAL_ENTRY = 3, (
        "Variant of option 2 with different performance. The numerical values of the "
        "original matrix (mumps par%A) must be provided by the user during the analysis"
        " phase."
    )
    MAXIMIZE_SUM_DIAGONAL_ENTRIES = 4, (
        "The sum of the diagonal entries of the permuted matrix is maximized. The "
        "numerical values of the original matrix (mumps par%A) must be provided by the "
        "user during the analysis phase."
    )
    MAXIMIZE_PRODUCT_DIAGONAL_ENTRIES = 5, (
        "The product of the diagonal entries of the permuted matrix is maximized. "
        "Scaling vectors are also computed and stored in COLSCA and ROWSCA, if "
        "ICNTL(8) is set to -2 or 77. With these scaling vectors, the nonzero diagonal "
        "entries in the permuted matrix are one in absolute value and all the "
        "off-diagonal entries less than or equal to one in absolute value. For "
        "unsymmetric matrices, COLSCA and ROWSCA are meaningful on the permuted matrix "
        "A Qc (see Equation (5)). For symmetric matrices, COLSCA and ROWSCA are "
        "meaningful on the original matrix A. The numerical values of the original "
        "matrix, mumps par%A, must be provided by the user during the analysis phase."
    )
    VARIANT_PRODUCT_DIAGONAL_ENTRIES = 6, (
        "Similar to 5 but with a more costly (time and memory footprint) algorithm. "
        "The numerical values of the original matrix, mumps par%A, must be provided by "
        "the user during the analysis phase."
    )
    AUTO = 7, (
        "Based on the structural symmetry of the input matrix and on the availability "
        "of the numerical values, the value of ICNTL(6) is automatically chosen by the "
        "software."
    )


class Orderings(IntEnum):
    AMD = 0
    USER_DEFINED = 1
    AMF = 2
    SCOTCH = 3
    PORD = 4
    METIS = 5
    QAMD = 6
    AUTO = 7


class ScalingStrategies(IntEnum):
    ANALYSIS = -2, (
        "Scaling done during analysis (see [24, 25] for the unsymmetric case and [26] "
        "for the symmetric case), under certain conditions. The original matrix has to "
        "be centralized (ICNTL(18)=0) and the user has to provide its numerical values "
        "(mumps par%A) on entry to the analysis, otherwise the scaling will not be "
        "computed. Also, ICNTL(6) must be activated for the scaling to be computed. The"
        " effective value of the scaling applied is reported in INFOG(33). We recommend"
        " that the user checks for INFOG(33) in order to know if the scaling was "
        "applied. When not applied, other scalings can be performed during the "
        "factorization."
    )
    USER_PROVIDED = -1, (
        "Scaling provided by the user. Scaling arrays must be provided in COLSCA and "
        "ROWSCA on entry to the numerical factorization phase by the user, who is then "
        "responsible for allocating and freeing them. If the input matrix is symmetric "
        "(SYM= 1 or 2), then the user should ensure that the array ROWSCA is equal to "
        "(or points to the same location as) the array COLSCA."
    )
    NO_SCALING = 0, "No scaling applied/computed."
    DIAGONAL = 1, "Diagonal scaling computed during the numerical factorization phase."
    COLUMN = 3, "Column scaling computed during the numerical factorization phase."
    ROW_COLUMN_INF_NORMS = 4, (
        "Row and column scaling based on infinite row/column norms, computed during the"
        " numerical factorization phase."
    )
    SIMULTANEOUS_ROW_COLUMN = 7, (
        "Simultaneous row and column iterative scaling (based on [45, 16, 36, 35]) "
        "computed during the numerical factorization phase."
    )
    RIGOROUS_SIMULTANEOUS_ROW_COLUMN = 8, (
        "Similar to 7 but more rigorous and expensive to compute; computed during the "
        "numerical factorization phase."
    )
    AUTOMATIC = 77, "Automatic choice of the value of ICNTL(8) done during analysis."


class SolveTransposed(IntEnum):
    TRANSPOSED = 0
    NOT_TRANSPOSED = 1


class ErrorAnalysis(IntEnum):
    NO_STATISTICS = 0, "No error analysis is performed (no statistics)."
    MAIN_STATISTICS = 1, ("Compute all the statistics (very expensive).")
    COMPUTE_MAIN_STATISTICS = 2, (
        "Compute main statistics (norms, residuals, componentwise backward errors), "
        "but not the most expensive ones like (condition number and forward error "
        "estimates)."
    )


class OrderingsSymmetric(IntEnum):
    AUTO = 0, "Automatic choice"
    USUAL = 1, "Usual ordering (nothing done)"
    COMPRESSED_GRAPH = 2, "Ordering on the compressed graph associated with the matrix"
    CONSTRAINED_ORDERING = 3, (
        "Constrained ordering, only available with AMF (ICNTL(7)=2)."
    )


class DistributedInputStrategies(IntEnum):
    CENTRALIZED = 0, (
        "The input matrix is centralized on the host (see Subsection 5.2.2.1)."
    )
    MAPPING = 1, (
        "The user provides the structure of the matrix on the host at analysis, MUMPS "
        "returns a mapping and the user should then provide the matrix entries "
        "distributed according to the mapping on entry to the numerical factorization "
        "phase (see Subsection 5.2.2.2)."
    )
    DISTRIBUTED_ENTRIES = 2, (
        "The user provides the structure of the matrix on the host at analysis, and the"
        " distributed matrix entries on all slave processors at factorization. Any "
        "distribution is allowed (see Subsection 5.2.2.2)."
    )
    DIRECTLY_PROVIDED = 3, (
        "User directly provides the distributed matrix, pattern and entries, input both"
        " for analysis and factorization (see Subsection 5.2.2.2)."
    )


class SchurComplement(IntEnum):
    NO_SCHUR = 0, "Complete factorization. No Schur complement is returned."
    CENTRALIZED = 1, (
        "The Schur complement matrix will be returned centralized by rows on the host "
        "after the factorization phase. On the host before the analysis phase, the user"
        " must set the integer variable SIZE SCHUR to the size of the Schur matrix, the"
        " integer pointer array LISTVAR SCHUR to the list of indices of the Schur "
        "matrix."
    )
    DISTRIBUTED = 2, (
        "The Schur complement matrix will be returned distributed by columns: the Schur"
        " will be returned on the slave processors in the form of a 2D block cyclic "
        "distributed matrix (ScaLAPACK style) after factorization. Workspace should be "
        "allocated by the user before the factorization phase in order for MUMPS to "
        "store the Schur complement (see SCHUR, SCHUR MLOC, SCHUR NLOC, and SCHUR LLD "
        "in Subsection 5.15). On the host before the analysis phase, the user must set "
        "the integer variable SIZE SCHUR to the size of the Schur matrix, the integer "
        "pointer array LISTVAR SCHUR to the list of indices of the Schur matrix. The "
        "integer variables NPROW, NPCOL, MBLOCK, NBLOCK may also be defined (default "
        "values will otherwise be provided)."
    )


class RHSFormats(IntEnum):
    DENSE = 0, (
        "The right-hand side is in dense format in the structure component RHS, NRHS, "
        "LRHS (see Subsection 5.14.1)"
    )
    SPARSE_AUTO = 1, (
        "The decision of exploiting sparsity of the right-hand side to accelerate the "
        "solution phase is done automatically."
    )
    SPARSE_NO = 2, (
        "Sparsity of the right-hand side is NOT exploited to improve solution phase."
    )
    SPARSE_YES = 3, (
        "Sparsity of the right-hand side is exploited during solution phase."
    )
    DISTRIBUTED = 10, (
        "The right-hand side is provided distributed in the structure components Nloc "
        "RHS, LRHS loc, IRHS loc, RHS loc (see Subsection 5.14.3)."
    )
    DISTRIBUTED_FILL = 11, (
        "Fill IRHS loc to match the distribution of the solution (imposed by MUMPS), "
        "in case of distributed solution (ICNTL(21)=1)."
    )


class SolutionDistributions(IntEnum):
    CENTRALIZED = 0, (
        "The solution vector is assembled and stored in the structure component RHS "
        "(gather phase), that must have been allocated earlier by the user (see "
        "Subsection 5.14.5)."
    )
    DISTRIBUTED = 1, (
        "The solution vector is kept distributed on each slave processor in the "
        "structure components ISOL loc and SOL loc. ISOL loc and SOL loc must then "
        "have been allocated by the user and must be of size at least INFO(23), where "
        "INFO(23) has been returned by MUMPS at the end of the factorization phase "
        "(see Subsection 5.14.6)."
    )


class OutOfCore(IntEnum):
    IN_CORE = 0, "In-core factorization and solution phases (default standard version)."
    OUT_OF_CORE = 1, (
        "Out-of-core factorization and solve phases. The complete matrix of factors is "
        "written to disk (see Subsection 3.14)."
    )


class DetectNullPivots(IntEnum):
    DONT_DETECT = 0, (
        "Nothing done. A null pivot row will result in error INFO(1)=-10."
    )
    DETECT = 1, "Null pivot row detection."


class SchurSolutionType(IntEnum):
    STANDARD = 0, (
        "Standard solution phase on the internal problem; referring to the notations "
        "from Subsection 3.17, only the system A1,1 x1 = b1 is solved and the entries "
        "of the right-hand side corresponding to the Schur are explicitly set to 0 on "
        "output."
    )
    CONDENSE = 1, (
        "Condense/reduce the right-hand side on the Schur. Only a forward elimination "
        "is performed. The solution corresponding to the ‘internal” (non-Schur) "
        "variables is returned together with the reduced/condensed right-hand-side. "
        "The reduced right-hand side is made available on the host in the pointer "
        "array REDRHS, that must be allocated by the user. Its leading dimension "
        "LREDRHS must be provide, too."
    )
    EXPAND = 2, (
        "Expand the Schur local solution on the complete solution variables. REDRHS is "
        "considered to be the solution corresponding to the Schur variables. It must "
        "be allocated by the user as well as its leading dimension LREDRHS must be "
        "provided. The backward substitution is then performed with the given right-"
        "hand side to compute the solution associated with the ”internal” variables. "
        "Note that the solution corresponding to the Schur variables is also made "
        "available in the main solution vector/matrix."
    )


class OrderingParallelism(IntEnum):
    AUTO = 0, "Automatic choice."
    SEQUENTIAL = 1, (
        "Sequential computation. In this case the ordering method is set by ICNTL(7) "
        "and the ICNTL(29) parameter is meaningless (choice of the parallel ordering "
        "tool)."
    )
    PARALLEL = 2, (
        "Parallel computation. A parallel ordering and parallel symbolic factorization "
        "is requested by the user. For that, one of the parallel ordering tools (or "
        "all) must be available, and the matrix should not be too small. The ordering "
        "method is set by ICNTL(29) and the ICNTL(7) parameter is meaningless."
    )


class ParallelOrderingTool(IntEnum):
    AUTO = 0, "Automatic choice."
    PT_SCOTCH = 1, "PT-SCOTCH is used to reorder the input matrix, if available."
    PAR_METIS = 2, "ParMetis is used to reorder the input matrix, if available."


class ComputeInverseEntries(IntEnum):
    NO_ENTRIES = 0, "No entries in A^-1 are computed."
    COMPUTE_ENTRIES = 1, "Computes entries in A^-1."


class DiscardFactors(IntEnum):
    KEEP_FACTORS = 0, (
        "The factors are kept during the factorization, phase except in the case of "
        "unsymmetric matrices when the forward elimination is performed during "
        "factorization (ICNTL(32) = 1). In this case, since it will not be used during "
        "the solve phase, the L factor is discarded."
    )
    DISCARD_ALL = 1, (
        "All factors are discarded during the factorization phase. The user is not "
        "interested in solving the linear system (Equations (3) or (4)) and will not "
        "call MUMPS solution phase (JOB=3). This option is meaningful when only a Schur"
        " complement is needed (see ICNTL(19)), or when only statistics from the "
        "factorization, such as (for example) definiteness, value of the determinant, "
        "number of entries in factors after numerical pivoting, number of negative or "
        "null pivots are required. In this case, the memory allocated for the "
        "factorization will rely on the out-of-core estimates (and factors will not be "
        "written to disk)."
    )
    KEEP_U = 2, (
        "This setting is meaningful only for unsymmetric matrices and has no impact on "
        "symmetric matrices: only the U factor is kept after factorization so that "
        "exclusively a backward substitution is possible during the solve phase "
        "(JOB=3). This can be useful when: the user is only interested in the "
        "computation of a null space basis (see ICNTL(25)) during the solve phase, "
        "or the forward elimination is performed during the factorization "
        "(ICNTL(32)=1). Note that for unsymmetric matrices, if the forward elimination "
        "is performed during the factorization (ICNTL(32) = 1) then the L factor is "
        "always discarded during factorization. In this case (ICNTL(32) = 1), both "
        "ICNTL(31) = 0 and ICNTL(31) = 2 have the same behaviour."
    )


class ForwardElimination(IntEnum):
    NO_ELIMINATION = 0, ("Standard factorization not involving right-hand sides.")
    ELIMINATION = 1, (
        "Forward elimination (Equation (3)) of the right-hand side vectors is performed"
        " during factorization (JOB=2). The solve phase (JOB=3) will then only involve "
        "backward substitution (Equation (4))."
    )


class Determinant(IntEnum):
    NO_DETERMINANT = 0, "The determinant of the input matrix is not computed."
    COMPUTE_DETERMINANT = 1, (
        "Computes the determinant of the input matrix. The determinant is obtained by "
        "computing (a + ib) × 2c where a =RINFOG(12), b =RINFOG(13) and c = INFOG(34). "
        "In real arithmetic b=RINFOG(13) is equal to 0."
    )


class DeleteFiles(IntEnum):
    DELETE = 0, "The out-of-core files are marked out for deletion"
    NO_DELETE = 1, (
        "The out-of-core files should not be deleted because another saved instance "
        "references them."
    )


class BLR(IntEnum):
    STANDARD = 0, "Standard analysis and factorization (BLR feature is not activated)."
    AUTO = 1, (
        "BLR feature is activated and automatic choice of BLR option is performed by "
        "the software."
    )
    FACTORIZE_SOLVE = 2, (
        "BLR feature is activated during both the factorization and solution phases, "
        "which allows for memory gains by storing the factors in low-rank."
    )
    FACTORIZE = 3, (
        "BLR feature is activated during the factorization phase but not the solution "
        "phase, which is still performed in full-rank. As a consequence, the full-rank "
        "factors must be kept and no memory gains can be obtained. In an OOC context, "
        "(ICNTL(22)=1) this option enables the user to write all factors to disk which "
        "is not the case with ICNTL(35)=2 since factors in low-rank form are not "
        "written to disk."
    )


class BLRVariant(IntEnum):
    UFSC = 0, ("Standard UFSC variant with low-rank updates accumulation (LUA)")
    UCFS = 1, (
        "UCFS variant with low-rank updates accumulation (LUA). This variant consists "
        "in performing the compression earlier in order to further reduce the number "
        "of operations. Although it may have a numerical impact, the current "
        "implementation is still compatible with numerical pivoting."
    )


class BLRCompression(IntEnum):
    NO_COMPRESSION = 0, "Contribution blocks are not compressed"
    COMPRESSION = 1, (
        "Contribution blocks are compressed, reducing the memory consumption at the "
        "cost of some additional operations"
    )


class CompactWorkarrayId(IntEnum):
    NOTHING = 0, "Nothing is done."
    SATISFY_MEMORY = 1, (
        "Compact workarray id%S(MAXS) at the end of the factorization phase while "
        "satisfying the memory constraint that might have been provided with ICNTL(23) "
        "feature."
    )
    NO_MEMORY_CONSTRAINT = 2, (
        "Compact workarray id%S(MAXS) at the end of the factorization phase. The memory"
        " constraint that might have been provided with ICNTL(23) feature does not "
        "apply to this process."
    )


class SymbolicFactorization(IntEnum):
    QUOTIENT_GRAPH = 1, (
        "Symbolic factorization based on quotient graph, mixing right looking and left "
        "looking updates"
    )
    COLUMN_COUNT = 2, ("Column count based symbolic factorization based on [30]")


class ICNTL(IntEnum):

    ERROR_STREAM = 1, "is the output stream for error messages", Stream

    DIAGNOSTIC_STREAM = (
        2,
        (
            "is the output stream for diagnostic printing, statistics, and"
            " warning message"
        ),
        Stream,
    )

    GLOBAL_STREAM = (
        3,
        "is the output stream for global information, collected on the host",
        Stream,
    )

    PRINT_LEVEL = (
        4,
        "is the level of printing for error, warning, and diagnostic messages",
        PrintLevels,
    )

    MATRIX_INPUT_FORMAT = 5, "controls the matrix input format", MatrixInputFormats

    PERMUTE_SCALE = (
        6,
        "permutes the matrix to a zero-free diagonal and/or scale the matrix",
        PermuteScale,
    )

    ORDERING = (
        7,
        "computes a symmetric permutation in case of sequential analysis",
        Orderings,
    )

    SCALING_STRATEGY = 8, "describes the scaling strategy", ScalingStrategies

    SOLVE_TRANSPOSED = 9, "computes the solution using A or AT", SolveTransposed

    ITERATIVE_REFINEMENT = (
        10,
        "applies the iterative refinement to the computed solution",
    )

    ERROR_ANALYSIS = (
        11,
        "computes statistics related to an error analysis",
        ErrorAnalysis,
    )

    ORDERING_SYMMETRIC = (
        12,
        "defines an ordering strategy for symmetric matrices",
        OrderingsSymmetric,
    )

    ROOT_PARALLELISM = 13, "controls the parallelism of the root node"

    WORKING_SPACE_PERCENTAGE = (
        14,
        "controls the percentage increase in the estimated working space",
    )

    EXPLOIT_COMPRESSION = (
        15,
        "exploits compression of the input matrix resulting from a block format",
    )

    OMP_NUM_THREADS = 16, "controls the setting of the number of OpenMP threads"

    DISTRIBUTED_INPUT_STRATEGY = (
        18,
        "defines the strategy for the distributed input matrix",
        DistributedInputStrategies,
    )

    SCHUR_COMPLEMENT = 19, "computes the Schur complement matrix", SchurComplement

    RHS_FORMAT = (
        20,
        (
            "determines the format (dense, sparse, or distributed) "
            "of the right-hand sides"
        ),
        RHSFormats,
    )

    SOLUTION_DISTRIBUTION = (
        21,
        (
            "determines the distribution (centralized or distributed) "
            "of the solution vectors"
        ),
        SolutionDistributions,
    )

    OUT_OF_CORE = (
        22,
        "controls the in-core/out-of-core (OOC) factorization and solve",
        OutOfCore,
    )

    WORKING_MEMORY_SIZE_MB = (
        23,
        "corresponds to the maximum size of the working memory in MegaBytes that MUMPS",
    )

    DETECT_NULL_PIVOTS = (
        24,
        "controls the detection of “null pivot rows”",
        DetectNullPivots,
    )

    DEFICIENT_MATRIX = 25, (
        "allows the computation of a solution of a deficient matrix and also "
        "of a null space basis"
    )

    SCHUR_SOLUTION_TYPE = (
        26,
        "drives the solution phase if a Schur complement matrix",
        SchurSolutionType,
    )

    BLOCK_SIZE = 27, "controls the blocking size for multiple right-hand sides"

    ORDERING_PARALLELISM = (
        28,
        (
            "determines whether a sequential or parallel computation of the "
            "ordering is performed"
        ),
        OrderingParallelism,
    )

    PARALLEL_ORDERING_TOOL = (
        29,
        (
            "defines the parallel ordering tool to be used to compute the fill-in "
            "reducing permutation"
        ),
        ParallelOrderingTool,
    )

    COMPUTE_INVERSE_ENTRIES = (
        30,
        (
            "computes a user-specified set of entries in the inverse A^-1 of the "
            "original matrix"
        ),
        ComputeInverseEntries,
    )

    DISCARD_FACTORS = (
        31,
        "indicates which factors may be discarded during the factorization",
        DiscardFactors,
    )

    FORWARD_ELIMINATION = (
        32,
        (
            "performs the forward elimination of the right-hand sides during the "
            "factorization"
        ),
        ForwardElimination,
    )

    DETERMINANT = 33, "computes the determinant of the input matrix", Determinant

    DELETE_FILES = (
        34,
        "controls the deletion of the files in case of save/restore",
        DeleteFiles,
    )

    BLR = 35, "controls the activation of the Block Low-Rank (BLR) feature", BLR

    BLR_VARIANT = 36, "controls the choice of BLR factorization variant", BLRVariant

    BLR_COMPRESSION = (
        37,
        "controls the BLR compression of the contribution blocks",
        BLRCompression,
    )

    LU_COMPRESSION_RATE = 38, "estimates compression rate of LU factors"

    BLOCK_COMPRESSION_RATE = 39, "estimates compression rate of contribution blocks"

    COMPACT_WORKARRAY_ID = (
        49,
        "compact workarray id%S at the end of factorization phase",
        CompactWorkarrayId,
    )

    SYMBOLIC_FACTORIZATION = (
        58,
        "defines options for symbolic factorization",
        SymbolicFactorization,
    )
