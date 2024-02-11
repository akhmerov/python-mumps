from enum import Enum


class IntEnum(int, Enum):
    def __new__(cls, value, description=None, validator=None):
        """IntEnum with description and validator values which can be Enum."""
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        obj.validator = validator
        return obj


class Orderings(IntEnum):
    AMD = 0
    USER_DEFINED = 1
    AMF = 2
    SCOTCH = 3
    PORD = 4
    METIS = 5
    QAMD = 6
    AUTO = 7


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
    STDERR = 0  # Fortran convention
    STDOUT = 6  # Fortran convention


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
    )

    MATRIX_INPUT_FORMAT = 5, "controls the matrix input format"

    PERMUTE_SCALE = (
        6,
        "permutes the matrix to a zero-free diagonal and/or scale the matrix",
    )

    ORDERING = 7, "computes a symmetric permutation in case of sequential analysis"

    SCALING_STRATEGY = 8, "_icntl_describes the scaling strategy"

    SOLVE_TRANSPOSED = 9, "computes the solution using A or AT"

    ITERATIVE_REFINEMENT = (
        10,
        "applies the iterative refinement to the computed solution",
    )

    ERROR_ANALYSIS = 11, "computes statistics related to an error analysis"

    ORDERING_SYMMETRIC = 12, "defines an ordering strategy for symmetric matrices"

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
    )

    SCHUR_COMPLEMENT = 19, "computes the Schur complement matrix"

    RHS_FORMAT = 20, (
        "determines the format (dense, sparse, or distributed) "
        "of the right-hand sides"
    )

    DISTRIBUTION = 21, (
        "determines the distribution (centralized or distributed) "
        "of the solution vectors"
    )

    OUT_OF_CORE = 22, "controls the in-core/out-of-core (OOC) factorization and solve"

    WORKING_MEMORY_SIZE_MB = (
        23,
        "corresponds to the maximum size of the working memory in MegaBytes that MUMPS",
    )

    DETECT_NULL_PIVOTS = 24, "controls the detection of “null pivot rows”"

    DEFICIENT_MATRIX = 25, (
        "allows the computation of a solution of a deficient matrix and also "
        "of a null space basis"
    )

    SCHUR_SOLUTION_TYPE = 26, "drives the solution phase if a Schur complement matrix"

    BLOCK_SIZE = 27, "controls the blocking size for multiple right-hand sides"

    ORDERING_COMPUTATION = 28, (
        "determines whether a sequential or parallel computation of the "
        "ordering is performed"
    )

    PARALLEL_ORDERING_TOOL = 29, (
        "defines the parallel ordering tool to be used to compute the fill-in "
        "reducing permutation"
    )

    COMPUTE_INVERSE_ENTRIES = 30, (
        "computes a user-specified set of entries in the inverse A^-1 of the "
        "original matrix"
    )

    DISCARD_FACTORS = (
        31,
        "indicates which factors may be discarded during the factorization",
    )

    FORWARD_ELIMINATION = 32, (
        "performs the forward elimination of the right-hand sides during the "
        "factorization"
    )

    DETERMINANT = 33, "computes the determinant of the input matrix"

    DELETE_FILES = 34, "controls the deletion of the files in case of save/restore"

    BLR = 35, "controls the activation of the Block Low-Rank (BLR) feature"

    BLR_VARIANT = 36, "controls the choice of BLR factorization variant"

    BLR_COMPRESSION = 37, "controls the BLR compression of the contribution blocks"

    LU_COMPRESSION_RATE = 38, "estimates compression rate of LU factors"

    BLOCK_COMPRESSION_RATE = 39, "estimates compression rate of contribution blocks"

    COMPACT_WORKARRAY_ID = (
        49,
        "compact workarray id%S at the end of factorization phase",
    )

    SYMBOLIC_FACTORIZATION = 58, "defines options for symbolic factorization"
