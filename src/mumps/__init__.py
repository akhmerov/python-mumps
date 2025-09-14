# Retrieve the installed package version via PEP 566 metadata.
from importlib.metadata import PackageNotFoundError, version

from .mumps import (
    AnalysisStatistics,
    Context,
    FactorizationStatistics,
    MUMPSError,
    complex_to_real,
    nullspace,
    orderings,
    possible_orderings,
    schur_complement,
)
from .enums import JOB, PAR, SYM

try:
    __version__ = version("python-mumps")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = [
    "Context",
    "MUMPSError",
    "AnalysisStatistics",
    "FactorizationStatistics",
    "schur_complement",
    "nullspace",
    "possible_orderings",
    "orderings",
    "complex_to_real",
    "JOB",
    "PAR",
    "SYM",
    "__version__",
    "test",
]


def test(verbose=True):
    try:
        from pytest import main
    except ImportError:
        raise RuntimeError("pytest is required to run the tests")
    import os.path

    return main(
        [os.path.dirname(os.path.abspath(__file__)), "-s"] + (["-v"] if verbose else [])
    )


test.__test__ = False
