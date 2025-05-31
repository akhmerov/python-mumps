from .mumps import (
    Context,
    MUMPSError,
    AnalysisStatistics,
    FactorizationStatistics,
    schur_complement,
    nullspace,
    possible_orderings,
    orderings,
    complex_to_real,
)
from ._version import __version__

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
