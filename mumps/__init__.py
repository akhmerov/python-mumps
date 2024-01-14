from .mumps import (
    Context,
    MUMPSError,
    AnalysisStatistics,
    FactorizationStatistics,
    schur_complement,
    nullspace,
    possible_orderings,
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
    "__version__",
    "test",
]


def test(verbose=True):
    from pytest import main
    import os.path

    return main(
        [os.path.dirname(os.path.abspath(__file__)), "-s"] + (["-v"] if verbose else [])
    )


test.__test__ = False
