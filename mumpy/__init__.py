from .mumpy import (
    MUMPSContext,
    MUMPSError,
    AnalysisStatistics,
    FactorizationStatistics,
    schur_complement,
    nullspace,
    possible_orderings,
)

__all__ = [
    "MUMPSContext",
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


# Infer the version. Prioritize trying to get the version from git, but if that
# fails, fall back to the version in the package metadata.
try:
    # Setuptools_scm is not a dependency, and this is mainly for convenience of
    # developers.
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
except (ImportError, LookupError):
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("mumpy")
    except PackageNotFoundError:
        __version__ = "0.0.0+unknown"
