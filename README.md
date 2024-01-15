# Python-MUMPS

Python bindings for the [MUMPS](http://mumps-solver.org/): a parallel sparse direct solver.

## Scope

This package targets MUMPS packaged by conda-forge using Cython bindings. It
aims to provide a full wrapper of the MUMPS sequential API. Its primary target
OS is Linux.

Next steps include:

- Support for Windows and OSX
- Support for distributed (MPI) MUMPS

## Installation

`python-mumps` works with Python 3.10 and higher on Linux, Windows and Mac.

The recommended way to install `python-mumps` is using `mamba`/`conda`.

```bash
mamba install -c conda-forge python-mumps
```

`python-mumps` can also be installed from PyPI, however this is a more involved procedure
that requires separately installing the MUMPS library and a C compiler.

## Usage example

The following example shows how Python-MUMPS can be used to implement sparse diagonalization
with Scipy.

```python
import scipy.sparse.linalg as sla
from scipy.sparse import identity
import mumps


def sparse_diag(matrix, k, sigma, **kwargs):
    """Call sla.eigsh with mumps support.

    See scipy.sparse.linalg.eigsh for documentation.
    """
    class LuInv(sla.LinearOperator):
        def __init__(self, A):
            inst = mumps.Context()
            inst.analyze(A, ordering='pord')
            inst.factor(A)
            self.solve = inst.solve
            sla.LinearOperator.__init__(self, A.dtype, A.shape)

        def _matvec(self, x):
            return self.solve(x.astype(self.dtype))

    opinv = LuInv(matrix - sigma * identity(matrix.shape[0]))
    return sla.eigsh(matrix, k, sigma=sigma, OPinv=opinv, **kwargs)
```

## Development

`python-mumps` recommends [Spin](https://github.com/scientific-python/spin/). Get spin with:

```bash
pip install spin
```

Then to build, test and install `python-mumps`:

```bash
spin build
spin test -- --lf  # (Pytest arguments go after --)
spin install
```
