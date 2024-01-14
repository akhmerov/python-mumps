# mumpy

Python bindings for the [MUMPS](http://mumps-solver.org/): a parallel sparse direct solver.

## Installation

`mumpy` works with Python 3.10 and higher on Linux, Windows and Mac.

The recommended way to install `mumpy` is using `mamba`/`conda`.

```bash
mamba install -c conda-forge mumpy
```

`mumpy` can also be installed from PyPI, however this is a more involved procedure
that requires separately installing the MUMPS library and a C compiler.

## Usage example

The following example shows how mumpy can be used to implement sparse diagonalization
with Scipy.

```python
import scipy.sparse.linalg as sla
from scipy.sparse import identity
import mumpy


def sparse_diag(matrix, k, sigma, **kwargs):
    """Call sla.eigsh with mumps support.

    See scipy.sparse.linalg.eigsh for documentation.
    """
    class LuInv(sla.LinearOperator):
        def __init__(self, A):
            inst = mumpy.Context()
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

`mumpy` recommends [Spin](https://github.com/scientific-python/spin/). Get spin with:

```bash
pip install spin
```

Then to build, test and install `mumpy`:

```bash
spin build
spin test -- --lf  # (Pytest arguments go after --)
spin install
```
