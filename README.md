# mumpy

Python bindings for the [MUMPS](http://mumps.enseeiht.fr/): a parallel sparse direct solver.


# Installation

`mumpy` works with Python 3.5 and higher on Linux, Windows and Mac.

The recommended way to install `mumpy` is using [`conda`](https://conda.io/):
```bash
conda install -c conda-forge mumpy
```

`mumpy` can also be installed from PyPI, however this is a more involved procedure
that requires separately installing the MUMPS library and a C compiler.


# Usage example

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
            inst = mumpy.MUMPSContext()
            inst.analyze(A, ordering='pord')
            inst.factor(A)
            self.solve = inst.solve
            sla.LinearOperator.__init__(self, A.dtype, A.shape)

        def _matvec(self, x):
            return self.solve(x.astype(self.dtype))

    opinv = LuInv(matrix - sigma * identity(matrix.shape[0]))
    return sla.eigsh(matrix, k, sigma=sigma, OPinv=opinv, **kwargs)
```
