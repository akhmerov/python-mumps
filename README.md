# mumpy
Python bindings for the MUMPS package MUMPS: a parallel sparse direct solver

# Installation
## On Linux
Either install `mumps` with:
```
sudo apt-get instal libmumps-scotch-dev
```

## Both on Linux and OSX
or use the `conda` package on `conda-forge`
```
conda config --add channels conda-forge
conda install cython mumps numpy scipy libgfortran
pip install git+git@github.com:basnijholt/mumpy.git
```
and use the following `build.conf` (if on Linux, not needed for OSX)

```
[mumps]
include_dirs = $CONDA_PREFIX/include
library_dirs = $CONDA_PREFIX/lib
libraries = zmumps mumps_common pord metis esmumps scotch scotcherr mpiseq gfortran
extra_link_args = -Wl,-rpath=$CONDA_PREFIX/lib
```

# Usage
```python
import scipy.sparse.linalg as sla
from scipy.sparse import identity
import mumpy.mumps as mumps


def sparse_diag(matrix, k, sigma, **kwargs):
    """Call sla.eigsh with mumps support.

    Please see scipy.sparse.linalg.eigsh for documentation.
    
    Notes
    -----
    mumpy only works with complex numbers at the moment.
    """
    class LuInv(sla.LinearOperator):
        def __init__(self, A):
            inst = mumps.MUMPSContext()
            inst.analyze(A, ordering='pord')
            inst.factor(A)
            self.solve = inst.solve
            sla.LinearOperator.__init__(self, A.dtype, A.shape)

        def _matvec(self, x):
            return self.solve(x.astype(self.dtype))

    opinv = LuInv(matrix - sigma * identity(matrix.shape[0]))
    return sla.eigsh(matrix, k, sigma=sigma, OPinv=opinv, **kwargs)
```
