# Copyright 2011-2016 Kwant authors.
# Copyright 2018 Python-MUMPS Authors.
#
# This file is part of Python-MUMPS. It is subject to the license terms in the file
# LICENSE found in the top-level directory of this distribution. A list of
# Python-MUMPS authors can be found in the file AUTHORS.md at the top-level
# directory of this distribution and at
# https://gitlab.kwant-project.org/kwant/python-mumps.

import numpy as np


class _Random:
    def __init__(self):
        self._x = 13

    def _set_seed(self, seed):
        self._x = seed

    def _randf(self):
        # A very bad random number generator returning numbers between -1 and
        # +1.  Just for making some matrices, and being sure that they are the
        # same on any architecture.
        m = 2**16
        a = 11929
        c = 36491

        self._x = (a * self._x + c) % m

        return (float(self._x) / m - 0.5) * 2

    def _randi(self):
        # A very bad random number generator returning number between 0 and 20.
        # Just for making some matrices, and being sure that they are the same
        # on any architecture.
        m = 2**16
        a = 11929
        c = 36491

        self._x = (a * self._x + c) % m

        return self._x % 21

    def randmat(self, n, m, dtype):
        mat = np.empty((n, m), dtype=dtype)

        if issubclass(dtype, np.complexfloating):
            for i in range(n):
                for j in range(m):
                    mat[i, j] = self._randf() + 1j * self._randf()
        elif issubclass(dtype, np.floating):
            for i in range(n):
                for j in range(m):
                    mat[i, j] = self._randf()
        elif issubclass(dtype, np.integer):
            for i in range(n):
                for j in range(m):
                    mat[i, j] = self._randi()

        return mat

    def randvec(self, n, dtype):
        vec = np.empty(n, dtype=dtype)

        if issubclass(dtype, np.complexfloating):
            for i in range(n):
                vec[i] = self._randf() + 1j * self._randf()
        elif issubclass(dtype, np.floating):
            for i in range(n):
                vec[i] = self._randf()
        elif issubclass(dtype, np.integer):
            for i in range(n):
                vec[i] = self._randi()

        return vec
