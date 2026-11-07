"""Runtime shims that let jsonctmctree 0.2.0 import and run on NumPy >= 2 and SciPy >= 1.14.

jsonctmctree is used unmodified. Two names it relies on were removed upstream:

* ``numpy.testing.Tester`` (removed in NumPy 1.25): ``jsonctmctree/__init__.py`` imports it
  only to attach ``test``/``bench`` hooks, so a no-op stub is enough.
* ``scipy.sparse.spmatrix.A`` (removed in SciPy 1.14): the live call site is
  ``jsonctmctree/pyexp/linear_system.py`` (small n); the others are in classes IGCexpansion
  does not use. ``.A`` is restored as ``toarray()``. ``numpy.matrix.A`` still exists and
  needs nothing.

``install()`` only adds what is missing, so it is a no-op on old library versions, and it is
idempotent. It runs when ``IGCexpansion`` is imported, before any module imports jsonctmctree.
"""

import numpy.testing as _npt
import scipy.sparse as _sp

APPLIED = []


class _Tester:
    """Stand-in for the removed ``numpy.testing.Tester``; jsonctmctree never calls it."""

    def __init__(self, *args, **kwargs):
        pass

    def test(self, *args, **kwargs):
        pass

    def bench(self, *args, **kwargs):
        pass


def install():
    if not hasattr(_npt, "Tester"):
        _npt.Tester = _Tester
        APPLIED.append("numpy.testing.Tester")
    if not hasattr(_sp.spmatrix, "A"):
        _sp.spmatrix.A = property(lambda self: self.toarray())
        APPLIED.append("scipy.sparse.spmatrix.A")
    return list(APPLIED)


install()
