"""dQ/dtheta and dpi/dtheta for the gradient (plan §6.2).

``scene_fd_jacobian`` works for any model that can produce a jsonctmctree scene at x: central
differences (h = 1e-5) of the generator and root prior. It is built on the full listed pattern of
each process, so entries that are zero at x (e.g. a rate forced to 0) still get their derivative.
Model-specific analytic Jacobians (e.g. PS-IGC in psjs.py) replace it where available.
"""

import numpy as np
import scipy.sparse as sp

from .scene import generator_triplets, scene_values


def _csr(trip, S):
    r, c, v = trip
    return sp.csr_matrix((v, (r, c)), shape=(S, S))


def triplets_of(M):
    """Triplets of a sparse matrix with explicit zeros removed (None if empty)."""
    M = sp.coo_matrix(M)
    keep = M.data != 0.0
    if not np.any(keep):
        return None
    return (np.ascontiguousarray(M.row[keep], dtype=np.int32), np.ascontiguousarray(M.col[keep], dtype=np.int32),
            np.ascontiguousarray(M.data[keep], dtype=np.float64))


def scene_fd_jacobian(scene_at, x, free, cs, h=1e-5):
    """D[k][g] (triplets or None) and dpi[k] (array or None) for x[free[k]], by central differences.

    ``scene_at(x)`` must return the jsonctmctree scene at parameter vector x (and may change model
    state; the caller restores it). ``cs`` is the compiled scene (shape, S).
    """
    x = np.asarray(x, dtype=float)
    S = cs.S
    D, dpi = [], []
    for i in free:
        mats, pis = [], []
        for sgn in (+1, -1):
            xx = x.copy()
            xx[i] += sgn * h
            sc = scene_at(xx)
            mats.append([_csr(generator_triplets(cs.shape, p), S) for p in sc["process_definitions"]])
            pis.append(scene_values(sc, cs)[1])
        D.append([triplets_of((a - b) / (2 * h)) for a, b in zip(*mats)])
        d = (pis[0] - pis[1]) / (2 * h)
        dpi.append(d if np.any(d != 0) else None)
    scene_at(x)
    return D, dpi
