"""Block-generator oracle for parameter gradients (plan §6.2, tests only).

For an edge with generator Q, differential matrix D, length t, bottom vector p and top vector q,
e^{M t} [0; p] with M = [[Q, D], [0, Q]] is [ int_0^t e^{Q(t-u)} D e^{Qu} p du ; e^{Qt} p ], so
q'top / q'bottom is the exact (to action tolerance) Frechet term of d lnL_s / d theta. This uses
a separate 2S-state action instance and is far more expensive than Gauss-Legendre; it is the
reference GL is tested against.
"""

import numpy as np

from ._lib import NONE
from .instance import ActionInstance


def block_gradient(lib, bs, D, dpi=None):
    """d lnL / d theta_k for each k, with the same inputs as BeagleScene.gl_gradient."""
    cs = bs.cs
    S, P, E, G = cs.S, bs.P, cs.E, cs.G
    vecs = bs.edge_vectors()
    rp = bs.inst.get_partials(cs.idx[cs.root_b])[0]
    den_root = rp @ bs.distn
    w = bs.data.weights
    out = np.zeros(len(D))
    for k, Dk in enumerate(D):
        edges = [e for e in range(E) if Dk[cs.procs[e]] is not None and bs.rates[e] > 0]
        total = np.zeros(P)
        if edges:
            nb = 2 * len(edges) + 2
            inst = ActionInstance(lib, tips=1, partials=nb, states=2 * S, patterns=P, generators=G,
                                  keys=nb, scale_buffers=0)
            try:
                inst.set_categories([1.0], [1.0])
                inst.set_pattern_weights(w)
                for g in range(G):
                    r, c, v = bs.trips[g]
                    if Dk[g] is not None:
                        dr, dc, dv = Dk[g]
                    else:
                        dr, dc, dv = np.zeros(0, int), np.zeros(0, int), np.zeros(0)
                    inst.set_generator(g, np.concatenate([r, r + S, dr]), np.concatenate([c, c + S, dc + S]),
                                       np.concatenate([v, v, dv]))
                ONES = nb - 1
                inst.set_partials(ONES, np.ones(P * 2 * S))
                for i, e in enumerate(edges):
                    x = np.zeros((P, 2 * S))
                    x[:, S:] = vecs[e][1]
                    inst.set_partials(1 + 2 * i, x.ravel())
                for g in range(G):
                    idx = [1 + 2 * i for i, e in enumerate(edges) if cs.procs[e] == g]
                    lens = [bs.rates[e] for e in edges if cs.procs[e] == g]
                    if g == 0:
                        idx.append(ONES)
                        lens.append(0.0)
                    inst.update_transitions(g, idx, lens)
                inst.update_partials([(2 + 2 * i, NONE, NONE, 1 + 2 * i, 1 + 2 * i, ONES, ONES)
                                      for i in range(len(edges))])
                for i, e in enumerate(edges):
                    o = inst.get_partials(2 + 2 * i)[0]
                    q = vecs[e][0]
                    total += np.sum(q * o[:, :S], axis=1) / np.sum(q * o[:, S:], axis=1)
            finally:
                inst.close()
        if dpi is not None and dpi[k] is not None:
            total += (rp @ dpi[k]) / den_root
        out[k] = float(w @ total)
    return out
