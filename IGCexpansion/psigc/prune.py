"""Batched pruning over many patterns on a tree with one duplication node, and its adjoint.

Post-duplication nodes use an n_post-state space, pre-duplication nodes an n_pre-state space; the
edge from a pre-duplication parent into the duplication node sees the duplication node's partial
restricted to the `diag` states (both paralogs identical right after duplication).
"""
import numpy as np


class Pruner:
    def __init__(self, tree, diag):
        self.tree, self.diag = tree, diag

    def _child_partial(self, v, ch, part):
        cp = part[ch]
        if ch == self.tree.dup and not self.tree.is_post[v]:
            cp = cp[self.diag]
        return cp

    def forward(self, P, tips, prior, keep=False):
        """P: edge -> matrix; tips: tip name -> TipData; prior: root distribution.
        Returns L (U,) and, with keep=True, the cache needed by backward()."""
        t = self.tree
        part, up = {}, {}
        for v in t.postorder:
            kids = t.children[v]
            if not kids:
                continue
            prod = None
            for ch in kids:
                e = (v, ch)
                u = tips[ch].up(P[e]) if not t.children[ch] else P[e] @ self._child_partial(v, ch, part)
                if keep:
                    up[e] = u
                prod = u.copy() if prod is None else prod * u
            part[v] = prod
            if not keep:
                for ch in kids:
                    part.pop(ch, None)
        rp = part[t.root][self.diag] if t.root_by_dup else part[t.root]
        L = prior @ rp
        return (L, (part, up, rp)) if keep else (L, None)

    def backward(self, P, tips, prior, cache, omega):
        """Gradient of sum_x omega_x L_x with respect to every edge matrix and the root prior."""
        t = self.tree
        part, up, rp = cache
        G = {}
        g_prior = rp @ omega
        bar_root = prior[:, None] * omega[None, :]
        if t.root_by_dup:
            full = np.zeros_like(part[t.root])
            full[self.diag] = bar_root
            bar_root = full
        bar = {t.root: bar_root}
        for v in reversed(t.postorder):
            kids = t.children[v]
            if not kids:
                continue
            bp = bar.pop(v)
            for i, ch in enumerate(kids):
                bu = bp
                for j, oc in enumerate(kids):
                    if j != i:
                        bu = bu * up[(v, oc)]
                e = (v, ch)
                if not t.children[ch]:
                    G[e] = tips[ch].grad(bu)
                else:
                    cp = self._child_partial(v, ch, part)
                    G[e] = bu @ cp.T
                    bc = P[e].T @ bu
                    if ch == t.dup and not t.is_post[v]:
                        full = np.zeros_like(part[ch])
                        full[self.diag] = bc
                        bc = full
                    bar[ch] = bc
        return G, g_prior
