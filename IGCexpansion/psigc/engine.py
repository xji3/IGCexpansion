"""Composite likelihood and exact gradient for HKY + PS-IGC (and single-site HKY + IS-IGC).

Parameters use the legacy PSJSGeneconv x layout (see model.unpack_x). The gradient is exact up to the
Chebyshev interpolation error: edge-matrix adjoints from pruning, one adjoint Frechet derivative per
(class, edge, node), and the derivative in p through the interpolant (design document section 4).
"""
import time

import numpy as np
import scipy.linalg as sla

from . import cheb
from .model import (PS256, IS16, PRE16, PRE4, DIAG256, DIAG16, N_PM, I_ETA, I_P,
                    unpack_x, site_Q_jacobian, pi_jacobian)
from .prune import Pruner


def _frechet_adj(X, G):
    """M = L_{X^T}(G), so that <G, L_X(E)> = <M, E> for every E."""
    return sla.expm_frechet(X.T, G, compute_expm=False)


class _Grad:
    """Accumulator for d lnL / d (site matrices, tau, edges, pi, p via w)."""

    def __init__(self):
        self.H = {s: np.zeros((4, 4)) for s in (1, 2, 3)}
        self.dtau = 0.0
        self.dedge = {}
        self.dpi = np.zeros(4)
        self.dlogp_w = 0.0

    def edge(self, e, v):
        self.dedge[e] = self.dedge.get(e, 0.0) + v

    def to_x(self, x, par, tree):
        g = np.zeros(len(x))
        xpm = np.asarray(x[:N_PM], float)
        for s in (1, 2, 3):
            if np.any(self.H[s]):
                g[:N_PM] += np.einsum('jab,ab->j', site_Q_jacobian(xpm, s), self.H[s])
        g[:N_PM] += pi_jacobian(xpm) @ self.dpi
        dlogtau = par.tau * self.dtau
        g[I_ETA] = dlogtau                           # tau = eta / p
        g[I_P] = -dlogtau + self.dlogp_w
        g[I_P + 1:] = tree.dlogblen_to_x(self.dedge)
        return g


class ISEngine:
    """Single-site HKY + IS-IGC over the alignment columns (also used for far pairs)."""

    def __init__(self, ds):
        self.ds, self.tree = ds, ds.tree
        self.pruner = Pruner(self.tree, DIAG16)

    def _site(self, par, s, col_weights=None, acc=None):
        """lnL per alignment column of codon position s; with col_weights also accumulates the
        gradient of sum_col col_weights[col] * lnL_col into acc."""
        sb = self.ds.site_blocks[s]
        Q = par.site_Q(s)
        Gpost = IS16.build([Q], par.tau)
        P = {}
        for e in self.tree.post_edges:
            P[e] = sla.expm(par.blen[e] * Gpost)
        for e in self.tree.pre_edges:
            P[e] = sla.expm(par.blen[e] * Q)
        L, cache = self.pruner.forward(P, sb['tips'], par.pi, keep=col_weights is not None)
        lncol = np.log(L)[sb['inv']]
        if col_weights is not None:
            wu = np.bincount(sb['inv'], weights=col_weights, minlength=sb['U'])
            G, gp = self.pruner.backward(P, sb['tips'], par.pi, cache, wu / L)
            Mpost = np.zeros((16, 16)); Mpre = np.zeros((4, 4))
            for e, Ge in G.items():
                post = self.tree.is_post[e[0]]
                X = par.blen[e] * (Gpost if post else Q)
                M = _frechet_adj(X, Ge)
                acc.edge(e, float(np.sum(M * X)))
                if post:
                    Mpost += par.blen[e] * M
                else:
                    Mpre += par.blen[e] * M
            H, sA, _ = IS16.fold(Mpost)
            acc.H[s] += H[0]
            acc.dtau += sA
            acc.H[s] += PRE4.fold(Mpre)[0][0]
            acc.dpi += gp
        return lncol

    def site_lnl(self, par, col_weights=None, acc=None):
        out = np.zeros(self.ds.L)
        for s, sb in self.ds.site_blocks.items():
            cw = None if col_weights is None else col_weights[sb['cols']]
            out[sb['cols']] = self._site(par, s, cw, acc)
        return out

    def lnl(self, x):
        return float(self.site_lnl(unpack_x(x, self.tree)).sum())

    def lnl_grad(self, x):
        """Single-site IS lnL and its exact gradient in the legacy x layout (p enters only via tau)."""
        x = np.asarray(x, float)
        par = unpack_x(x, self.tree)
        acc = _Grad()
        f = float(self.site_lnl(par, np.ones(self.ds.L), acc).sum())
        return f, acc.to_x(x, par, self.tree)


class PSEngine:
    K_LADDER = (16, 24, 32, 48)
    TAIL_TOL = 1e-11

    def __init__(self, ds, K=20, w_tol=1e-15, exact_fallback=True):
        """K: fixed number of Chebyshev nodes, or 'auto': per class, move up the ladder 16/24/32/48 (starting
        at the task floor, which only rises until begin_task()) until the Chebyshev tail relative to each used
        pattern's smallest node value passes cheb.tail_min (rounding-aware threshold, 1e-11 to 1e-10) and no
        interpolated pair likelihood is non-positive; a class
        that fails at K = 48 is evaluated exactly (one node per distinct distance) when exact_fallback."""
        self.ds, self.tree, self.w_tol = ds, ds.tree, w_tol
        self.auto = (K == 'auto')
        self.K = self.K_LADDER[0] if self.auto else K
        self.exact_fallback = exact_fallback
        self.pruner = Pruner(self.tree, DIAG256)
        self.is_engine = ISEngine(ds)
        self.last = {}
        self._offsets, off = {}, 0
        for c, blk in ds.blocks.items():          # pair order of the `pairs` output: block by block
            self._offsets[c] = off
            off += len(blk.n)

    # ------------------------------------------------------------------ public API
    def begin_task(self):
        if self.auto:
            self.K = self.K_LADDER[0]

    def lnl(self, x, pairs=False):
        out = self._run(np.asarray(x, float), grad=False, pairs=pairs)
        return (out['lnL'], out['pairs']) if pairs else out['lnL']

    def lnl_grad(self, x):
        out = self._run(np.asarray(x, float), grad=True)
        return out['lnL'], out['grad']

    # ------------------------------------------------------------------ core
    def _run(self, x, grad, pairs=False):
        t_all = time.perf_counter()
        par = unpack_x(x, self.tree)
        p, tau = par.p, par.tau
        acc = _Grad() if grad else None
        timing = {'expm': 0.0, 'prune': 0.0, 'adjoint': 0.0, 'frechet': 0.0}
        ds = self.ds
        # far / near split
        lw = np.log1p(-p) if p < 1 else -np.inf
        near_masks, w_by_block = {}, {}
        far_cnt = np.zeros(ds.L)
        wmin, wmax = np.inf, -np.inf
        for c, blk in ds.blocks.items():
            w = np.exp(blk.n * lw) if p < 1 else np.zeros(len(blk.n))
            nm = w >= self.w_tol
            near_masks[c], w_by_block[c] = nm, w
            far_cnt += np.bincount(blk.pair_I[~nm], minlength=ds.L) + np.bincount(blk.pair_J[~nm], minlength=ds.L)
            if nm.any():
                wmin, wmax = min(wmin, w[nm].min()), max(wmax, w[nm].max())
        # far pairs: single-site IS likelihoods
        lnL_far = 0.0
        lncol = None
        if far_cnt.any():
            lncol = self.is_engine.site_lnl(par, far_cnt if grad else None, acc)
            lnL_far = float(far_cnt @ lncol)
        lnL_near = 0.0
        pair_out = np.full(ds.n_pairs, np.nan) if pairs else None
        tails = {}
        Ks, fallback = {}, []
        if np.isfinite(wmin):
            lo, hi = float(wmin), float(wmax)
            prior16 = np.outer(par.pi, par.pi).ravel()
            Qs = {s: par.site_Q(s) for s in (1, 2, 3)}
            ladder = [k for k in self.K_LADDER if k >= self.K] if self.auto else [self.K]
            for c, blk in ds.blocks.items():
                nm = near_masks[c]
                if not nm.any():
                    continue
                sa, sb = c
                A, B = PS256.build_AB([Qs[sa], Qs[sb]], tau)
                Gpre = PRE16.build([Qs[sa], Qs[sb]])
                Ppre = {e: sla.expm(par.blen[e] * Gpre) for e in self.tree.pre_edges}
                inv, n = blk.inv[nm], blk.n[nm]
                wn = w_by_block[c][nm]
                attempts = [(k, False) for k in ladder]
                if self.auto and self.exact_fallback:
                    attempts.append((None, True))
                for K, exact in attempts:
                    if exact:
                        n_u = np.unique(n)
                        K = len(n_u)
                        wk = np.exp(n_u * lw)
                        pos = np.full(int(n_u.max()) + 1, -1); pos[n_u] = np.arange(K)
                        Bn = np.zeros((len(n), K)); Bn[np.arange(len(n)), pos[n]] = 1.0
                    else:
                        if hi - lo < 1e-14:
                            K = 1
                        wk = cheb.nodes(K, lo, hi)[0] if K > 1 else np.array([hi])
                        Bn = cheb.basis(wn, K, lo, hi) if K > 1 else np.ones((nm.sum(), 1))
                    Ps, Lk = [], np.empty((blk.U, K))
                    for k in range(K):
                        t0 = time.perf_counter()
                        Gk = A + wk[k] * B
                        P = dict(Ppre)
                        for e in self.tree.post_edges:
                            P[e] = sla.expm(par.blen[e] * Gk)
                        timing['expm'] += time.perf_counter() - t0
                        t0 = time.perf_counter()
                        Lk[:, k] = self.pruner.forward(P, blk.tips, prior16)[0]
                        timing['prune'] += time.perf_counter() - t0
                        if grad:
                            Ps.append((Gk, P))
                    Lpair = np.einsum('pk,pk->p', Lk[inv], Bn)
                    if exact:
                        tails[str(c)] = 0.0
                        fallback.append(str(c))
                        break
                    used = np.unique(inv)
                    tails[str(c)] = cheb.tail_min(Lk[used]) if K > 2 else 0.0
                    ok = tails[str(c)] <= 1.0 and (Lpair > 0).all()   # tail_min: tail / threshold
                    if ok or not self.auto or K == 1:
                        break
                Ks[str(c)] = K
                if self.auto and not exact:
                    self.K = max(self.K, K)                    # ratchet: K only rises within a task
                if (Lpair <= 0).any():
                    raise FloatingPointError('non-positive interpolated pair likelihood; increase K')
                lnL_near += float(np.log(Lpair).sum())
                if pairs:
                    idx_all = np.flatnonzero(nm)
                    pair_out[self._offsets[c] + idx_all] = np.log(Lpair)
                if grad:
                    t0 = time.perf_counter()
                    # p through w_n inside the interpolant (exact mode: through the nodes, below)
                    dwk_dlogp = -n_u * p * wk / (1 - p) if exact else None
                    if K > 1 and not exact:
                        Bd = cheb.basis_deriv(wn, K, lo, hi)
                        dLdw = np.einsum('pk,pk->p', Lk[inv], Bd)
                        dw_dlogp = -n * p * wn / (1 - p)
                        acc.dlogp_w += float(np.sum(dLdw / Lpair * dw_dlogp))
                    Om = np.zeros((blk.U, K))
                    np.add.at(Om, inv, Bn / Lpair[:, None])
                    MA = np.zeros((256, 256)); MB = np.zeros((256, 256))
                    Gpre_acc = {e: 0.0 for e in self.tree.pre_edges}
                    gprior = np.zeros(16)
                    timing['adjoint'] += time.perf_counter() - t0
                    for k in range(K):
                        Gk, P = Ps[k]
                        t0 = time.perf_counter()
                        _, cache = self.pruner.forward(P, blk.tips, prior16, keep=True)
                        Gm, gp = self.pruner.backward(P, blk.tips, prior16, cache, Om[:, k])
                        del cache
                        gprior += gp
                        timing['adjoint'] += time.perf_counter() - t0
                        t0 = time.perf_counter()
                        Mk = np.zeros((256, 256))
                        for e in self.tree.post_edges:
                            X = par.blen[e] * Gk
                            M = _frechet_adj(X, Gm[e])
                            acc.edge(e, float(np.sum(M * X)))
                            Mk += par.blen[e] * M
                        for e in self.tree.pre_edges:
                            Gpre_acc[e] = Gpre_acc[e] + Gm[e]
                        MA += Mk
                        MB += wk[k] * Mk
                        if exact:
                            acc.dlogp_w += float(np.sum(Mk * B)) * dwk_dlogp[k]
                        timing['frechet'] += time.perf_counter() - t0
                    H, sA, _ = PS256.fold(MA)
                    _, _, sB = PS256.fold(MB)
                    acc.H[sa] += H[0]; acc.H[sb] += H[1]
                    acc.dtau += sA + sB
                    Mpre = np.zeros((16, 16))
                    for e in self.tree.pre_edges:
                        X = par.blen[e] * Gpre
                        M = _frechet_adj(X, Gpre_acc[e])
                        acc.edge(e, float(np.sum(M * X)))
                        Mpre += par.blen[e] * M
                    Hp = PRE16.fold(Mpre)[0]
                    acc.H[sa] += Hp[0]; acc.H[sb] += Hp[1]
                    P4 = gprior.reshape(4, 4)
                    acc.dpi += P4 @ par.pi + P4.T @ par.pi
                    del Ps
        if pairs:
            for c, blk in ds.blocks.items():
                idx = np.flatnonzero(~near_masks[c])
                if len(idx):
                    pair_out[self._offsets[c] + idx] = lncol[blk.pair_I[idx]] + lncol[blk.pair_J[idx]]
        lnL = lnL_far + lnL_near
        timing['total'] = time.perf_counter() - t_all
        self.last = {'timing': timing, 'K': self.K, 'K_used': Ks, 'fallback': fallback,
                     'interval': (wmin, wmax), 'cheb_tail': tails,
                     'far_pairs': int(sum((~m).sum() for m in near_masks.values()))}
        out = {'lnL': lnL}
        if grad:
            out['grad'] = acc.to_x(x, par, self.tree)
        if pairs:
            out['pairs'] = pair_out
        return out
