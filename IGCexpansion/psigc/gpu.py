"""psigc.gpu: the HKY + PS-IGC composite likelihood and its exact gradient as one float64 PyTorch program
(GPU backend; also runs on the CPU for tests). Same inputs and outputs as psigc.PSEngine (legacy x layout,
per-pair lnL in psigc block order), so the CPU engine is its reference.

    eng = GPUEngine(ds, device='cuda:0')     # K='auto': adaptive Chebyshev K with a per-task ratchet
    eng.begin_task()                         # start of an M2 task: K may fall back to 16
    f = eng.lnl(x)
    f, g = eng.lnl_grad(x)
    f, per_pair = eng.lnl(x, pairs=True)
    lp, S = eng.pair_scores(x)               # per-pair lnL and d lnL_ab / dx (forward mode)
    lp = eng.pair_lnl_direct(x)              # exact per-pair lnL, no interpolation (oracle)
    eng.last, eng.eval_log                   # per-evaluation diagnostics (K, tails, fallback, ...)

Design (PSIGC_GPU_PLAN.md section 4):
* generators by one GEMM against fixed dense basis tensors; Q_ck = A_c + w_k B;
* transition matrices by uniformization, exp(tQ) = sum_j Pois(j; lam t) (I + Q/lam)^j, with the powers
  shared by all edges (UnifSeries: hand-written reverse and forward rules); lam is held fixed, which is exact
  for the series; d/dt is exact through a zero-valued (t - t.detach()) Q P term, avoiding cancellation in
  the derivative of the Poisson weights; (m, squarings) come from a host-side plan (no device sync);
* pruning per class with the K nodes as batch, site-repeat compression (each internal node holds the unique
  patterns of the tips below it), one-hot GEMM adjoints for tip gathers, the duplication node in its 16
  diagonal rows only;
* interpolation on unique (pattern, distance) keys with counts; the basis is evaluated per key and its
  derivative in w goes through the Chebyshev differentiation matrix (B D), which is stable near nodes;
* far pairs (w_n < w_tol) through the single-site IS model with per-column far counts looked up from a
  suffix table, so every shape is static: one evaluation (lnL, gradient, diagnostics) is one CUDA-graph
  replay, keyed by the plan (K, m, squarings);
* K ladder 16/24/32/48 with the Chebyshev tail relative to each used pattern's SMALLEST node value, tested
  against a rounding-aware threshold (tail_threshold); K only rises within a task; a failing call is replayed at the next rung before it returns;
  classes that still fail at K = 48 are evaluated exactly (one node per distinct distance) and flagged.
"""
import math
import time

import numpy as np
import torch

from . import model as M
from .model import PS256, IS16, PRE16, DIAG256, DIAG16, unpack_x

DT = torch.float64
EPS = float(np.finfo(float).eps)
K_LADDER = (16, 24, 32, 48)
TAIL_TOL = 1e-11                 # tail test, see tail_threshold()
NOISE_C = 64.0
ERR_CAP = 1e-10
UNIF_TOL = 1e-17
MU_CAP = 7.0
M_ROUND = 4
NX_EXTRA = 3                     # packed buffer: x, then lo, hi, ncut

_TR = np.zeros((4, 4))
for _a, _b in M.TRANSITIONS:
    _TR[_a, _b] = 1.0
_OFF = 1.0 - np.eye(4)


# ============================================================================ host-side helpers
def spec_tensors(spec):
    """Dense basis of a generator family: G = sum_h Qflat[h] Epm[h] + tau (EA + w EB); zero row sums."""
    n = spec.n
    Epm = np.zeros((spec.nslot * 16, n, n))
    np.add.at(Epm, (spec.pm_flat_h, spec.pm_src, spec.pm_dst), 1.0)
    np.add.at(Epm, (spec.pm_flat_h, spec.pm_src, spec.pm_src), -1.0)
    EA, EB = np.zeros((n, n)), np.zeros((n, n))
    if len(spec.ig_src):
        np.add.at(EA, (spec.ig_src, spec.ig_dst), spec.ig_cA)
        np.add.at(EA, (spec.ig_src, spec.ig_src), -spec.ig_cA)
        np.add.at(EB, (spec.ig_src, spec.ig_dst), spec.ig_cB)
        np.add.at(EB, (spec.ig_src, spec.ig_src), -spec.ig_cB)
    return Epm, EA, EB


def poisson_m(mu, tol=UNIF_TOL, m_min=4):
    """Smallest m with P(N > m) < tol for N ~ Poisson(mu)."""
    if mu <= 0:
        return m_min
    term = math.exp(-mu)
    j = 0
    while True:
        j += 1
        term *= mu / j
        if j + 2 > mu:
            tail = term * (mu / (j + 1)) / (1 - mu / (j + 2))
            if tail < tol and j >= m_min:
                return j


def unif_plan(mu, tol=UNIF_TOL, mu_cap=MU_CAP, round_to=M_ROUND):
    s = 0 if mu <= mu_cap else int(math.ceil(math.log2(mu / mu_cap)))
    m = poisson_m(mu / 2 ** s, tol)
    return int(round_to * math.ceil(m / round_to)), s


def cheb_tables(K):
    k = np.arange(K)
    x = np.cos((2 * k + 1) * np.pi / (2 * K))
    lam = (-1.0) ** k * np.sin((2 * k + 1) * np.pi / (2 * K))
    T = np.cos(np.pi * np.outer(np.arange(K), 2 * k + 1) / (2 * K)) * (2.0 / K)
    T[0] /= 2
    X = x[:, None] - x[None, :]
    np.fill_diagonal(X, 1.0)
    D = (lam[None, :] / lam[:, None]) / X
    np.fill_diagonal(D, 0.0)
    np.fill_diagonal(D, -D.sum(1))
    return x, lam, T, D                  # D in x units on [-1, 1]


def tail_threshold(ratio, tol=TAIL_TOL):
    """Acceptance threshold for a pattern's min-normalised Chebyshev tail: max(tol, min(NOISE_C eps r, ERR_CAP)),
    r = max_k L_k / min_k L_k over the nodes. Once the interpolant has converged in w, the tail stops falling at a
    rounding floor of about 10-45 eps r (gpu2/results/tail_floor.json: 72 points, 4 arm-L datasets, tau x 1-16);
    unconverged tails were >= 659 eps r. A fixed 1e-11 sits on that floor for patterns with r ~ 5e3 and escalated
    K on noise. The cap keeps the certified per-pair error <= 1e-10: where the floor itself is higher (tau x 16),
    the test fails and the exact fallback is used."""
    return np.maximum(tol, np.minimum(NOISE_C * np.finfo(float).eps * ratio, ERR_CAP))


def _tip_keys(td):
    if td.idx is not None:
        return td.idx
    _, inv = np.unique(td.mask.T, axis=0, return_inverse=True)
    return inv.ravel()


# ============================================================================ autograd functions
class UnifSeries(torch.autograd.Function):
    """P[g, e] = sum_j wts[g, e, j] Pm[g]^j (wts not differentiated). Returns (P, Pw) with Pw the stored
    powers (non-differentiable). Reverse: A_m = Wb_m, A_j = Wb_j + A_{j+1} Pm^T, dPm = sum_j (Pm^{j-1})^T A_j
    (one GEMM with inner dimension m n). Forward mode: D_1 = dPm, D_j = D_{j-1} Pm + Pm^{j-1} dPm."""
    generate_vmap_rule = True

    @staticmethod
    def forward(Pm, wts):
        G, n, _ = Pm.shape
        m = wts.shape[-1] - 1
        Pw = torch.empty(G, m + 1, n, n, dtype=Pm.dtype, device=Pm.device)
        Pw[:, 0] = torch.eye(n, dtype=Pm.dtype, device=Pm.device)
        Pw[:, 1] = Pm
        for j in range(2, m + 1):
            torch.bmm(Pw[:, j - 1], Pm, out=Pw[:, j])
        P = torch.bmm(wts, Pw.view(G, m + 1, n * n)).view(G, wts.shape[1], n, n)
        return P, Pw

    @staticmethod
    def setup_context(ctx, inputs, output):
        Pm, wts = inputs
        P, Pw = output
        ctx.mark_non_differentiable(Pw)
        ctx.save_for_backward(Pm, wts, Pw)
        ctx.save_for_forward(Pm, wts, Pw)

    @staticmethod
    def backward(ctx, Pbar, _Pwbar):
        Pm, wts, Pw = ctx.saved_tensors
        G, E, n, _ = Pbar.shape
        m = wts.shape[-1] - 1
        A = torch.bmm(wts.transpose(1, 2), Pbar.reshape(G, E, n * n)).view(G, m + 1, n, n)
        PmT = Pm.transpose(1, 2)
        for j in range(m - 1, 0, -1):
            A[:, j] = torch.baddbmm(A[:, j], A[:, j + 1], PmT)
        Pmbar = torch.bmm(Pw[:, :m].reshape(G, m * n, n).transpose(1, 2), A[:, 1:].reshape(G, m * n, n))
        return Pmbar, None

    @staticmethod
    def jvp(ctx, dPm, _dwts):
        Pm, wts, Pw = ctx.saved_tensors
        G, n, _ = Pm.shape
        m = wts.shape[-1] - 1
        D = dPm
        acc = wts[:, :, 1, None, None] * D[:, None]
        for j in range(2, m + 1):
            D = torch.bmm(D, Pm) + torch.bmm(Pw[:, j - 1], dPm)
            acc = acc + wts[:, :, j, None, None] * D[:, None]
        return acc, None


SCORE_K_EXTRA = (64, 96)          # pair_scores: node counts tried when the ladder needs the exact fallback
SCORE_PAIR_TOL = 1e-10            # pair_scores: largest per-pair lnL error accepted against the exact values


def poisson_weights(mu, m):
    """Pois(j; mu), j = 0..m, by w_j = w_{j-1} mu / j (no cumprod: its backward syncs). (...) -> (..., m+1)."""
    w = [torch.exp(-mu)]
    for j in range(1, m + 1):
        w.append(w[-1] * (mu / j))
    return torch.stack(w, -1)


def expm_unif(Q, t, m, s):
    """exp(t_e Q_g) for generators Q (G, n, n) and times t (E,) -> (G, E, n, n); uniformization with powers
    shared by all E times, s squarings of t / 2^s. Differentiable in Q (reverse and forward mode) and exactly
    in t (zero-valued Q P term)."""
    G, n, _ = Q.shape
    E = t.shape[0]
    Qd = Q.detach()
    lam = (-torch.diagonal(Qd, dim1=-2, dim2=-1)).amax(-1).clamp_min(1e-30)
    Pm = torch.eye(n, dtype=Q.dtype, device=Q.device) + Q / lam[:, None, None]
    mu = lam[:, None] * t.detach()[None, :] * (0.5 ** s)
    P, _ = UnifSeries.apply(Pm, poisson_weights(mu, m))
    if s:
        P = P.reshape(G * E, n, n)
        for _ in range(s):
            P = torch.bmm(P, P)
        P = P.reshape(G, E, n, n)
    # d exp(tQ)/dt = Q exp(tQ), carried by a term whose value is exactly zero
    QP = torch.matmul(Qd[:, None], P.detach())
    return P + (t - t.detach())[None, :, None, None] * QP


class BaryBasis(torch.autograd.Function):
    """Barycentric basis B (N, K) at w; dB/dw = BD (Chebyshev differentiation matrix), not autograd of the
    barycentric quotient (whose error grows like eps / distance to the nearest node)."""
    generate_vmap_rule = True

    @staticmethod
    def forward(w, B, BD):
        return B.clone()

    @staticmethod
    def setup_context(ctx, inputs, output):
        ctx.save_for_backward(inputs[2])
        ctx.save_for_forward(inputs[2])

    @staticmethod
    def backward(ctx, gB):
        BD, = ctx.saved_tensors
        return (gB * BD).sum(-1), None, None

    @staticmethod
    def jvp(ctx, dw, _dB, _dBD):
        BD, = ctx.saved_tensors
        return BD * dw[..., None]


class ColGather(torch.autograd.Function):
    """P[..., idx] with a dense one-hot GEMM as reverse rule (no atomics: deterministic)."""
    generate_vmap_rule = True

    @staticmethod
    def forward(P, idx, onehot):
        return P.index_select(-1, idx)

    @staticmethod
    def setup_context(ctx, inputs, output):
        ctx.save_for_backward(inputs[2])
        ctx.idx = inputs[1]

    @staticmethod
    def backward(ctx, g):
        onehot, = ctx.saved_tensors
        return g @ onehot, None, None

    @staticmethod
    def jvp(ctx, dP, _di, _do):
        return dP.index_select(-1, ctx.idx)


class SortedGather(torch.autograd.Function):
    """X.index_select(dim, idx) whose reverse rule sums rows in a fixed order (deterministic): the target
    rows are visited in sorted order and summed with segment_reduce. Used for site-repeat maps and keys."""
    generate_vmap_rule = True

    @staticmethod
    def forward(X, idx, dim, perm, lengths, nrow):
        return X.index_select(dim, idx)

    @staticmethod
    def setup_context(ctx, inputs, output):
        X, idx, dim, perm, lengths, nrow = inputs
        ctx.dim, ctx.nrow = dim, nrow
        ctx.idx = idx
        ctx.save_for_backward(perm, lengths)

    @staticmethod
    def backward(ctx, g):
        perm, lengths = ctx.saved_tensors
        gs = g.index_select(ctx.dim, perm).movedim(ctx.dim, 0)
        red = torch.segment_reduce(gs, 'sum', lengths=lengths, axis=0, unsafe=True)
        return red.movedim(0, ctx.dim), None, None, None, None, None

    @staticmethod
    def jvp(ctx, dX, *_):
        return dX.index_select(ctx.dim, ctx.idx)


# ============================================================================ plan
class Plan:
    """Host-side plan of one evaluation: K, the near/far cut and the uniformization (m, s) per family.
    `key` identifies a CUDA graph; lo, hi, ncut travel in the packed input buffer."""
    __slots__ = ('K', 'near', 'lo', 'hi', 'ncut', 'ps', 'pre', 'is16', 'is4', 'mu_ps')

    def key(self, grad, parts):
        """Graph family: everything that changes the program except the numbers of series terms m."""
        return (self.K, self.near, self.ps[1], self.pre[1], self.is16[1], self.is4[1], bool(grad), parts)

    def ms(self):
        return (self.ps[0], self.pre[0], self.is16[0], self.is4[0])

    def with_ms(self, ms):
        q = Plan()
        for k in self.__slots__:
            setattr(q, k, getattr(self, k))
        q.ps, q.pre, q.is16, q.is4 = ((ms[0], self.ps[1]), (ms[1], self.pre[1]), (ms[2], self.is16[1]),
                                      (ms[3], self.is4[1]))
        return q


# ============================================================================ engine
class GPUEngine:
    def __init__(self, ds, K='auto', w_tol=1e-15, device='cuda:0', graphs=None, deterministic=True,
                 K_ladder=K_LADDER, tail_tol=TAIL_TOL, exact_fallback=True, exact_chunk=48,
                 parts='auto', log_evals=True, graph_mem_gb=12.0, max_patterns=None):
        """K: 'auto' (ladder with ratchet) or a fixed int. graphs: CUDA graphs (default: on for CUDA).
        deterministic: fixed-order reductions for every gather adjoint (bitwise-repeatable gradients; +3.5% on
        the A100, gpu2/results/bench_gpu.json), the default.
        graph_mem_gb: if the reserved GPU memory exceeds this before a new capture, the graph cache is cleared
        (bounded memory per worker; None disables).
        max_patterns: pattern chunking for large datasets (Dataset.chunk_blocks): no class block holds more
        than this many unique pair patterns; each chunk is evaluated as its own part (parts forced on), so
        memory scales with the chunk, not the dataset. Classes are then keyed (pos_a, pos_b, chunk).
        parts: 'auto' builds one graph up to K = 32 (<= 9.6 GB on the real data) and one per class (and the far
        part) above (K = 48: 6.1 GB instead of 14 GB, +7% time), or True/False."""
        self.ds, self.tree, self.w_tol = ds, ds.tree, float(w_tol)
        self.dev = dev = torch.device(device)
        self.auto = (K == 'auto')
        self.K_ladder = tuple(K_ladder)
        self.K_fixed = None if self.auto else int(K)
        self.K_floor = self.K_ladder[0] if self.auto else int(K)
        self.tail_tol = tail_tol
        self.exact_fallback, self.exact_chunk = exact_fallback, exact_chunk
        self.graphs = (dev.type == 'cuda') if graphs is None else bool(graphs)
        self.deterministic = deterministic
        self.parts_mode = True if max_patterns else parts
        self.max_patterns = max_patterns
        self.log_evals = log_evals
        self.last, self.eval_log = {}, []
        self._graphs, self._pool, self._cap_stream = {}, None, None
        self.n_capture, self.n_graph_clears = 0, 0
        self.graph_mem_gb = graph_mem_gb
        t = self.tree
        self.nx = 8 + t.n_rate
        T = lambda a, dt=DT: torch.as_tensor(np.asarray(a), dtype=dt, device=dev)
        L = lambda a: torch.as_tensor(np.asarray(a), dtype=torch.long, device=dev)
        self._T, self._L = T, L
        # generator bases (dense) and their diagonals (host, for the uniformization plan)
        Epm, EA, EB = spec_tensors(PS256)
        self.Epm, self.EA, self.EB = T(Epm.reshape(32, -1)), T(EA.reshape(-1)), T(EB.reshape(-1))
        self.dEpm, self.dEA, self.dEB = Epm.diagonal(0, 1, 2).copy(), EA.diagonal().copy(), EB.diagonal().copy()
        Epre, _, _ = spec_tensors(PRE16)
        self.Epre, self.dEpre = T(Epre.reshape(32, -1)), Epre.diagonal(0, 1, 2).copy()
        Eis, EAis, _ = spec_tensors(IS16)
        self.Eis, self.EAis = T(Eis.reshape(16, -1)), T(EAis.reshape(-1))
        self.dEis, self.dEAis = Eis.diagonal(0, 1, 2).copy(), EAis.diagonal().copy()
        self.TR, self.OFF = T(_TR), T(_OFF)
        self.diag256, self.diag16 = L(DIAG256), L(DIAG16)
        self.ie_post = {e: i for i, e in enumerate(t.post_edges)}
        self.ie_pre = {e: i for i, e in enumerate(t.pre_edges)}
        # classes: repeat maps, unique (pattern, distance) keys with counts, pair -> key
        blocks = ds.chunk_blocks(max_patterns) if max_patterns else ds.blocks
        self.classes = list(blocks.keys())
        self.blk = []
        off = 0
        for c in self.classes:
            b = blocks[c]
            keys, kinv, cnt = np.unique(np.stack([b.inv, b.n], 1), axis=0, return_inverse=True, return_counts=True)
            ku = keys[:, 0]
            d = dict(c=c, U=b.U, rep=self._repeat_maps(b.tips, b.U), off=off, npair=len(b.n),
                     key_u=L(ku), key_n=T(keys[:, 1]), key_n_np=keys[:, 1], cnt=T(cnt),
                     pair_key=L(kinv.ravel()), pair_n=T(b.n), I=L(b.pair_I), J=L(b.pair_J))
            d['key_gather'] = self._sorted_gather_meta(ku, b.U)
            self.blk.append(d)
            off += len(b.n)
        self.n_pairs = off
        # single-site blocks (IS model, far pairs)
        self.sites = []
        for s, sb in ds.site_blocks.items():
            self.sites.append(dict(s=s, rep=self._repeat_maps(sb['tips'], sb['U']), inv=L(sb['inv'])))
        cols = np.concatenate([ds.site_blocks[s]['cols'] for s in ds.site_blocks])
        self.col_perm = L(np.argsort(cols))
        # far counts: suf[col, j] = number of pairs of col with distance >= j
        H = np.zeros((ds.L, ds.n_max + 2))
        for b in ds.blocks.values():
            np.add.at(H, (b.pair_I, b.n), 1.0)
            np.add.at(H, (b.pair_J, b.n), 1.0)
        self.suf = T(np.cumsum(H[:, ::-1], 1)[:, ::-1].copy())
        self.n_present = np.unique(ds.pair_n_all)
        self._cheb = {}

    # ------------------------------------------------------------------ data structures
    def _sorted_gather_meta(self, idx, nrow):
        """perm (stable sort of idx) and segment lengths for SortedGather; rows with no entry get length 0."""
        idx = np.asarray(idx)
        perm = np.argsort(idx, kind='stable')
        lengths = np.bincount(idx, minlength=nrow)
        return self._L(idx), self._L(perm), self._L(lengths), int(nrow)

    def _repeat_maps(self, tips, U):
        """For every internal node v: the unique patterns S_v of the tips below v, tip columns on S_v, and for
        internal children the index of S_v in S_child. The root's set is the full pattern set (reordered)."""
        t = self.tree
        keys = {k: _tip_keys(td) for k, td in tips.items()}
        below, first, inv = {}, {}, {}
        for v in t.postorder:
            if not t.children[v]:
                below[v] = [v]
                continue
            below[v] = [x for c in t.children[v] for x in below[c]]
            kk = np.stack([keys[x] for x in below[v]], 1)
            _, f, iv = np.unique(kk, axis=0, return_index=True, return_inverse=True)
            first[v], inv[v] = f, iv.ravel()
        out = {}
        for v in t.postorder:
            if not t.children[v]:
                continue
            d = {'tips': {}, 'maps': {}, 'U': len(first[v])}
            for ch in t.children[v]:
                if not t.children[ch]:
                    td = tips[ch]
                    if td.idx is not None:
                        ii = td.idx[first[v]]
                        oh = torch.zeros(len(ii), td.n, dtype=DT, device=self.dev)
                        oh[torch.arange(len(ii), device=self.dev), self._L(ii)] = 1.0
                        d['tips'][ch] = ('idx', self._L(ii), oh)
                    else:
                        d['tips'][ch] = ('mask', self._T(td.mask[:, first[v]]))
                else:
                    d['maps'][ch] = self._sorted_gather_meta(inv[ch][first[v]], len(first[ch]))
            out[v] = d
        assert len(first[t.root]) == U
        out['_order'] = self._sorted_gather_meta(np.argsort(first[t.root]), U)
        return out

    def _cheb_dev(self, K):
        if K not in self._cheb:
            x, lam, T, D = cheb_tables(K)
            self._cheb[K] = tuple(self._T(a) for a in (x, lam, T, D))
        return self._cheb[K]

    def _gather(self, X, meta, dim):
        idx, perm, lengths, nrow = meta
        if self.deterministic:
            return SortedGather.apply(X, idx, dim, perm, lengths, nrow)
        return X.index_select(dim, idx)

    # ------------------------------------------------------------------ parameters
    def _unpack(self, x):
        a, b, c = torch.exp(x[0]), torch.exp(x[1]), torch.exp(x[2])
        pi = torch.stack([a * b, (1 - a) * c, a * (1 - b), (1 - a) * (1 - c)])
        kappa, r2, r3 = torch.exp(x[3]), torch.exp(x[4]), torch.exp(x[5])
        p = torch.exp(x[7])
        tau = torch.exp(x[6]) / p
        t = self.tree
        xr = x[8:]
        blen = {}
        if t.tie_root:
            tot = torch.exp(xr[0])
            for e in t.edges:
                if e[0] == t.root:
                    blen[e] = tot * (0.1 if e[1] == t.dup else 0.9)
            rest = xr[1:]
        else:
            rest = xr
        for i, e in enumerate(t.free_edges):
            blen[e] = torch.exp(rest[i])
        t_post = torch.stack([blen[e] for e in t.post_edges])
        t_pre = torch.stack([blen[e] for e in t.pre_edges]) if t.pre_edges else None
        S = self.OFF + (kappa - 1.0) * self.TR
        Q0 = S * pi[None, :]
        Q0 = Q0 / torch.sum(pi * Q0.sum(1))
        Qs = {}
        for s_, r in ((1, None), (2, r2), (3, r3)):
            Qr = Q0 if r is None else Q0 * r
            Qs[s_] = Qr - torch.diag_embed(Qr.sum(1))
        return pi, Qs, tau, p, t_post, t_pre

    def plan(self, x, K):
        """Host-side plan from exact arithmetic on x (no device sync)."""
        x = np.asarray(x, float)
        par = unpack_x(x, self.tree)
        pl = Plan()
        p = par.p
        lw = np.log1p(-p) if p < 1 else -np.inf
        w = np.exp(self.n_present * lw) if p < 1 else np.zeros(len(self.n_present))
        near = w >= self.w_tol
        if near.any():
            pl.lo, pl.hi = float(w[near].min()), float(w[near].max())
            pl.ncut = int(self.n_present[near].max())
            pl.near = True
            pl.K = int(K) if pl.hi - pl.lo >= 1e-14 else 1
        else:
            pl.lo, pl.hi, pl.ncut, pl.near, pl.K = 0.0, 0.5, 0, False, 0
        Qs = {s: par.site_Q(s) for s in (1, 2, 3)}
        tpost = max(par.blen[e] for e in self.tree.post_edges)
        tpre = max([par.blen[e] for e in self.tree.pre_edges], default=0.0)
        lam_ps = lam_pre = 0.0
        for cl in self.classes:
            sa, sb = cl[0], cl[1]
            qf = np.concatenate([Qs[sa].ravel(), Qs[sb].ravel()])
            dA = qf @ self.dEpm + par.tau * self.dEA
            dB = par.tau * self.dEB
            lam_ps = max(lam_ps, np.abs(dA + pl.lo * dB).max(), np.abs(dA + pl.hi * dB).max(), np.abs(dA).max())
            lam_pre = max(lam_pre, np.abs(qf @ self.dEpre).max())
        lam_is = max(np.abs(Qs[s['s']].ravel() @ self.dEis + par.tau * self.dEAis).max() for s in self.sites)
        lam_4 = max(np.abs(np.diag(Qs[s['s']])).max() for s in self.sites)
        pl.mu_ps = lam_ps * tpost
        pl.ps = unif_plan(lam_ps * tpost * (1 + 1e-12))
        pl.pre = unif_plan(lam_pre * tpre * (1 + 1e-12))
        pl.is16 = unif_plan(lam_is * tpost * (1 + 1e-12))
        pl.is4 = unif_plan(lam_4 * tpre * (1 + 1e-12))
        return pl

    def _buffer(self, x, pl):
        return np.concatenate([np.asarray(x, float), [pl.lo, pl.hi, float(pl.ncut)]])

    # ------------------------------------------------------------------ pieces
    def _prune(self, rep, Ppost, Ppre, prior, diag):
        """Pruning on repeat-compressed patterns. Ppost / Ppre: lists over edges of (..., n, n)."""
        t = self.tree
        part = {}
        for v in t.postorder:
            kids = t.children[v]
            if not kids:
                continue
            rv = rep[v]
            prod = None
            for ch in kids:
                e = (v, ch)
                if t.is_post[v]:
                    Pe = Ppost[self.ie_post[e]]
                    if v == t.dup:
                        Pe = Pe.index_select(-2, diag)
                else:
                    Pe = Ppre[self.ie_pre[e]]
                if not t.children[ch]:
                    tip = rv['tips'][ch]
                    u = ColGather.apply(Pe, tip[1], tip[2]) if tip[0] == 'idx' else torch.matmul(Pe, tip[1])
                else:
                    u = self._gather(torch.matmul(Pe, part.pop(ch)), rv['maps'][ch], -1)
                prod = u if prod is None else prod * u
            part[v] = prod
        Lr = torch.einsum('s,...su->...u', prior, part[t.root])
        return self._gather(Lr, rep['_order'], -1)

    def _site_lncol(self, pi, Qs, tau, t_post, t_pre, pl):
        """Single-site IS lnL per alignment column (L,)."""
        Qf = torch.stack([Qs[s['s']].reshape(-1) for s in self.sites])
        Pis = expm_unif((Qf @ self.Eis + tau * self.EAis).reshape(-1, 16, 16), t_post, *pl.is16)
        P4 = expm_unif(torch.stack([Qs[s['s']] for s in self.sites]), t_pre, *pl.is4) if t_pre is not None else None
        parts = []
        for si, s in enumerate(self.sites):
            Lc = self._prune(s['rep'], Pis[si].unbind(0), [] if P4 is None else P4[si].unbind(0), pi, self.diag16)
            parts.append(torch.log(Lc).index_select(0, s['inv']))
        return torch.cat(parts).index_select(0, self.col_perm)

    def _class_mats(self, Qs, tau, wk, t_post, t_pre, pl, cis):
        """Post-duplication matrices (C, K, E, 256, 256) and pre-duplication (C, E_pre, 16, 16)."""
        Qflat = torch.stack([torch.cat([Qs[self.classes[ci][0]].reshape(-1), Qs[self.classes[ci][1]].reshape(-1)])
                             for ci in cis])
        A = Qflat @ self.Epm + tau * self.EA
        B = tau * self.EB
        Kn = wk.shape[0]
        Gk = (A[:, None, :] + wk[None, :, None] * B[None, None, :]).reshape(-1, 256, 256)
        P = expm_unif(Gk, t_post, *pl.ps).reshape(len(cis), Kn, t_post.shape[0], 256, 256)
        Ppre = expm_unif((Qflat @ self.Epre).reshape(-1, 16, 16), t_pre, *pl.pre) if t_pre is not None else None
        return P, Ppre

    def _basis(self, w, near, lo, hi, K):
        """Barycentric basis (N, K) at w (differentiable through BaryBasis); out-of-interval rows (far keys)
        are evaluated at lo so that nothing divides by an underflowed denominator."""
        xk, lamk, _, Dk = self._cheb_dev(K)
        wd = w.detach()
        ws = torch.where(near, wd, lo)
        xq = 2 * (ws - lo) / (hi - lo) - 1
        Dd = xq[:, None] - xk[None, :]
        hit = Dd.abs() < 1e-15
        Wt = torch.where(hit, torch.zeros_like(Dd), lamk[None, :] / torch.where(hit, torch.ones_like(Dd), Dd))
        Bm = Wt / Wt.sum(1, keepdim=True)
        Bm = torch.where(hit.any(1, keepdim=True), hit.to(DT), Bm)
        BD = Bm @ (Dk * (2.0 / (hi - lo)))
        wsafe = torch.where(near, w, w.detach() * 0 + lo)
        return BaryBasis.apply(wsafe, Bm, BD)

    # ------------------------------------------------------------------ one evaluation
    def _forward(self, xb, pl, cis=None, far=True, want_pairs=False, want_keys=False):
        """lnL over the far part (if far) and classes cis (default all), plus diagnostics, as device tensors
        with no host synchronisation (CUDA-graph safe)."""
        nx = self.nx
        x = xb[:nx]
        lo, hi, ncut = xb[nx].detach(), xb[nx + 1].detach(), xb[nx + 2].detach()
        pi, Qs, tau, p, t_post, t_pre = self._unpack(x)
        C = len(self.classes)
        cis = list(range(C)) if cis is None else list(cis)
        f = torch.zeros((), dtype=DT, device=self.dev)
        tails = torch.zeros(C, dtype=DT, device=self.dev)
        minL = torch.full((), float('inf'), dtype=DT, device=self.dev)
        nonpos = torch.zeros((), dtype=DT, device=self.dev)
        lncol = None
        if far or want_pairs:
            lncol = self._site_lncol(pi, Qs, tau, t_post, t_pre, pl)
            if far:
                fc = self.suf.index_select(1, (ncut + 1).long()).squeeze(1)
                f = f + (fc * lncol).sum()
        lk_all = {}
        if pl.near and cis:
            K = pl.K
            if K > 1:
                xk = self._cheb_dev(K)[0]
                wk = lo + (hi - lo) * (xk + 1) / 2
            else:
                wk = hi.reshape(1)
            P, Ppre = self._class_mats(Qs, tau, wk, t_post, t_pre, pl, cis)
            prior16 = torch.outer(pi, pi).reshape(-1)
            lw = torch.log1p(-p)
            P_l = P.unbind(0)
            Ppre_l = Ppre.unbind(0) if Ppre is not None else None
            for j, ci in enumerate(cis):
                b = self.blk[ci]
                Lk = self._prune(b['rep'], P_l[j].unbind(1), [] if Ppre_l is None else Ppre_l[j].unbind(0),
                                 prior16, self.diag256)                                  # (K, U)
                near = b['key_n'] <= ncut
                if K > 1:
                    Bk = self._basis(torch.exp(b['key_n'] * lw), near, lo, hi, K)        # (Nk, K)
                    Lg = self._gather(Lk.t(), b['key_gather'], 0)                         # (Nk, K)
                    Lkey = (Lg * Bk).sum(1)
                    with torch.no_grad():
                        Ld = Lk.detach()
                        used = torch.zeros(b['U'], dtype=DT, device=self.dev).scatter_reduce(
                            0, b['key_u'], near.to(DT), 'amax') > 0
                        T2 = self._cheb_dev(K)[2][-2:]
                        tl = (T2 @ Ld).abs().sum(0)
                        Lmin = Ld.amin(0)
                        Lpos = torch.where(Lmin > 0, Lmin, torch.ones_like(Lmin))
                        thr = torch.clamp(NOISE_C * EPS * Ld.amax(0) / Lpos, max=ERR_CAP).clamp_min(self.tail_tol)
                        tu = torch.where(used, tl / Lpos / thr, torch.zeros_like(tl))
                        bad = (used & ~(Lmin > 0)).any()
                        tails[ci] = torch.where(bad, torch.full_like(tu[0], float('inf')), tu.max())
                else:
                    Lkey = self._gather(Lk[0], b['key_gather'], 0)
                lk = torch.log(torch.where(near, Lkey, torch.ones_like(Lkey)))
                f = f + (b['cnt'] * lk).sum()
                with torch.no_grad():
                    Lkd = Lkey.detach()
                    minL = torch.minimum(minL, torch.where(near, Lkd, torch.full_like(Lkd, float('inf'))).min())
                    nonpos = nonpos + (near & ~(Lkd > 0)).sum().to(DT)
                lk_all[ci] = lk
        out = {'f': f, 'tails': tails, 'minL': minL, 'nonpos': nonpos}
        if want_pairs:
            pp = []
            for ci, b in enumerate(self.blk):
                farv = lncol.index_select(0, b['I']) + lncol.index_select(0, b['J'])
                if ci in lk_all:
                    nearp = b['pair_n'] <= ncut
                    pp.append(torch.where(nearp, lk_all[ci].index_select(0, b['pair_key']), farv))
                else:
                    pp.append(farv)
            out['pairs'] = torch.cat(pp)
        if want_keys:
            out['lk'] = lk_all
        return out

    def _pack(self, o, g):
        parts = [o['f'].reshape(1)]
        if g is not None:
            parts.append(g.reshape(-1))
        parts += [o['tails'], o['minL'].reshape(1), o['nonpos'].reshape(1)]
        return torch.cat(parts)

    def _fg_tensor(self, xb, pl, grad, cis, far):
        if grad:
            x = xb[:self.nx].detach().requires_grad_(True)
            xb2 = torch.cat([x, xb[self.nx:].detach()])
            o = self._forward(xb2, pl, cis, far)
            g, = torch.autograd.grad(o['f'], x, allow_unused=True)
            if g is None:
                g = torch.zeros_like(x)
            return self._pack(o, g)
        with torch.no_grad():
            return self._pack(self._forward(xb, pl, cis, far), None)

    # ------------------------------------------------------------------ graphs
    def _graph_run(self, buf, pl, grad, cis, far):
        """Replay a cached graph of the same family whose series lengths m all cover the plan's (extra terms
        only add accuracy); otherwise capture one with m = elementwise max(plan, largest cached), so that the
        number of captures per family stays small."""
        key = pl.key(grad, (tuple(cis) if cis is not None else None, far))
        fam = self._graphs.setdefault(key, {})
        need = pl.ms()
        cover = [m for m in fam if all(a >= b for a, b in zip(m, need))]
        if cover:
            g = fam[min(cover, key=sum)]
        else:
            ms = tuple(max([n] + [m[i] for m in fam]) for i, n in enumerate(need))
            plc = pl.with_ms(ms)
            if self.graph_mem_gb and torch.cuda.memory_reserved(self.dev) > self.graph_mem_gb * 1e9 and self._graphs:
                self.clear_graphs()                          # memory budget: start the cache afresh
                self.n_graph_clears += 1
                fam = self._graphs.setdefault(key, {})
            sb = torch.as_tensor(buf, dtype=DT, device=self.dev).clone()
            if self._cap_stream is None:                     # one warm-up stream: the caching allocator keeps
                self._cap_stream = torch.cuda.Stream(self.dev)   # blocks per stream, a new stream per capture leaks
            st = self._cap_stream
            st.wait_stream(torch.cuda.current_stream(self.dev))
            with torch.cuda.stream(st):
                for _ in range(2):
                    self._fg_tensor(sb, plc, grad, cis, far)
            torch.cuda.current_stream(self.dev).wait_stream(st)
            graph = torch.cuda.CUDAGraph()
            with torch.cuda.graph(graph, pool=self._pool):
                out = self._fg_tensor(sb, plc, grad, cis, far)
            if self._pool is None:
                self._pool = graph.pool()
            g = fam[ms] = (graph, sb, out)
            self.n_capture += 1
            torch.cuda.empty_cache()                         # release warm-up blocks outside the graph pool
        graph, sb, out = g
        sb.copy_(torch.as_tensor(buf, dtype=DT))
        graph.replay()
        return out.detach().cpu().numpy().copy()

    def clear_graphs(self):
        self._graphs.clear()
        self._pool = None
        if self.dev.type == 'cuda':
            torch.cuda.empty_cache()

    def _eval_raw(self, x, pl, grad, cis=None, far=True):
        buf = self._buffer(x, pl)
        if self.graphs:
            try:
                return self._graph_run(buf, pl, grad, cis, far)
            except RuntimeError as e:                  # capture failure (driver/torch change): fall back to eager
                import warnings
                warnings.warn(f'psigc.gpu: CUDA-graph capture failed ({str(e)[:120]}); continuing eagerly')
                self.graphs = False
                self.capture_error = str(e)
                self.clear_graphs()
        xb = torch.as_tensor(buf, dtype=DT, device=self.dev)
        return self._fg_tensor(xb, pl, grad, cis, far).detach().cpu().numpy()

    def _unpack_out(self, v, grad):
        C = len(self.classes)
        f = float(v[0])
        o = 1
        g = None
        if grad:
            g = v[o:o + self.nx].copy(); o += self.nx
        tails = v[o:o + C]; o += C
        return f, g, tails, float(v[o]), int(round(v[o + 1]))

    def _use_parts(self, K):
        if isinstance(self.parts_mode, bool):
            return self.parts_mode
        return K > 32

    def _eval(self, x, pl, grad):
        """One evaluation at plan pl: one graph, or one per class (+ far) when parts are on."""
        if self._use_parts(pl.K) and pl.near:
            C = len(self.classes)
            f, g, tails, minL, nonpos = 0.0, (np.zeros(self.nx) if grad else None), np.zeros(C), np.inf, 0
            for part in [((), True)] + [((ci,), False) for ci in range(C)]:
                fi, gi, ti, mi, ni = self._unpack_out(self._eval_raw(x, pl, grad, part[0], part[1]), grad)
                f += fi
                if grad:
                    g += gi
                tails = np.maximum(tails, ti)
                minL, nonpos = min(minL, mi), nonpos + ni
            return f, g, tails, minL, nonpos
        return self._unpack_out(self._eval_raw(x, pl, grad), grad)

    # ------------------------------------------------------------------ exact evaluation of classes
    def _exact_class(self, x, pl, ci, grad, want_keys=False):
        """lnL of class ci with one node per distinct near distance (no interpolation), chunked; with grad the
        gradient is accumulated chunk by chunk (bounded memory). Returns (f, g or None[, per-key lnL])."""
        b = self.blk[ci]
        kn = b['key_n_np']
        nu = np.unique(kn[kn <= pl.ncut])
        f_tot, g_tot = 0.0, (np.zeros(self.nx) if grad else None)
        lk_keys = np.full(len(kn), np.nan) if want_keys else None
        for c0 in range(0, len(nu), self.exact_chunk):
            nc = nu[c0:c0 + self.exact_chunk]
            sel = np.flatnonzero(np.isin(kn, nc))
            pos = np.full(int(nc.max()) + 1, -1)
            pos[nc] = np.arange(len(nc))
            xt = torch.as_tensor(np.asarray(x, float), dtype=DT, device=self.dev)
            ctx = torch.enable_grad() if grad else torch.no_grad()
            with ctx:
                if grad:
                    xt.requires_grad_(True)
                pi, Qs, tau, p, t_post, t_pre = self._unpack(xt)
                wk = torch.exp(self._T(nc) * torch.log1p(-p))
                P, Ppre = self._class_mats(Qs, tau, wk, t_post, t_pre, pl, [ci])
                Lk = self._prune(b['rep'], P[0].unbind(1), [] if Ppre is None else Ppre[0].unbind(0),
                                 torch.outer(pi, pi).reshape(-1), self.diag256)       # (len(nc), U)
                Lsel = Lk[self._L(pos[kn[sel]]), b['key_u'][self._L(sel)]]
                lk = torch.log(Lsel)
                fc = (b['cnt'][self._L(sel)] * lk).sum()
                if grad:
                    gc, = torch.autograd.grad(fc, xt)
                    g_tot += gc.detach().cpu().numpy()
            f_tot += float(fc.detach())
            if want_keys:
                lk_keys[sel] = lk.detach().cpu().numpy()
        return (f_tot, g_tot, lk_keys) if want_keys else (f_tot, g_tot)

    def _eval_mixed(self, x, pl, grad, exact_cis):
        """Interpolated classes and the far part at pl, exact for exact_cis (eager)."""
        rest = [ci for ci in range(len(self.classes)) if ci not in exact_cis]
        xb = torch.as_tensor(self._buffer(x, pl), dtype=DT, device=self.dev)
        if self._use_parts(pl.K):                            # bounded memory: far part, then class by class
            C = len(self.classes)
            f, g, tails, minL, nonpos = 0.0, (np.zeros(self.nx) if grad else None), np.zeros(C), np.inf, 0
            for cis, far in [([], True)] + [([ci], False) for ci in rest]:
                fi, gi, ti, mi, ni = self._unpack_out(self._fg_tensor(xb, pl, grad, cis, far).detach().cpu().numpy(), grad)
                f += fi
                if grad:
                    g += gi
                tails, minL, nonpos = np.maximum(tails, ti), min(minL, mi), nonpos + ni
        else:
            f, g, tails, minL, nonpos = self._unpack_out(
                self._fg_tensor(xb, pl, grad, rest, True).detach().cpu().numpy(), grad)
        for ci in exact_cis:
            fe, ge = self._exact_class(x, pl, ci, grad)
            f += fe
            if grad:
                g += ge
            tails[ci] = 0.0
        return f, g, tails, minL, nonpos

    # ------------------------------------------------------------------ K policy
    def begin_task(self):
        """Start of an M2 task (or any sequence that may use a lower K again)."""
        if self.auto:
            self.K_floor = self.K_ladder[0]

    def _run(self, x, grad, kind):
        t0 = time.perf_counter()
        K = self.K_floor if self.auto else self.K_fixed
        tried, fallback = [], []
        while True:
            pl = self.plan(x, K)
            f, g, tails, minL, nonpos = self._eval(x, pl, grad)
            ok_cls = (tails <= 1.0)                      # tails: max over used patterns of tail / threshold
            ok = bool(ok_cls.all()) and nonpos == 0 and np.isfinite(f)
            tried.append((K, float(np.max(tails)) if len(tails) else 0.0, nonpos))
            if ok or not self.auto or not pl.near or pl.K <= 1:
                break
            i = self.K_ladder.index(K) if K in self.K_ladder else len(self.K_ladder) - 1
            if i + 1 < len(self.K_ladder):
                K = self.K_ladder[i + 1]
                self.K_floor = max(self.K_floor, K)            # ratchet: K only rises within a task
                continue
            if self.exact_fallback:
                fallback = [ci for ci in range(len(self.classes)) if not ok_cls[ci]] or list(range(len(self.classes)))
                f, g, tails, minL, nonpos = self._eval_mixed(x, pl, grad, fallback)
            break
        rec = {'kind': kind, 'K': pl.K, 'tried': tried, 'tail_max': float(np.max(tails)) if len(tails) else 0.0,
               'tails': {str(c): float(t) for c, t in zip(self.classes, tails)},
               'tail_ok': bool(np.all(tails <= 1.0)), 'fallback': [str(self.classes[c]) for c in fallback],
               'minL': minL, 'nonpos': nonpos, 'plan': (pl.ps, pl.pre, pl.is16, pl.is4), 'mu_ps': pl.mu_ps,
               'ncut': pl.ncut, 'interval': (pl.lo, pl.hi), 'time': time.perf_counter() - t0}
        self.last = rec
        if self.log_evals:
            xa = np.asarray(x, float)
            self.eval_log.append({**{k: rec[k] for k in ('kind', 'K', 'tail_max', 'fallback', 'nonpos', 'time')},
                                  'n_try': len(tried), 'psi': float(-xa[7]), 'log_tau': float(xa[6] - xa[7]),
                                  'K_need': self._K_need(tried)})
        return f, g

    def node_values(self, x, K):
        """Diagnostics: per class, the node values Lk (K, U) at x and the mask of used patterns (eager)."""
        pl = self.plan(x, K)
        out = {}
        with torch.no_grad():
            xb = torch.as_tensor(self._buffer(x, pl), dtype=DT, device=self.dev)
            x_, lo, hi, ncut = xb[:self.nx], xb[self.nx], xb[self.nx + 1], xb[self.nx + 2]
            pi, Qs, tau, p, t_post, t_pre = self._unpack(x_)
            wk = lo + (hi - lo) * (self._cheb_dev(pl.K)[0] + 1) / 2
            P, Ppre = self._class_mats(Qs, tau, wk, t_post, t_pre, pl, list(range(len(self.classes))))
            for ci, b in enumerate(self.blk):
                Lk = self._prune(b['rep'], P[ci].unbind(1), [] if Ppre is None else Ppre[ci].unbind(0),
                                 torch.outer(pi, pi).reshape(-1), self.diag256)
                used = np.zeros(b['U'], bool)
                used[b['key_u'].cpu().numpy()[b['key_n_np'] <= pl.ncut]] = True
                out[self.classes[ci]] = (Lk.cpu().numpy(), used)
        return out

    def _K_need(self, tried):
        """Smallest K tried at this point whose tails passed (None if none did): K escalation data (E7)."""
        ok = [k for k, t, n in tried if t <= 1.0 and n == 0]
        return min(ok) if ok else None

    # ------------------------------------------------------------------ public API
    def lnl(self, x, pairs=False):
        if pairs:
            f, _ = self._run(x, False, 'lnl')
            return f, self._pairs_at(x, self.last['K'])
        return self._run(x, False, 'lnl')[0]

    def lnl_grad(self, x):
        return self._run(x, True, 'lnl_grad')

    def _pairs_at(self, x, K):
        """Per-pair lnL at K (eager); classes that needed the exact fallback are filled exactly."""
        pl = self.plan(x, max(K, 1) if K else self.K_floor)
        with torch.no_grad():
            xb = torch.as_tensor(self._buffer(x, pl), dtype=DT, device=self.dev)
            o = self._forward(xb, pl, want_pairs=True)
            pp = o['pairs'].cpu().numpy()
        fb = self.last.get('fallback', [])
        for ci, b in enumerate(self.blk):
            if str(self.classes[ci]) in fb:
                _, _, lk = self._exact_class(x, pl, ci, False, want_keys=True)
                seg = slice(b['off'], b['off'] + b['npair'])
                nearp = b['pair_n'].cpu().numpy() <= pl.ncut
                kv = lk[b['pair_key'].cpu().numpy()]
                pp[seg] = np.where(nearp, kv, pp[seg])
        self.last['n_nonfinite_pairs'] = int(np.sum(~np.isfinite(pp)))
        return pp

    def pair_lnl_direct(self, x):
        """Exact per-pair lnL (every distinct near distance its own node); same order as lnl(pairs=True)."""
        pl = self.plan(x, self.K_ladder[0])
        with torch.no_grad():
            xb = torch.as_tensor(self._buffer(x, pl), dtype=DT, device=self.dev)
            o = self._forward(xb, pl, cis=[], far=False, want_pairs=True)
            pp = o['pairs'].cpu().numpy()
        if pl.near:
            for ci, b in enumerate(self.blk):
                _, _, lk = self._exact_class(x, pl, ci, False, want_keys=True)
                seg = slice(b['off'], b['off'] + b['npair'])
                nearp = b['pair_n'].cpu().numpy() <= pl.ncut
                pp[seg] = np.where(nearp, lk[b['pair_key'].cpu().numpy()], pp[seg])
        return pp

    def _scores_plan_checked(self, x):
        """The exact fallback has no forward mode. When the tail test sent classes to it, check the interpolant
        directly instead: per-pair lnL at K against the exact per-pair lnL of those classes, at the ladder's top K
        and then SCORE_K_EXTRA; use the first K whose largest error is <= SCORE_PAIR_TOL. (The tail test can fail
        on rounding alone: armL_UQ2 tract 100 sim 16, class (1, 3), tail 7-18 x threshold at K = 24-96 with
        per-pair errors <= 1.3e-14; gpu2/check_pair_err.py.)"""
        fb = set(self.last['fallback'])
        top = self.last['K'] or self.K_ladder[-1]
        exact, Ks = {}, [top] + [K for K in SCORE_K_EXTRA if K > top]
        tried = []
        for K in Ks:
            pl = self.plan(x, K)
            with torch.no_grad():
                xb = torch.as_tensor(self._buffer(x, pl), dtype=DT, device=self.dev)
                pp = self._forward(xb, pl, want_pairs=True)['pairs'].cpu().numpy()
            err = 0.0
            for ci, b in enumerate(self.blk):
                if str(self.classes[ci]) not in fb:
                    continue
                if ci not in exact:
                    exact[ci] = self._exact_class(x, pl, ci, False, want_keys=True)[2]
                seg = slice(b['off'], b['off'] + b['npair'])
                nearp = b['pair_n'].cpu().numpy() <= pl.ncut
                kv = exact[ci][b['pair_key'].cpu().numpy()]
                d = np.abs(pp[seg] - kv)[nearp]
                err = max(err, float(d.max()) if d.size else 0.0)
            tried.append((K, err))
            if np.isfinite(err) and err <= SCORE_PAIR_TOL:
                self.last['scores_check'] = {'K': K, 'tried': tried}
                return pl
        self.last['scores_check'] = {'K': None, 'tried': tried}
        raise RuntimeError(f'pair_scores: exact fallback needed and per-pair errors {tried} exceed {SCORE_PAIR_TOL}')

    def pair_scores(self, x, chunk=5):
        """Per-pair lnL (n_pairs,) and scores d lnL_ab / dx (n_pairs, nx) by forward mode, `chunk` directions
        per pass (torch.func.vmap over jvp). Evaluated at the current K after an lnl call has settled it."""
        self._run(x, False, 'pair_scores')
        pl = self.plan(x, self.last['K'] or self.K_floor)
        if self.last['fallback']:
            pl = self._scores_plan_checked(x)
        t0 = time.perf_counter()
        extra = torch.as_tensor(self._buffer(x, pl)[self.nx:], dtype=DT, device=self.dev)
        xt = torch.as_tensor(np.asarray(x, float), dtype=DT, device=self.dev)

        def fpp(xv):
            return self._forward(torch.cat([xv, extra]), pl, want_pairs=True)['pairs']

        prim, S = _jvp_columns(fpp, xt, self.nx, chunk, self.dev)
        if self.dev.type == 'cuda':
            torch.cuda.synchronize(self.dev)
        self.last['timing_scores'] = time.perf_counter() - t0
        return prim.detach().cpu().numpy(), S.detach().cpu().numpy()


    def site_scores(self, x, chunk=5):
        """Per-column single-site IS lnL (L,) and scores d lnL_a / dx (L, nx) by forward mode. The IS model
        depends on p only through tau, so in theta coordinates these scores have no psi component."""
        pl = self.plan(x, self.K_floor)
        xt = torch.as_tensor(np.asarray(x, float), dtype=DT, device=self.dev)

        def fs(xv):
            pi, Qs, tau, _, t_post, t_pre = self._unpack(xv)
            return self._site_lncol(pi, Qs, tau, t_post, t_pre, pl)

        prim, S = _jvp_columns(fs, xt, self.nx, chunk, self.dev)
        return prim.detach().cpu().numpy(), S.detach().cpu().numpy()

    def pair_sites(self):
        """Site indices (I, J) of every pair, in the order of lnl(pairs=True) and pair_scores."""
        I = np.concatenate([b['I'].cpu().numpy() for b in self.blk])
        J = np.concatenate([b['J'].cpu().numpy() for b in self.blk])
        return I, J


def _jvp_columns(fn, x, nx, chunk, dev):
    """(primal, Jacobian (m, nx)) by forward mode, `chunk` directions per pass; on CUDA out-of-memory the pass
    is retried with half as many directions (down to 1)."""
    eye = torch.eye(nx, dtype=DT, device=dev)
    cols, prim, c0 = [], None, 0
    while c0 < nx:
        oom = False
        try:
            pr, tg = _vmap_jvp(fn, x, eye[c0:c0 + chunk])
        except torch.OutOfMemoryError:
            if chunk == 1:
                raise
            oom = True
        if oom:                                   # outside the except block, so the traceback's tensors are freed
            pr = tg = None
            torch.cuda.empty_cache()
            chunk = max(1, chunk // 2)
            continue
        prim = pr if prim is None else prim
        cols.append(tg)
        c0 += tg.shape[0]
    return prim, torch.cat(cols, 0).t()


def _vmap_jvp(fn, x, tans):
    """(primal, tangents (n_dir, ...)) for several directions at once."""
    def one(tv):
        return torch.func.jvp(fn, (x,), (tv,))
    pr, tg = torch.func.vmap(one, randomness='error')(tans)
    return pr[0], tg
