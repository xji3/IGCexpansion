"""HKY + PS-IGC model: parameters, legacy x layout, and generator specs.

A generator spec lists the off-diagonal entries of a CTMC generator as two kinds of terms:

* point-mutation (PM) terms: Q[src, dst] += Qsite[slot][cur, y], where Qsite[slot] is the 4 x 4 HKY
  matrix of the site that `slot` stands for;
* IGC terms: Q[src, dst] += tau * (cA + w * cB).

Diagonals are minus the row sums. The same lists give the adjoint fold: for a matrix M = d lnL / dQ
(full, including the diagonal), the derivative with respect to an off-diagonal entry is
M[src, dst] - M[src, src], and summing it over the term lists gives d lnL / d Qsite[slot][cur, y] and
d lnL / d tau.

State conventions (as in the manuscript and the legacy PSJSModel):
* PS256: (i_a, i_b, j_a, j_b), index 64 i_a + 16 i_b + 4 j_a + j_b; i, j are the two paralogs,
  a, b the two sites. Slot 0 = site a, slot 1 = site b.
* IS16:  (i, j) at one site, index 4 i + j. Slot 0.
* PRE16: (a, b) single copy at two independent sites, index 4 a + b. Slots 0, 1.
* PRE4:  one nucleotide. Slot 0.
"""
from dataclasses import dataclass, field

import numpy as np

NT = 'ACGT'
TRANSITIONS = {(0, 2), (2, 0), (1, 3), (3, 1)}


def hky_Q(pi, kappa, r=1.0):
    """HKY generator normalised to one expected substitution per unit time, times r.

    Written without abs/max so that it works with complex inputs (complex-step derivatives)."""
    pi = np.asarray(pi)
    dtype = np.result_type(pi, kappa, r, float)
    S = np.ones((4, 4), dtype=dtype)
    for x, y in TRANSITIONS:
        S[x, y] = kappa
    np.fill_diagonal(S, 0.0)
    Q = S * pi[None, :]
    Q = Q / np.sum(pi * Q.sum(1))
    Q = Q * r
    Q = Q - np.diag(Q.sum(1))
    return Q


def freq_from_x(x3):
    """Legacy PMModel.unpack_frequency: x = log(%AG, %A among AG, %C among CT)."""
    a, b, c = np.exp(x3[0]), np.exp(x3[1]), np.exp(x3[2])
    return np.array([a * b, (1 - a) * c, a * (1 - b), (1 - a) * (1 - c)])


def x_from_freq(pi):
    pi = np.asarray(pi, float)
    ag = pi[0] + pi[2]
    return np.log([ag, pi[0] / ag, pi[1] / (1 - ag)])


@dataclass
class Params:
    pi: np.ndarray
    kappa: float
    r: dict                      # codon position -> rate multiplier, r[1] = 1
    eta: float                   # IGC initiation rate (legacy init_rate)
    p: float                     # tract continuation parameter; mean tract length 1/p
    blen: dict = field(default_factory=dict)   # edge (parent, child) -> length

    @property
    def tau(self):
        return self.eta / self.p

    def site_Q(self, s):
        return hky_Q(self.pi, self.kappa, self.r[s])


N_PM = 6        # x[0:6] = log(%AG, %A|AG, %C|CT), log kappa, log r2, log r3
I_ETA, I_P = 6, 7


def unpack_pm(xpm):
    """Complex-step friendly: returns (pi, kappa, r) from x[0:6]."""
    pi = freq_from_x(xpm[:3])
    return pi, np.exp(xpm[3]), {1: 1.0, 2: np.exp(xpm[4]), 3: np.exp(xpm[5])}


def unpack_x(x, tree):
    """Legacy PSJSGeneconv x layout -> Params. x = x_pm (6) + log(eta, p) + log edge rates."""
    x = np.asarray(x, float)
    pi, kappa, r = unpack_pm(x[:N_PM])
    blen = tree.blen_from_x(x[I_P + 1:])
    return Params(pi=pi, kappa=kappa, r=r, eta=float(np.exp(x[I_ETA])), p=float(np.exp(x[I_P])), blen=blen)


def pack_x(par, tree):
    xpm = np.concatenate([x_from_freq(par.pi), np.log([par.kappa, par.r[2], par.r[3]])])
    return np.concatenate([xpm, np.log([par.eta, par.p]), tree.x_from_blen(par.blen)])


# ---------------------------------------------------------------- generator specs
class Spec:
    """Sparse description of a generator family (see module docstring)."""

    def __init__(self, n, pm, igc, nslot):
        self.n, self.nslot = n, nslot
        pm = np.array(pm, dtype=np.int64).reshape(-1, 5)
        self.pm_src, self.pm_dst, self.pm_slot, self.pm_cur, self.pm_y = pm.T
        self.pm_flat_q = self.pm_cur * 4 + self.pm_y
        self.pm_flat_h = self.pm_slot * 16 + self.pm_cur * 4 + self.pm_y
        self.pm_flat_g = self.pm_src * n + self.pm_dst
        if igc:
            ig = np.array(igc, dtype=float).reshape(-1, 4)
            self.ig_src, self.ig_dst = ig[:, 0].astype(np.int64), ig[:, 1].astype(np.int64)
            self.ig_cA, self.ig_cB = ig[:, 2], ig[:, 3]
        else:
            self.ig_src = self.ig_dst = np.zeros(0, np.int64)
            self.ig_cA = self.ig_cB = np.zeros(0)
        self.ig_flat_g = self.ig_src * n + self.ig_dst

    def build(self, Qslots, tau=0.0, w=0.0):
        """Dense generator for site matrices Qslots (list of 4 x 4), tau, w."""
        n = self.n
        vals = np.empty(len(self.pm_src))
        for s in range(self.nslot):
            m = self.pm_slot == s
            vals[m] = np.asarray(Qslots[s]).ravel()[self.pm_flat_q[m]]
        G = np.bincount(self.pm_flat_g, weights=vals, minlength=n * n)
        if len(self.ig_src) and tau != 0.0:
            G = G + np.bincount(self.ig_flat_g, weights=tau * (self.ig_cA + w * self.ig_cB), minlength=n * n)
        G = G.reshape(n, n)
        G[np.diag_indices(n)] -= G.sum(1)
        return G

    def build_AB(self, Qslots, tau):
        """(A, B) with G(w) = A + w B."""
        A = self.build(Qslots, tau, 0.0)
        Bv = np.bincount(self.ig_flat_g, weights=tau * self.ig_cB, minlength=self.n * self.n).reshape(self.n, self.n)
        Bv[np.diag_indices(self.n)] -= Bv.sum(1)
        return A, Bv

    def fold(self, M):
        """M: d lnL / dG (n x n). Returns H (nslot, 4, 4) = d lnL / d Qsite off-diagonals and the
        IGC sums (sA, sB) = (sum Mt cA, sum Mt cB), so that d lnL / d tau = sA + w sB."""
        Mt = M - np.diag(M)[:, None]
        h = np.bincount(self.pm_flat_h, weights=Mt.ravel()[self.pm_flat_g], minlength=self.nslot * 16)
        mig = Mt.ravel()[self.ig_flat_g]
        return h.reshape(self.nslot, 4, 4), float(mig @ self.ig_cA), float(mig @ self.ig_cB)


def idx256(ia, ib, ja, jb):
    return 64 * ia + 16 * ib + 4 * ja + jb


def _spec_ps256():
    pm, igc = [], []
    other = {0: 2, 1: 3, 2: 0, 3: 1}          # same site, other paralog
    for s in range(256):
        st = [s // 64, (s // 16) % 4, (s // 4) % 4, s % 4]
        ia, ib, ja, jb = st
        for pos in range(4):
            cur = st[pos]
            site_a = pos in (0, 2)
            slot = 0 if site_a else 1
            other_site_identical = (ib == jb) if site_a else (ia == ja)
            for y in range(4):
                if y == cur:
                    continue
                t = list(st); t[pos] = y
                d = idx256(*t)
                pm.append((s, d, slot, cur, y))
                if y == st[other[pos]]:
                    # tau (1 - w) + tau w if the partner site is identical between paralogs, else tau (1 - w)
                    igc.append((s, d, 1.0, 0.0 if other_site_identical else -1.0))
        if ia != ja and ib != jb:
            igc.append((s, idx256(ja, jb, ja, jb), 0.0, 1.0))   # i copies j at both sites
            igc.append((s, idx256(ia, ib, ia, ib), 0.0, 1.0))   # j copies i at both sites
    return Spec(256, pm, igc, 2)


def _spec_is16():
    pm, igc = [], []
    for i in range(4):
        for j in range(4):
            s = 4 * i + j
            for y in range(4):
                if y != i:
                    pm.append((s, 4 * y + j, 0, i, y))
                    if y == j:
                        igc.append((s, 4 * y + j, 1.0, 0.0))
                if y != j:
                    pm.append((s, 4 * i + y, 0, j, y))
                    if y == i:
                        igc.append((s, 4 * i + y, 1.0, 0.0))
    return Spec(16, pm, igc, 1)


def _spec_pre16():
    pm = []
    for a in range(4):
        for b in range(4):
            s = 4 * a + b
            for y in range(4):
                if y != a:
                    pm.append((s, 4 * y + b, 0, a, y))
                if y != b:
                    pm.append((s, 4 * a + y, 1, b, y))
    return Spec(16, pm, [], 2)


def _spec_pre4():
    return Spec(4, [(x, y, 0, x, y) for x in range(4) for y in range(4) if x != y], [], 1)


PS256 = _spec_ps256()
IS16 = _spec_is16()
PRE16 = _spec_pre16()
PRE4 = _spec_pre4()

DIAG256 = np.array([idx256(a, b, a, b) for a in range(4) for b in range(4)])   # -> PRE16 index 4a + b
DIAG16 = np.array([4 * a + a for a in range(4)])                             # -> PRE4 index a


def site_Q_jacobian(xpm, s):
    """Complex-step Jacobian d Qsite_s / d x_pm, shape (6, 4, 4); exact to rounding."""
    h = 1e-30
    J = np.zeros((N_PM, 4, 4))
    for j in range(N_PM):
        xc = np.array(xpm, dtype=complex)
        xc[j] += 1j * h
        pi, kappa, r = unpack_pm(xc)
        J[j] = hky_Q(pi, kappa, r[s]).imag / h
    return J


def pi_jacobian(xpm):
    """d pi / d x_pm, shape (6, 4)."""
    h = 1e-30
    J = np.zeros((N_PM, 4))
    for j in range(N_PM):
        xc = np.array(xpm, dtype=complex)
        xc[j] += 1j * h
        J[j] = freq_from_x(xc[:3]).imag / h
    return J
