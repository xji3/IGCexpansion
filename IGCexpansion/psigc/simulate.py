"""Simulator for HKY + PS-IGC on a tree with one duplication node (port of IGCexpansion.Simulator, HKY case).

Model (identical to the legacy Simulator, which generated the archived replicates):
* The simulated span is positions 1 .. max(nt_index) of seq_index; columns not listed in seq_index (e.g.
  238-240 of YDR418W_YEL054C) are simulated but not output. Their rate multiplier is 1 when
  legacy_missing_rate=True (legacy behaviour), otherwise their codon position continues the reading frame.
* Root sequence: i.i.d. draws from pi. Before the duplication node there is one copy; at the duplication node
  it is copied into two paralogs.
* Point mutations: site i of copy k mutates at rate -Q[s, s] * m_i (Q = normalised HKY, m_i = r of the codon
  position, r_1 = 1); the new state is drawn from Q[s, .] / -Q[s, s].
* IGC (post-duplication edges only): each ordered paralog pair (j -> k) initiates tracts at total rate
  eta * (L - 1 + 1/p); the start position is 0 with weight 1/p and 1..L-1 with weight 1 each; the tract
  length is Geometric(p) on {1, 2, ...}, truncated at L - 1; positions start..stop of copy k are overwritten
  by copy j. Per site and direction this gives the copying rate tau = eta / p of the likelihood.
* Exact Gillespie simulation; the total point-mutation rate is maintained incrementally.

Randomness: a numpy Generator; use rng_for(...) to give every task its own PCG64 stream derived from a
SeedSequence (study entropy, stage, cell, replicate), as SIMULATION_REDO_PLAN.md section 4 requires.
"""
from dataclasses import dataclass, field

import numpy as np

from .model import hky_Q


def rng_for(entropy, *key):
    """Independent PCG64 stream for one task: SeedSequence(entropy, spawn_key=key). key: small non-negative
    integers, e.g. (stage, cell_index, rep)."""
    ss = np.random.SeedSequence(entropy=entropy, spawn_key=tuple(int(k) for k in key))
    return np.random.Generator(np.random.PCG64(ss)), ss


@dataclass
class SimResult:
    tips: dict                    # species -> (2, L) post-duplication tips or (1, L) outgroup (states 0-3)
    igc_log: list                 # (edge, time, from_copy, to_copy, tract_length, start, stop, num_diff)
    n_pm: int                     # number of point mutations
    L: int
    columns: np.ndarray           # simulated positions (0-based) that form the output alignment
    extra: dict = field(default_factory=dict)

    def summary(self):
        lens = [e[6] - e[5] + 1 for e in self.igc_log]
        return {'N': len(self.igc_log), 'N_inf': int(sum(e[7] > 0 for e in self.igc_log)),
                'm_real': float(np.mean(lens)) if lens else None, 'n_pm': self.n_pm,
                'n_start0': int(sum(e[5] == 0 for e in self.igc_log)),
                'n_dir01': int(sum(e[2] == 0 for e in self.igc_log))}


def _site_multipliers(seq_index, r, legacy_missing_rate):
    nt = np.asarray(seq_index)[:, 0].astype(int)
    cp = np.asarray(seq_index)[:, 2].astype(int)
    L = int(nt.max())
    pos = np.full(L, -1, dtype=int)
    pos[nt - 1] = cp
    if not legacy_missing_rate:
        # continue the reading frame through unlisted positions
        last_nt, last_cp = None, None
        for i in range(L):
            if pos[i] >= 1:
                last_nt, last_cp = i, pos[i]
            elif last_cp is not None:
                pos[i] = (last_cp + (i - last_nt) - 1) % 3 + 1
    m = np.ones(L)
    for s in (2, 3):
        m[pos == s] = r[s]
    return m, nt - 1, L


def simulate(par, tree, seq_index, rng, legacy_missing_rate=True):
    """One replicate under psigc.Params `par` on psigc.Tree `tree`. Returns SimResult."""
    Q = np.asarray(hky_Q(par.pi, par.kappa, 1.0), float)
    qrow = -np.diag(Q).copy()
    jump = Q.copy()
    np.fill_diagonal(jump, 0.0)
    jump /= qrow[:, None]
    jump_cum = np.cumsum(jump, 1)
    m, cols, L = _site_multipliers(seq_index, par.r, legacy_missing_rate)
    p, eta = float(par.p), float(par.eta)
    W = L - 1 + 1.0 / p                                   # start-position weight total
    igc_rate_dir = eta * W                                # per ordered pair
    seqs = {tree.root: rng.choice(4, size=(1, L), p=np.asarray(par.pi, float))}
    igc_log, n_pm = [], 0

    def evolve(S, t, edge, igc):
        nonlocal n_pm
        S = S.copy()
        k = S.shape[0]
        R = qrow[S] * m[None, :]
        tot_pm = R.sum()
        tot_igc = 2 * igc_rate_dir if (igc and k == 2) else 0.0
        clock, since = 0.0, 0
        while True:
            total = tot_pm + tot_igc
            clock += rng.exponential(1.0 / total)
            if clock > t:
                break
            u = rng.random() * total
            if u < tot_pm:
                flat = np.cumsum(R.ravel())
                j = int(np.searchsorted(flat, u, side='right'))
                j = min(j, flat.size - 1)
                c, i = divmod(j, L)
                old = S[c, i]
                new = int(np.searchsorted(jump_cum[old], rng.random(), side='right'))
                new = min(new, 3)
                S[c, i] = new
                nr = qrow[new] * m[i]
                tot_pm += nr - R[c, i]
                R[c, i] = nr
                n_pm += 1
            else:
                fr = int(rng.integers(2))
                to = 1 - fr
                v = rng.random() * W
                start = 0 if v < 1.0 / p else min(L - 1, int(v - 1.0 / p) + 1)
                length = int(rng.geometric(p))
                stop = min(start + length - 1, L - 1)
                nd = int(np.sum(S[to, start:stop + 1] != S[fr, start:stop + 1]))
                S[to, start:stop + 1] = S[fr, start:stop + 1]
                R[to, start:stop + 1] = qrow[S[to, start:stop + 1]] * m[start:stop + 1]
                tot_pm = R.sum()
                igc_log.append(('_'.join(edge), clock, fr, to, length, start, stop, nd))
            since += 1
            if since % 256 == 0:
                tot_pm = R.sum()                          # guard against drift
        return S

    for v in tree.postorder[::-1]:                        # preorder: parents first
        for ch in tree.children[v]:
            e = (v, ch)
            S = evolve(seqs[v], par.blen[e], e, tree.is_post[v])
            if ch == tree.dup and not tree.is_post[v]:
                S = np.vstack([S, S])                     # duplication: two identical paralogs
            seqs[ch] = S
    tips = {sp: seqs[sp] for sp in tree.tips}
    return SimResult(tips=tips, igc_log=igc_log, n_pm=n_pm, L=L, columns=cols)


def to_dataset(res, tree, seq_index):
    """psigc.Dataset from a SimResult (columns ordered as Dataset expects: post tips (2 copies), then pre tips)."""
    from .data import Dataset
    si = np.asarray(seq_index)
    cols = [res.tips[sp][k, res.columns] for sp in tree.post_tips for k in (0, 1)]
    cols += [res.tips[sp][0, res.columns] for sp in tree.pre_tips]
    return Dataset(np.array(cols).T, si[:, 0], si[:, 2], tree)


def write_fasta(res, tree, paralogs, path, outgroup_gene=None):
    nt = 'ACGT'
    with open(path, 'w') as f:
        for sp in tree.post_tips:
            for k, g in enumerate(paralogs):
                f.write(f'>{sp}__{g}\n' + ''.join(nt[s] for s in res.tips[sp][k, res.columns]) + '\n')
        for sp in tree.pre_tips:
            f.write(f'>{sp}__{outgroup_gene or paralogs[0]}\n' + ''.join(nt[s] for s in res.tips[sp][0, res.columns]) + '\n')
