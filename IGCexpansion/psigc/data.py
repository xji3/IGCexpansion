"""Alignment -> site pairs, codon-position classes, unique (class, pattern) keys and tip states."""
import numpy as np
import scipy.sparse as sp

NT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
CLASSES = [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            name = line[1:].strip()
            seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: ''.join(v).upper() for k, v in seqs.items()}


def _find(seqs, species, gene):
    exact = f'{species}__{gene}'
    if exact in seqs:
        return seqs[exact]
    hits = [k for k in seqs if k.startswith(species) and gene in k]
    if len(hits) != 1:
        raise KeyError(f'cannot find a unique sequence for {species} / {gene}: {hits}')
    return seqs[hits[0]]


class TipData:
    """Tip observations for U patterns in an n-state space: one-hot indices when unambiguous,
    otherwise an indicator mask (n, U)."""

    def __init__(self, comps, n):
        """comps: list of (U,) int arrays, most significant first; -1 = missing."""
        comps = [np.asarray(c, np.int64) for c in comps]
        U = len(comps[0])
        self.n, self.U = n, U
        if all((c >= 0).all() for c in comps):
            idx = np.zeros(U, np.int64)
            for c in comps:
                idx = idx * 4 + c
            self.idx = idx
            self.mask = None
            self.onehot = sp.csr_matrix((np.ones(U), (np.arange(U), idx)), shape=(U, n))
        else:
            self.idx = None
            ind = []
            for c in comps:
                m = np.zeros((U, 4))
                ok = c >= 0
                m[ok, c[ok]] = 1.0
                m[~ok] = 1.0
                ind.append(m)
            mask = ind[0]
            for m in ind[1:]:
                mask = (mask[:, :, None] * m[:, None, :]).reshape(U, -1)
            self.mask = np.ascontiguousarray(mask.T)

    def up(self, P):
        """P @ tip partial, shape (P.shape[0], U)."""
        return P[:, self.idx] if self.mask is None else P @ self.mask

    def grad(self, bu):
        """d lnL / dP contribution for up = P @ tip, given bu = d lnL / d up."""
        return np.asarray(bu @ self.onehot) if self.mask is None else bu @ self.mask.T


class ClassBlock:
    """All near-or-far pairs of one codon-position class with their unique patterns."""

    def __init__(self, cls, inv, n, pair_I, pair_J, tips):
        self.cls, self.inv, self.n = cls, inv, n
        self.pair_I, self.pair_J = pair_I, pair_J
        self.tips = tips
        self.U = next(iter(tips.values())).U


class Dataset:
    def __init__(self, cols, nt_index, codon_pos, tree, names=None, locus=None):
        """cols: (L, 2 * n_post_tips + n_pre_tips) states 0-3, -1 missing; columns are (paralog 1,
        paralog 2) for each tree.post_tips entry, then one column per tree.pre_tips entry.
        codon_pos: 1-3, or 0 for non-coding sites (legacy seq_index convention); non-coding sites get rate 1
        (they are treated as position 1), so a fully non-coding dataset does not depend on r2, r3.
        locus: optional (L,) locus ids; pairs are formed within a locus only (multi-locus composite
        likelihood, SIMULATION_REDO_PLAN.md section 6), with distances from nt_index inside the locus."""
        self.cols = np.asarray(cols, np.int64)
        self.nt_index = np.asarray(nt_index, np.int64)
        cp = np.asarray(codon_pos, np.int64)
        assert ((cp >= 0) & (cp <= 3)).all()
        self.noncoding = bool((cp == 0).all())
        if (cp == 0).any() and not self.noncoding:
            raise ValueError('mixed coding / non-coding sites are not supported')
        self.codon_pos = np.where(cp == 0, 1, cp)
        self.locus = np.zeros(len(cp), np.int64) if locus is None else np.asarray(locus, np.int64)
        assert len(self.locus) == len(cp) and (np.diff(self.locus) >= 0).all(), 'loci must be contiguous'
        self.n_loci = int(len(np.unique(self.locus)))
        self.tree = tree
        self.L = self.cols.shape[0]
        assert self.cols.shape[1] == 2 * len(tree.post_tips) + len(tree.pre_tips)
        self._build_pairs()
        self._build_sites()

    @classmethod
    def from_fasta(cls, fasta, seq_index, tree, paralogs, outgroup_gene=None):
        seqs = read_fasta(fasta) if isinstance(fasta, str) else fasta
        si = np.loadtxt(seq_index, dtype=np.int64) if isinstance(seq_index, str) else np.asarray(seq_index)
        columns = []
        for spn in tree.post_tips:
            for g in paralogs:
                columns.append(_find(seqs, spn, g))
        for spn in tree.pre_tips:
            genes = [outgroup_gene] if outgroup_gene else list(paralogs)
            for g in genes:
                try:
                    columns.append(_find(seqs, spn, g)); break
                except KeyError:
                    continue
            else:
                raise KeyError(f'no sequence for outgroup {spn}')
        L = len(columns[0])
        assert all(len(c) == L for c in columns), 'unaligned sequences'
        assert si.shape[0] == L, f'seq_index has {si.shape[0]} rows, alignment {L} columns'
        cols = np.array([[NT.get(ch, -1) for ch in c] for c in columns], dtype=np.int64).T
        return cls(cols, si[:, 0], si[:, 2], tree)

    # ------------------------------------------------------------------ pairs
    def _tips_ps(self, A, B):
        tips = {}
        for k, spn in enumerate(self.tree.post_tips):
            tips[spn] = TipData([A[:, 2 * k], B[:, 2 * k], A[:, 2 * k + 1], B[:, 2 * k + 1]], 256)
        off = 2 * len(self.tree.post_tips)
        for k, spn in enumerate(self.tree.pre_tips):
            tips[spn] = TipData([A[:, off + k], B[:, off + k]], 16)
        return tips

    @classmethod
    def concat(cls, datasets):
        """Multi-locus dataset from single-locus datasets on the same tree (locus ids 0..K-1)."""
        tree = datasets[0].tree
        return cls(np.concatenate([d.cols for d in datasets]), np.concatenate([d.nt_index for d in datasets]),
                   np.concatenate([np.where(np.full(d.L, d.noncoding), 0, d.codon_pos) for d in datasets]), tree,
                   locus=np.concatenate([np.full(d.L, k) for k, d in enumerate(datasets)]))

    def pair_count_per_site(self):
        """Number of pairs each site belongs to (L_locus - 1)."""
        sizes = np.bincount(self.locus)
        return sizes[self.locus] - 1

    def _build_pairs(self):
        starts = np.r_[0, np.flatnonzero(np.diff(self.locus)) + 1, self.L]
        Is, Js = [], []
        for a, b in zip(starts[:-1], starts[1:]):
            i, j = np.triu_indices(b - a, 1)
            Is.append(i + a); Js.append(j + a)
        I, J = np.concatenate(Is), np.concatenate(Js)
        n = self.nt_index[J] - self.nt_index[I]
        assert (n > 0).all()
        sa, sb = self.codon_pos[I], self.codon_pos[J]
        swap = sa > sb
        Ic, Jc = np.where(swap, J, I), np.where(swap, I, J)      # canonical: codon position of a <= b
        self.n_pairs = len(I)
        self.blocks = {}
        for c in CLASSES:
            m = (self.codon_pos[Ic] == c[0]) & (self.codon_pos[Jc] == c[1])
            if not m.any():
                continue
            a, b = Ic[m], Jc[m]
            rows = np.concatenate([self.cols[a], self.cols[b]], axis=1).astype(np.int8)
            uniq, first, inv = np.unique(rows, axis=0, return_index=True, return_inverse=True)
            inv = inv.ravel()
            A, B = self.cols[a[first]], self.cols[b[first]]
            blk = ClassBlock(c, inv, n[m], a, b, self._tips_ps(A, B))
            blk.U = len(first)
            blk.A, blk.B = A, B
            self.blocks[c] = blk
        self.n_max = int(n.max())
        self.pair_n_all = n

    def chunk_blocks(self, max_patterns):
        """Blocks split so that no block has more than max_patterns unique pair patterns (pattern chunking for
        large datasets, e.g. the primate introns: 570k patterns). Keys are (pos_a, pos_b, chunk); a block that
        fits keeps chunk 0. Pairs of one chunk are those whose pattern falls in its pattern range."""
        out = {}
        for c, b in self.blocks.items():
            nk = max(1, -(-b.U // max_patterns))
            for k in range(nk):
                if nk == 1:
                    out[(c[0], c[1], 0)] = b
                    continue
                lo, hi = k * max_patterns, min((k + 1) * max_patterns, b.U)
                m = (b.inv >= lo) & (b.inv < hi)
                blk = ClassBlock(c, b.inv[m] - lo, b.n[m], b.pair_I[m], b.pair_J[m], self._tips_ps(b.A[lo:hi], b.B[lo:hi]))
                blk.U, blk.A, blk.B = hi - lo, b.A[lo:hi], b.B[lo:hi]
                out[(c[0], c[1], k)] = blk
        return out

    def _build_sites(self):
        """Single-site patterns per codon position (for the IS model and far pairs)."""
        self.site_blocks = {}
        npost = len(self.tree.post_tips)
        for s in (1, 2, 3):
            colsel = np.flatnonzero(self.codon_pos == s)
            if len(colsel) == 0:
                continue
            C = self.cols[colsel]
            uniq, first, inv = np.unique(C.astype(np.int8), axis=0, return_index=True, return_inverse=True)
            Cu = C[first]
            tips = {}
            for k, spn in enumerate(self.tree.post_tips):
                tips[spn] = TipData([Cu[:, 2 * k], Cu[:, 2 * k + 1]], 16)
            for k, spn in enumerate(self.tree.pre_tips):
                tips[spn] = TipData([Cu[:, 2 * npost + k]], 4)
            self.site_blocks[s] = dict(cols=colsel, inv=inv.ravel(), tips=tips, U=len(first))

    def summary(self):
        return {'columns': self.L, 'pairs': self.n_pairs,
                'unique_patterns': {f'{c}': b.U for c, b in self.blocks.items()},
                'unique_patterns_total': int(sum(b.U for b in self.blocks.values()))}
