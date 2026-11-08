"""jsonctmctree scene -> BEAGLE action instance (plan §4).

``compile_scene`` does everything that depends only on the tree and the data (topology,
operation lists, buffer and key maps, compressed tip partials); ``BeagleScene`` owns one
instance for a compiled scene and re-sends only values (generators, root prior, edge
lengths) on each call.

Topology mapping onto BEAGLE's binary operations:
* Unary node: second child is the shared ONES buffer on the zero-length key Z (exact: e^{0}=I).
* Multifurcation: left fold through virtual nodes joined by key Z.
* Observed internal node (or observed root): its observations become a pseudo-tip child on key Z.
* Unobserved leaf: a tip of ones.
Buffers: [tips + pseudo-tips | internal + virtual post | ONES | pre-order for every node].
Keys: key(v) = post buffer of v for every non-root node; Z = the ONES buffer index (no node
uses it as an edge).
Scaling: every post op and every pre op writes its own scale buffer (SCALERS_LOG); post scales
are accumulated into CUM for the likelihood, pre scales cancel in the derivative ratios.
"""

import hashlib
from collections import defaultdict

import numpy as np

from ._lib import NONE, f64, i32
from .instance import ActionInstance


class UnsupportedScene(ValueError):
    """The scene uses something the BEAGLE path does not handle; the engine falls back."""


def generator_triplets(shape, proc):
    """Off-diagonal rates (duplicates summed, zeros dropped) plus the diagonal -row sums."""
    S = int(np.prod(shape))
    rs, cs = np.asarray(proc["row_states"]), np.asarray(proc["column_states"])
    v = np.asarray(proc["transition_rates"], dtype=float)
    if len(v) == 0:
        return i32([]), i32([]), f64([])
    if rs.ndim != 2 or rs.shape[1] != len(shape) or cs.shape != rs.shape or len(v) != len(rs):
        raise UnsupportedScene("process definition does not match the state space shape")
    r = np.ravel_multi_index(rs.T, shape)
    c = np.ravel_multi_index(cs.T, shape)
    if np.any(r == c):
        raise UnsupportedScene("diagonal entry in a process definition")
    if np.any(~np.isfinite(v)) or np.any(v < 0):
        raise UnsupportedScene("transition rates must be finite and >= 0")
    key = r.astype(np.int64) * S + c
    uk, inv = np.unique(key, return_inverse=True)
    vs = np.bincount(np.asarray(inv).reshape(-1), weights=v)
    r, c = uk // S, uk % S
    keep = vs != 0.0
    r, c, vs = r[keep], c[keep], vs[keep]
    rowsum = np.bincount(r, weights=vs, minlength=S)
    d = np.flatnonzero(rowsum)
    return i32(np.concatenate([r, d])), i32(np.concatenate([c, d])), f64(np.concatenate([vs, -rowsum[d]]))


class CompiledScene:
    pass


def compile_scene(scene, compress=True):
    cs = CompiledScene()
    shape = cs.shape = [int(s) for s in scene["state_space_shape"]]
    S = cs.S = int(np.prod(shape))
    # jsonctmctree uses the process definitions themselves; process_count may disagree.
    G = cs.G = len(scene["process_definitions"])
    tr = scene["tree"]
    rows, cols = [int(u) for u in tr["row_nodes"]], [int(v) for v in tr["column_nodes"]]
    procs = [int(p) for p in tr["edge_processes"]]
    E = cs.E = len(cols)
    if E == 0:
        raise UnsupportedScene("tree without edges")
    if not (len(rows) == E == len(procs) == len(tr["edge_rate_scaling_factors"])):
        raise UnsupportedScene("tree arrays differ in length")
    if G == 0 or any(p < 0 or p >= G for p in procs):
        raise UnsupportedScene("edge process out of range")
    cs.rows, cs.cols, cs.procs = rows, cols, procs
    children, parent = defaultdict(list), {}
    for u, v in zip(rows, cols):
        if v in parent:
            raise UnsupportedScene("node with two parents")
        children[u].append(v)
        parent[v] = u
    nodes = set(rows) | set(cols)
    roots = nodes - set(parent)
    if len(roots) != 1:
        raise UnsupportedScene("tree must have exactly one root")
    root = cs.root = roots.pop()
    if "node_count" in scene and any(n < 0 or n >= int(scene["node_count"]) for n in nodes):
        raise UnsupportedScene("node index out of range")

    # Observation layout (which nodes/axes are observed); the rows themselves are pattern data.
    od = scene["observed_data"]
    obs_nodes = [int(n) for n in od["nodes"]]
    obs_vars = [int(a) for a in od["variables"]]
    if len(obs_vars) != len(obs_nodes) or not obs_nodes:
        raise UnsupportedScene("observed_data arrays are inconsistent or empty")
    if any(n not in nodes for n in obs_nodes) or any(a < 0 or a >= len(shape) for a in obs_vars):
        raise UnsupportedScene("observation on an unknown node or axis")
    cs.obs_nodes, cs.obs_vars = obs_nodes, obs_vars
    observed = set(obs_nodes)

    # Binary operation tree: ('real', n), ('virt', k), ('pseudo', n).
    leaves = sorted(n for n in nodes if not children[n])
    internal = sorted(n for n in nodes if children[n])
    pseudo = sorted(n for n in internal if n in observed)
    tips = [("real", n) for n in leaves] + [("pseudo", n) for n in pseudo]
    bchildren = {}           # binary node -> [(child, is_edge_with_key_Z)]
    virt = []
    for n in internal:
        ch = [(("real", c), False) for c in children[n]]
        if n in observed:
            ch.append((("pseudo", n), True))
        while len(ch) > 2:
            vnode = ("virt", len(virt))
            virt.append(vnode)
            bchildren[vnode] = ch[:2]
            ch = [(vnode, True)] + ch[2:]
        bchildren[("real", n)] = ch
    posts = [("real", n) for n in internal] + virt
    order = tips + posts
    idx = {b: i for i, b in enumerate(order)}
    cs.T = len(tips)
    ONES = len(order)
    npost = ONES + 1
    # Pre-order buffers for every binary node except pseudo-tips (never an edge of interest).
    pre_nodes = [b for b in order if b[0] != "pseudo"]
    pre = {b: npost + i for i, b in enumerate(pre_nodes)}
    cs.nbuf = npost + len(pre_nodes)
    cs.ONES, cs.Z = ONES, ONES
    cs.keys = cs.nbuf
    cs.idx, cs.pre, cs.root_b = idx, pre, ("real", root)
    cs.tips = tips
    # key of a binary node as a child: its own post index for real edges, Z otherwise.
    edge_of = {v: e for e, v in enumerate(cols)}

    def key(b, z):
        return cs.Z if z else idx[b]

    # Post-order over binary nodes.
    post_order, stack = [], [(cs.root_b, False)]
    while stack:
        b, done = stack.pop()
        if b not in bchildren:
            continue
        if done:
            post_order.append(b)
            continue
        stack.append((b, True))
        for c, _ in bchildren[b]:
            stack.append((c, False))
    cs.scale_of = {b: i for i, b in enumerate(post_order)}
    cs.CUM = len(post_order)
    ops = []
    for b in post_order:
        ch = bchildren[b]
        c0, z0 = ch[0]
        c1, k1 = (idx[ch[1][0]], key(*ch[1])) if len(ch) == 2 else (ONES, cs.Z)
        ops.append((idx[b], cs.scale_of[b], NONE, idx[c0], key(c0, z0), c1, k1))
    cs.post_ops = ops
    # Pre-order: parent before child, skipping pseudo-tips.
    pre_ops, top_inputs = [], []
    for b in reversed(post_order):
        ch = bchildren[b]
        for i, (c, z) in enumerate(ch):
            if c[0] == "pseudo":
                continue
            if len(ch) == 2:
                s, sz = ch[1 - i]
                sib, sk = idx[s], key(s, sz)
            else:
                sib, sk = ONES, cs.Z
            # Each pre op rescales its output (own write-scale buffer). Without it deep trees
            # underflow in the pre-order pass; the factor cancels in every derivative ratio.
            pre_ops.append((pre[c], cs.CUM + 1 + len(pre_ops), NONE, pre[b], key(c, z), sib, sk))
            top_inputs.append((c, sib, sk, b))
    cs.pre_ops = pre_ops
    cs.nscale = cs.CUM + 1 + len(pre_ops)
    cs.root_pre = pre[cs.root_b]
    # Real edges in scene order.
    cs.edge_post = [idx[("real", v)] for v in cols]
    cs.edge_pre = [pre[("real", v)] for v in cols]
    cs.edge_key = [idx[("real", v)] for v in cols]
    # Top-of-edge inputs for gradients: q_e = pre(binary parent) o e^{Q_sib t_sib} post(sibling).
    top = {}
    for (c, sb, sk, b) in top_inputs:
        top[c] = (pre[b], sb, sk)
    cs.edge_top = [top[("real", v)] for v in cols]
    h = hashlib.sha1()
    for part in (shape, [G], rows, cols, procs, obs_nodes, obs_vars):
        h.update(np.asarray(part, dtype=np.int64).tobytes())
        h.update(b"|")
    cs.topology = h.hexdigest()
    cs.data = None
    if "iid_observations" in od:
        cs.data = pattern_data(cs, od["iid_observations"], compress=compress)
        cs.P = cs.data.cap
        cs.signature = cs.topology + cs.data.digest
        cs.nbytes = 2 * cs.nbuf * cs.P * S * 8
    return cs


class PatternData:
    """Tip partials for one set of observation rows (compressed, optionally padded to ``cap``)."""


def pattern_data(cs, iid, compress=True, capacity=None):
    """Compress rows to unique patterns and build tip partials; pad to ``capacity`` with weight-0,
    all-ones patterns so one instance can serve every scene of a bucket."""
    iid = np.asarray(iid, dtype=np.int64)
    if iid.ndim != 2 or iid.shape[1] != len(cs.obs_nodes) or iid.shape[0] == 0:
        raise UnsupportedScene("observation rows do not match the observed nodes")
    for j, a in enumerate(cs.obs_vars):
        col = iid[:, j]
        if np.any(col < -1) or np.any(col >= cs.shape[a]):
            raise UnsupportedScene("observed state out of range")
    if compress:
        uniq, inv, w = np.unique(iid, axis=0, return_inverse=True, return_counts=True)
        inv = np.asarray(inv).reshape(-1)
    else:
        uniq, inv, w = iid, np.arange(len(iid)), np.ones(len(iid), dtype=np.int64)
    P = len(uniq)
    cap = P if capacity is None else int(capacity)
    if cap < P:
        raise ValueError("pattern capacity %d < %d patterns" % (cap, P))
    d = PatternData()
    d.nobs, d.P, d.cap, d.inv = len(iid), P, cap, inv
    d.weights = f64(np.concatenate([w, np.zeros(cap - P)]))
    coords = np.unravel_index(np.arange(cs.S), cs.shape)
    masks = {}
    for j, (n, a) in enumerate(zip(cs.obs_nodes, cs.obs_vars)):
        vals = uniq[:, j]
        m = (coords[a][:, None] == vals[None, :]) | (vals[None, :] == -1)
        masks[n] = masks[n] * m if n in masks else m.astype(float)
    d.tip_partials = []
    for b in cs.tips:
        a = np.ones((cap, cs.S))
        if b[1] in masks:
            a[:P] = masks[b[1]].T
        d.tip_partials.append((cs.idx[b], f64(a)))
    h = hashlib.sha1(np.ascontiguousarray(iid).tobytes())
    h.update(np.asarray([int(compress), cap], dtype=np.int64).tobytes())
    d.digest = h.hexdigest()
    return d


def scene_values(scene, cs):
    """Generators (triplets per process), root distribution, edge lengths."""
    trips = [generator_triplets(cs.shape, p) for p in scene["process_definitions"]]
    rp = scene["root_prior"]
    distn = np.zeros(cs.S)
    st = np.asarray(rp["states"])
    pr = np.asarray(rp["probabilities"], dtype=float)
    if st.ndim != 2 or st.shape[1] != len(cs.shape) or len(pr) != len(st):
        raise UnsupportedScene("root prior does not match the state space shape")
    np.add.at(distn, np.ravel_multi_index(st.T, cs.shape), pr)
    rates = f64(scene["tree"]["edge_rate_scaling_factors"])
    if np.any(~np.isfinite(rates)) or np.any(rates < 0):
        raise UnsupportedScene("edge rate scaling factors must be finite and >= 0")
    return trips, distn, rates


class BeagleScene:
    """One BEAGLE instance for a compiled topology.

    ``capacity``: patterns per instance (default: the compiled scene's own data). Pattern data of
    any scene with the same topology and at most ``capacity`` patterns can be loaded with
    ``load_patterns``. ``gl``: dict(m=<parameters>, K=6, J_max=4, rho0=2.0) adds the buffers, keys
    and spare slots for exact parameter gradients (``gl_gradient``).
    """

    def __init__(self, lib, cs, resource="cpu", device=None, threads=1, capacity=None, gl=None):
        self.cs = cs
        P = capacity or (cs.data.cap if cs.data is not None else None)
        if P is None:
            raise ValueError("no pattern capacity: compile with data or pass capacity")
        self.P = P
        self.gl = None
        partials, keys, spare = cs.nbuf, cs.keys, 0
        if gl:
            # K nodes on single-panel edges (t ||B|| <= rho0); K_graded nodes per panel on graded edges.
            # Defaults reach ~1e-13 on short edges and <= 5e-11 relative at log tau = 4..7 with
            # 10x branch lengths (tests/test_beagle_l3.py).
            # K="auto": on a single-panel edge use the smallest K in [K_min, K_graded] whose a-priori
            # bound (K!)^4 (2 t||B||)^(2K) / ((2K+1) ((2K)!)^3) is below ``tol`` (plan §6.2).
            g = dict(GL_DEFAULTS)
            g.update(gl)
            J = g["J_max"]
            Ks = list(range(g["K_min"], g["K_graded"] + 1)) if g["K"] == "auto" else [g["K"], g["K_graded"]]
            K = max(Ks)
            g["rules"] = {}
            for k in set(Ks):
                x, w = np.polynomial.legendre.leggauss(k)
                g["rules"][k] = ((x + 1.0) / 2.0, w / 2.0)
            b = cs.nbuf
            g["QT"] = b                         # top vector of the current edge
            g["BP"] = b + 1                     # boundary posts 1..J-1   (J_max + 1 slots)
            g["BR"] = b + 1 + (J + 1)           # boundary pres  1..J-1   (J_max + 1 slots)
            g["CP"] = b + 1 + 2 * (J + 1)       # panel post chain        (K)
            g["CR"] = g["CP"] + K               # panel pre chain         (K)
            partials = g["CR"] + K
            g["kB"] = cs.keys                   # boundary segment keys   (J_max)
            g["kP"] = cs.keys + J               # panel segment keys      (K + 1)
            keys = cs.keys + J + K + 1
            spare = g["m"] * cs.G
            self.gl = g
            self.gl_clipped = 0
        self.inst = ActionInstance(lib, tips=cs.T, partials=partials, states=cs.S, patterns=P,
                                   generators=cs.G, spare_slots=spare, keys=keys,
                                   categories=1, scale_buffers=cs.nscale, resource=resource,
                                   device=device, threads=threads, scratch_keys_from=cs.keys)
        self.inst.set_partials(cs.ONES, np.ones(P * cs.S))
        self.inst.set_categories([1.0], [1.0])
        self.rates = None
        self.data = None
        if cs.data is not None and cs.data.cap == P:
            self.load_patterns(cs.data)

    def load_patterns(self, data):
        if data.cap != self.P:
            raise ValueError("pattern data capacity %d != instance capacity %d" % (data.cap, self.P))
        for b, a in data.tip_partials:
            self.inst.set_tip_partials(b, a)
        self.inst.set_pattern_weights(data.weights)
        self.data = data

    # -------------------------------------------------------------- values and passes
    def set_values(self, trips, distn, rates):
        cs, inst = self.cs, self.inst
        for g, (r, c, v) in enumerate(trips):
            inst.set_generator(g, r, c, v)
        inst.set_state_frequencies(distn)
        self.distn, self.rates, self.trips = distn, np.asarray(rates, dtype=float), trips
        for g in range(cs.G):
            keys = [cs.edge_key[e] for e in range(cs.E) if cs.procs[e] == g]
            lens = [self.rates[e] for e in range(cs.E) if cs.procs[e] == g]
            if g == 0:
                keys.append(cs.Z)
                lens.append(0.0)
            inst.update_transitions(g, keys, lens)

    def loglik(self):
        """(total lnL, per-observation site lnL in the scene's original row order)."""
        cs, inst = self.cs, self.inst
        inst.update_partials(cs.post_ops)
        inst.accumulate_scaling(list(range(cs.CUM)), cs.CUM)
        total, site = inst.root_log_likelihood(cs.idx[cs.root_b], cs.CUM)
        return total, site[: self.data.P][self.data.inv]

    def pre_pass(self):
        cs = self.cs
        self.inst.set_partials(cs.root_pre, np.tile(self.distn, self.P))
        self.inst.update_pre_partials(cs.pre_ops)

    def edge_derivatives(self, per_site=False):
        """d lnL / d log t_e (jsonctmctree's ``deri``) for every edge, optionally per observation."""
        cs, inst = self.cs, self.inst
        self.pre_pass()
        tot, per = inst.edge_derivatives(cs.edge_post, cs.edge_pre, cs.procs, per_pattern=per_site)
        tot = self.rates * tot
        if per_site:
            per = (self.rates[:, None] * per[:, : self.data.P])[:, self.data.inv].T
        return tot, per

    # -------------------------------------------------------------- Gauss-Legendre gradients
    @staticmethod
    def _norm_inf_shifted(trip, S):
        """||Q - mu I||_inf with mu = trace/S, from generator triplets (diagonal = -row sums)."""
        r, c, v = trip
        d = np.zeros(S)
        diag = r == c
        np.add.at(d, r[diag], v[diag])
        mu = d.mean()
        return float(np.max(-d + np.abs(d - mu))) if S else 0.0

    def _panels(self, t, norm):
        bp, clipped = gl_panels(t, norm, self.gl)
        self.gl_clipped += clipped
        return bp

    def _k_single(self, L):
        return gl_k_single(L, self.gl)

    def gl_gradient(self, D, dpi=None, per_site=False, edges=None):
        """Exact-to-quadrature d lnL / d theta_k for k = 0..m-1.

        ``edges``: integrate only these edges and skip the root-prior term when it is given (the
        caller splits one gradient over several instances and adds dpi once).

        D[k][g]: triplets of dQ_g/dtheta_k (diagonal included) or None; dpi[k]: d pi/d theta_k or None.
        Needs ``loglik`` at the current values (the pre-order pass is run here if needed).
        Each edge is integrated over graded panels with K Gauss-Legendre nodes; boundary vectors
        and panel chains live in scratch buffers reused across edges and panels.
        Returns (grad[m], per-observation scores[nobs, m] or None).
        """
        cs, G, inst = self.cs, self.gl, self.inst
        if G is None:
            raise RuntimeError("BeagleScene was created without gl=...")
        m = len(D)
        if m > G["m"]:
            raise ValueError("more parameters than spare slots")
        if inst.pre_epoch != inst.values_epoch:
            self.pre_pass()
        slot = {}
        for k in range(m):
            for g in range(cs.G):
                if D[k][g] is not None and len(D[k][g][2]):
                    slot[(k, g)] = cs.G + k * cs.G + g
                    inst.set_differential(slot[(k, g)], *D[k][g])
        norms = [self._norm_inf_shifted(t, cs.S) for t in self.trips]
        grad = np.zeros(m)
        per_k = np.zeros((m, self.P)) if per_site else None
        Z, ONES = cs.Z, cs.ONES
        self._gl_nodes = [None] * cs.E
        for e in (range(cs.E) if edges is None else edges):
            g = cs.procs[e]
            t = float(self.rates[e])
            ks = [k for k in range(m) if (k, g) in slot]
            if t == 0.0 or not ks:
                continue
            bp = self._panels(t, norms[g])
            J = len(bp) - 1
            K = self._k_single(t * norms[g]) if J == 1 else G["K_graded"]
            s_, w_ = G["rules"][K]
            CP = [G["CP"] + i for i in range(K)]
            CR = [G["CR"] + i for i in range(K)]
            self._gl_nodes[e] = (J * K, bp)
            # top vector, boundary posts (bottom-up) and boundary pres (top-down)
            BP = [cs.edge_post[e]] + [G["BP"] + j for j in range(1, J)]
            BR = [None] + [G["BR"] + j for j in range(1, J)] + [G["QT"]]
            if J > 1:
                inst.update_transitions(g, [G["kB"] + j for j in range(J)], np.diff(bp))
            post = [(G["QT"], NONE, NONE, cs.edge_top[e][0], Z, cs.edge_top[e][1], cs.edge_top[e][2])]
            post += [(BP[j], NONE, NONE, BP[j - 1], G["kB"] + j - 1, ONES, Z) for j in range(1, J)]
            inst.update_partials(post)          # also writes cache(ONES) = ONES, read by pre ops
            if J > 1:
                inst.update_pre_partials([(BR[j], NONE, NONE, BR[j + 1], G["kB"] + j, ONES, Z)
                                          for j in range(J - 1, 0, -1)])
            for j in range(1, J + 1):
                a, b = bp[j - 1], bp[j]
                u = a + (b - a) * s_
                inst.update_transitions(g, [G["kP"] + i for i in range(K + 1)],
                                        np.maximum(np.diff(np.concatenate([[a], u, [b]])), 0.0))
                inst.update_partials([(CP[0], NONE, NONE, BP[j - 1], G["kP"], ONES, Z)] +
                                     [(CP[i], NONE, NONE, CP[i - 1], G["kP"] + i, ONES, Z) for i in range(1, K)])
                inst.update_pre_partials([(CR[K - 1], NONE, NONE, BR[j], G["kP"] + K, ONES, Z)] +
                                         [(CR[i], NONE, NONE, CR[i + 1], G["kP"] + i + 1, ONES, Z)
                                          for i in range(K - 2, -1, -1)])
                wts = (b - a) * w_
                tot, per = inst.edge_derivatives(CP * len(ks), CR * len(ks),
                                                 [slot[(k, g)] for k in ks for _ in range(K)],
                                                 per_pattern=per_site)
                tot = tot.reshape(len(ks), K)
                grad[ks] += tot @ wts
                if per_site:
                    per_k[ks] += np.einsum("i,kip->kp", wts, per.reshape(len(ks), K, self.P))
        if edges is None and dpi is not None and any(d is not None for d in dpi):
            rp = inst.get_partials(cs.idx[cs.root_b])[0]
            den = rp @ self.distn
            w = self.data.weights
            for k in range(m):
                if dpi[k] is not None:
                    r = (rp @ dpi[k]) / den
                    grad[k] += float(w @ r)
                    if per_site:
                        per_k[k] += r
        if per_site:
            per_k = per_k[:, : self.data.P][:, self.data.inv].T
        return grad, per_k

    def edge_vectors(self, edges=None):
        """(q_top, p_bottom) arrays [P x S] for each edge (for the block-generator oracle)."""
        cs = self.cs
        if self.gl is None:
            raise RuntimeError("needs gl buffers")
        if self.inst.pre_epoch != self.inst.values_epoch:
            self.pre_pass()
        out = []
        for e in (range(cs.E) if edges is None else edges):
            self.inst.update_partials([(self.gl["QT"], NONE, NONE, cs.edge_top[e][0], cs.Z, cs.edge_top[e][1],
                                        cs.edge_top[e][2])])
            out.append((self.inst.get_partials(self.gl["QT"])[0], self.inst.get_partials(cs.edge_post[e])[0]))
        return out


GL_DEFAULTS = dict(K="auto", K_min=2, K_graded=8, J_max=32, rho0=2.0, grade=2.0, tol=1e-13)


def gl_panels(t, norm, G):
    """Panel breakpoints on [0, t] and whether the J_max cap clipped them: one panel if
    t*norm <= rho0, else panels graded geometrically (ratio ``grade``) from both ends, the first of
    width rho0/norm. The integrand q'e^{Q(t-u)} D e^{Qu} p has boundary layers of width ~1/||B||
    at both ends."""
    L = t * norm
    if L <= G["rho0"]:
        return np.array([0.0, t]), 0
    h0 = G["rho0"] / norm
    left = [0.0]
    w = h0
    while left[-1] + w < t / 2:
        left.append(left[-1] + w)
        w *= G["grade"]
    bp = np.array(left + [t / 2] + [t - v for v in reversed(left)])
    if len(bp) - 1 > G["J_max"]:
        # keep the graded ends, merge the middle panels
        keep = G["J_max"] // 2
        bp = np.concatenate([bp[: keep + 1], bp[-keep - 1:]]) if keep > 0 else np.array([0.0, t])
        return bp, 1
    return bp, 0


def gl_k_single(L, G):
    """K on a single-panel edge: fixed, or ("auto") the smallest K in [K_min, K_graded] whose
    a-priori bound (K!)^4 (2L)^(2K) / ((2K+1) ((2K)!)^3) is below ``tol``."""
    if G["K"] != "auto":
        return G["K"]
    from math import factorial
    for k in range(G["K_min"], G["K_graded"] + 1):
        b = factorial(k) ** 4 * (2.0 * L) ** (2 * k) / ((2 * k + 1) * float(factorial(2 * k)) ** 3)
        if b <= G["tol"]:
            return k
    return G["K_graded"]


def gl_nodes(t, norm, gl=None):
    """(panels J, nodes per panel K) that ``gl_gradient`` uses for an edge of length t."""
    G = dict(GL_DEFAULTS)
    G.update(gl or {})
    if t == 0.0:
        return 0, 0
    bp, _ = gl_panels(t, norm, G)
    J = len(bp) - 1
    return J, (gl_k_single(t * norm, G) if J == 1 else G["K_graded"])
