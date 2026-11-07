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

    # Observations, compressed to unique patterns.
    od = scene["observed_data"]
    obs_nodes = [int(n) for n in od["nodes"]]
    obs_vars = [int(a) for a in od["variables"]]
    iid = np.asarray(od["iid_observations"], dtype=np.int64)
    if iid.ndim != 2 or iid.shape[1] != len(obs_nodes) or len(obs_vars) != len(obs_nodes) or iid.shape[0] == 0:
        raise UnsupportedScene("observed_data arrays are inconsistent or empty")
    if any(n not in nodes for n in obs_nodes) or any(a < 0 or a >= len(shape) for a in obs_vars):
        raise UnsupportedScene("observation on an unknown node or axis")
    for j, a in enumerate(obs_vars):
        col = iid[:, j]
        if np.any(col < -1) or np.any(col >= shape[a]):
            raise UnsupportedScene("observed state out of range")
    if compress:
        uniq, inv, w = np.unique(iid, axis=0, return_inverse=True, return_counts=True)
        inv = np.asarray(inv).reshape(-1)
    else:
        uniq, inv, w = iid, np.arange(len(iid)), np.ones(len(iid), dtype=np.int64)
    cs.nobs, cs.P, cs.inv, cs.weights = len(iid), len(uniq), inv, f64(w)
    coords = np.unravel_index(np.arange(S), shape)
    masks = {}
    for j, (n, a) in enumerate(zip(obs_nodes, obs_vars)):
        vals = uniq[:, j]
        m = (coords[a][:, None] == vals[None, :]) | (vals[None, :] == -1)
        masks[n] = masks[n] * m if n in masks else m.astype(float)

    # Binary operation tree: ('real', n), ('virt', k), ('pseudo', n).
    leaves = sorted(n for n in nodes if not children[n])
    internal = sorted(n for n in nodes if children[n])
    pseudo = sorted(n for n in internal if n in masks)
    tips = [("real", n) for n in leaves] + [("pseudo", n) for n in pseudo]
    bchildren = {}           # binary node -> [(child, is_edge_with_key_Z)]
    virt = []
    for n in internal:
        ch = [(("real", c), False) for c in children[n]]
        if n in masks:
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
    cs.tip_partials = []
    for b in tips:
        m = masks.get(b[1]) if b[0] == "pseudo" or b[1] in masks else None
        a = np.ones((S, cs.P)) if m is None else m
        cs.tip_partials.append((idx[b], f64(a.T)))
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
    pre_ops = []
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
    cs.pre_ops = pre_ops
    cs.nscale = cs.CUM + 1 + len(pre_ops)
    cs.root_pre = pre[cs.root_b]
    # Real edges in scene order.
    cs.edge_post = [idx[("real", v)] for v in cols]
    cs.edge_pre = [pre[("real", v)] for v in cols]
    cs.edge_key = [idx[("real", v)] for v in cols]
    cs.zero_keys = sorted({cs.Z})
    h = hashlib.sha1()
    for part in (shape, [G], rows, cols, procs, obs_nodes, obs_vars, [int(compress)]):
        h.update(np.asarray(part, dtype=np.int64).tobytes())
        h.update(b"|")
    h.update(np.ascontiguousarray(iid).tobytes())
    cs.signature = h.hexdigest()
    cs.nbytes = 2 * cs.nbuf * cs.P * S * 8
    return cs


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
    """One BEAGLE instance for a compiled scene."""

    def __init__(self, lib, cs, resource="cpu", device=None, threads=1, spare_slots=0):
        self.cs = cs
        self.inst = ActionInstance(lib, tips=cs.T, partials=cs.nbuf, states=cs.S, patterns=cs.P,
                                   generators=cs.G, spare_slots=spare_slots, keys=cs.keys,
                                   categories=1, scale_buffers=cs.nscale, resource=resource,
                                   device=device, threads=threads)
        for b, a in cs.tip_partials:
            self.inst.set_tip_partials(b, a)
        self.inst.set_partials(cs.ONES, np.ones(cs.P * cs.S))
        self.inst.set_pattern_weights(cs.weights)
        self.inst.set_categories([1.0], [1.0])
        self.rates = None

    def set_values(self, trips, distn, rates):
        cs, inst = self.cs, self.inst
        for g, (r, c, v) in enumerate(trips):
            inst.set_generator(g, r, c, v)
        inst.set_state_frequencies(distn)
        self.distn, self.rates = distn, rates
        for g in range(cs.G):
            keys = [cs.edge_key[e] for e in range(cs.E) if cs.procs[e] == g]
            lens = [rates[e] for e in range(cs.E) if cs.procs[e] == g]
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
        return total, site[cs.inv]

    def edge_derivatives(self, per_site=False):
        """d lnL / d log t_e (jsonctmctree's ``deri``) for every edge, optionally per observation."""
        cs, inst = self.cs, self.inst
        inst.set_partials(cs.root_pre, np.tile(self.distn, cs.P))
        inst.update_pre_partials(cs.pre_ops)
        tot, per = inst.edge_derivatives(cs.edge_post, cs.edge_pre, cs.procs, per_pattern=per_site)
        tot = self.rates * tot
        if per_site:
            per = (self.rates[:, None] * per)[:, cs.inv].T
        return tot, per
