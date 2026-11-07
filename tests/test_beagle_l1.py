"""L1: the BEAGLE engine against jsonctmctree (impl_v2) on synthetic and real scenes (plan §11).

Tolerances: site lnL 1e-11 relative; edge derivatives 1e-9 * max(1, |d|).
Synthetic scenes cover unary chains, polytomies, observed internal nodes and root, -1 masks,
zero-length edges, duplicate and zero triplets, several processes, multi-axis state spaces,
shuffled node labels, and a deep tree whose site likelihoods underflow without scaling.

The reference for ``grad`` requests is jsonctmctree's ``deri``: its O(E) ``grad`` pass is wrong
on edges below an observed internal node (``test_jsonctmctree_grad_bug_observed_internal``).
"""

import copy
import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from beagle_util import lib_or_skip  # noqa: E402
import cases  # noqa: E402

from IGCexpansion.beagle import engine  # noqa: E402
from IGCexpansion.beagle.scene import compile_scene  # noqa: E402

RTOL_SITE = 1e-11
ATOL_DERI = 1e-9
PROPS = ["snnlogl", "dnnlogl", "sdnderi", "ddnderi", "sdngrad", "ddngrad"]


@pytest.fixture(autouse=True)
def _beagle_mode():
    lib_or_skip()
    old = engine.get_mode()
    engine.set_mode("beagle")
    yield
    engine.set_mode(old)


def both(scene, props=PROPS, extra=()):
    reqs = [{"property": p} for p in props] + list(extra)
    j = {"scene": scene, "requests": reqs}
    before = engine.STATS["fallback_calls"]
    got = engine.process_json_in(copy.deepcopy(j))
    assert engine.STATS["fallback_calls"] == before, "BEAGLE engine fell back to jsonctmctree"
    jref = copy.deepcopy(j)
    for r in jref["requests"]:
        r["property"] = r["property"].replace("grad", "deri")
    try:
        ref = engine.ORIGINAL(jref)
    except Exception as e:  # jsonctmctree raises for some infeasible inputs
        ref = {"status": "error: %s" % type(e).__name__, "responses": None}
    return reqs, got, ref


def check(scene, props=PROPS, extra=()):
    reqs, got, ref = both(scene, props, extra)
    assert got["status"] == ref["status"]
    if ref["responses"] is None:
        return got, ref
    for r, a, b in zip(reqs, got["responses"], ref["responses"]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        assert a.shape == b.shape, r["property"]
        if r["property"].endswith("logl"):
            np.testing.assert_allclose(a, b, rtol=RTOL_SITE, atol=1e-13, err_msg=r["property"])
        else:
            np.testing.assert_array_less(np.abs(a - b), ATOL_DERI * np.maximum(1.0, np.abs(b)) + 1e-300,
                                         err_msg=r["property"])
    return got, ref


# ----------------------------------------------------------------------------- scenes
def naive_scene():
    """jsonctmctree/tests/test_vs_naive.py: observed internal node, 2 axes, 3 processes."""
    a, b, x = 0.2, 0.3, 0.4
    rs = [[0, 0], [0, 0], [0, 1], [0, 1], [1, 0], [1, 0], [1, 1], [1, 1]]
    cs = [[0, 1], [1, 0], [0, 0], [1, 1], [0, 0], [1, 1], [0, 1], [1, 0]]
    return dict(
        node_count=5, process_count=2, state_space_shape=[2, 2],
        tree=dict(row_nodes=[0, 0, 2, 2], column_nodes=[1, 2, 3, 4],
                  edge_rate_scaling_factors=[1.0, 2.0, 3.0, 4.0], edge_processes=[0, 1, 1, 2]),
        root_prior=dict(states=[[0, 0], [0, 1], [1, 0]], probabilities=[0.25, 0.25, 0.5]),
        process_definitions=[
            dict(row_states=rs, column_states=cs, transition_rates=[a, a, a, b, b, a, b, b]),
            dict(row_states=rs, column_states=cs, transition_rates=[a, a, a + x, b + x, b + x, a + x, b, b]),
            dict(row_states=[[0, 0], [0, 1], [1, 1], [1, 0]], column_states=[[0, 1], [1, 1], [1, 0], [0, 0]],
                 transition_rates=[1, 1, 1, 1])],
        observed_data=dict(nodes=[1, 1, 3, 3, 2, 4], variables=[0, 1, 0, 1, 1, 1],
                           iid_observations=[[0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1], [0, 1, 0, 1, 0, 1],
                                             [1, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]))


def random_scene(seed, deep=False, n_deep=700):
    rng = np.random.default_rng(seed)
    shape = [[4], [3, 2], [2, 2, 2], [5]][seed % 4]
    S = int(np.prod(shape))
    states = np.array(np.unravel_index(np.arange(S), shape)).T
    G = 1 + seed % 3
    procs = []
    for _ in range(G):
        r, c, v = [], [], []
        for i in range(S):
            for j in range(S):
                if i != j and rng.random() < 0.6:
                    r.append(states[i].tolist())
                    c.append(states[j].tolist())
                    v.append(float(rng.uniform(0.05, 1.5)))
        # duplicate and zero triplets
        for _ in range(3):
            k = int(rng.integers(len(v)))
            r.append(r[k])
            c.append(c[k])
            v.append(float(rng.uniform(0, 0.5)))
        r.append(r[0])
        c.append(c[0])
        v.append(0.0)
        procs.append(dict(row_states=r, column_states=c, transition_rates=v))
    if deep:
        n = n_deep                                 # caterpillar
        parent = [None] + [max(0, k - 2 + (k % 2)) for k in range(1, n)]
    else:
        n = int(rng.integers(5, 16))
        parent = [None] + [int(rng.integers(0, k)) for k in range(1, n)]   # unary nodes and polytomies
    label = rng.permutation(n)
    rows = [int(label[parent[k]]) for k in range(1, n)]
    cols = [int(label[k]) for k in range(1, n)]
    t = rng.uniform(0.0, 3.0 if deep else 1.5, size=n - 1)
    if not deep:
        t[rng.random(n - 1) < 0.15] = 0.0
    eproc = rng.integers(0, G, size=n - 1).tolist()
    children = {k: [] for k in range(n)}
    for k in range(1, n):
        children[parent[k]].append(k)
    obs_nodes, obs_vars = [], []
    for k in range(n):
        leaf = not children[k]
        if leaf or (not deep and rng.random() < 0.3):
            for a in range(len(shape)):
                if leaf or rng.random() < 0.6:
                    obs_nodes.append(int(label[k]))
                    obs_vars.append(a)
    nrow = 40 if deep else 25
    pool = int(rng.integers(3, 12))
    base = np.array([[int(rng.integers(-1 if rng.random() < 0.15 else 0, shape[a])) for a in obs_vars]
                     for _ in range(pool)])
    iid = base[rng.integers(0, pool, size=nrow)]                       # duplicate patterns
    support = rng.choice(S, size=int(rng.integers(1, S + 1)), replace=False)
    p = rng.random(len(support)) + 0.05
    return dict(
        node_count=n, process_count=G, state_space_shape=shape,
        tree=dict(row_nodes=rows, column_nodes=cols, edge_rate_scaling_factors=t.tolist(), edge_processes=eproc),
        root_prior=dict(states=states[support].tolist(), probabilities=(p / p.sum()).tolist()),
        process_definitions=procs,
        observed_data=dict(nodes=obs_nodes, variables=obs_vars, iid_observations=iid.tolist()))


# ----------------------------------------------------------------------------- tests
def test_naive_scene_all_supported_properties():
    w = {"observation_reduction": {"observation_indices": [0, 1, 2, 4, 3, 2],
                                   "weights": [0.1, 0.1, 0.2, 0.3, 0.5, 0.8]}}
    check(naive_scene(), extra=[dict(property="wnnlogl", **w), dict(property="wdnderi", **w),
                                dict(property="wdngrad", **w)])


@pytest.mark.parametrize("seed", range(40))
def test_random_scene(seed):
    sc = random_scene(seed)
    got, ref = check(sc)
    # the generator must actually produce the features it claims, across seeds
    cs = compile_scene(sc)
    test_random_scene.features = getattr(test_random_scene, "features", set())
    if any(b[0] == "pseudo" for b in cs.idx):
        test_random_scene.features.add("observed internal")
    if any(b[0] == "virt" for b in cs.idx):
        test_random_scene.features.add("polytomy")
    if any(op[5] == cs.ONES for op in cs.post_ops):
        test_random_scene.features.add("unary")
    if got["status"] == "infeasible":
        test_random_scene.features.add("infeasible")


def test_random_scene_features_covered():
    for seed in range(40):
        test_random_scene(seed)
    assert {"observed internal", "polytomy", "unary"} <= test_random_scene.features


def deep_scene(n=2000, seed=5):
    """Caterpillar with dense HKY-like generators, all states at the root: lnL far below log(1e-308)."""
    sc = random_scene(seed, deep=True, n_deep=n)
    rng = np.random.default_rng(seed)
    S = 4
    st = [[i] for i in range(S)]
    r, c, v = zip(*[([i], [j], float(rng.uniform(0.2, 1.2))) for i in range(S) for j in range(S) if i != j])
    sc.update(state_space_shape=[S], process_count=1,
              process_definitions=[dict(row_states=list(r), column_states=list(c), transition_rates=list(v))],
              root_prior=dict(states=st, probabilities=[0.1, 0.2, 0.3, 0.4]))
    sc["tree"]["edge_processes"] = [0] * len(sc["tree"]["edge_processes"])
    leaves = sorted(set(sc["observed_data"]["nodes"]))     # deep scenes observe leaves only
    sc["observed_data"] = dict(nodes=leaves, variables=[0] * len(leaves),
                               iid_observations=rng.integers(-1, S, size=(30, len(leaves))).tolist())
    return sc


def scaled_pruning(sc):
    """Per-site lnL by dense expm pruning with per-node rescaling (no underflow)."""
    import scipy.linalg
    S = sc["state_space_shape"][0]
    Q = np.zeros((S, S))
    p = sc["process_definitions"][0]
    for (i,), (j,), v in zip(p["row_states"], p["column_states"], p["transition_rates"]):
        Q[i, j] += v
    Q -= np.diag(Q.sum(axis=1))
    tr = sc["tree"]
    kids = {}
    for u, v_, t in zip(tr["row_nodes"], tr["column_nodes"], tr["edge_rate_scaling_factors"]):
        kids.setdefault(u, []).append((v_, t))
    root = (set(tr["row_nodes"]) - set(tr["column_nodes"])).pop()
    od = sc["observed_data"]
    iid = np.asarray(od["iid_observations"])
    obs = {}
    for j, nd in enumerate(od["nodes"]):
        m = (np.arange(S)[:, None] == iid[None, :, j]) | (iid[None, :, j] == -1)
        obs[nd] = obs.get(nd, np.ones((S, len(iid)))) * m
    order, stack = [], [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack += [v_ for v_, _ in kids.get(u, [])]
    vec, logs = {}, {}
    for u in reversed(order):
        a = obs.get(u, np.ones((S, len(iid)))).copy()
        lg = np.zeros(len(iid))
        for v_, t in kids.get(u, []):
            a *= scipy.linalg.expm(Q * t) @ vec[v_]
            lg += logs[v_]
        mx = a.max(axis=0)
        vec[u], logs[u] = a / mx, lg + np.log(mx)
    pi = np.zeros(S)
    for (i,), pr in zip(sc["root_prior"]["states"], sc["root_prior"]["probabilities"]):
        pi[i] = pr
    return np.log(pi @ vec[root]) + logs[root]


def test_deep_tree_underflow():
    """2000-node caterpillar: site likelihoods far below 1e-308. jsonctmctree (unscaled) calls it
    infeasible; BEAGLE with post-order scaling matches a log-scaled dense pruning, and its edge
    derivatives (unscaled pre-order pass) match central differences of its own lnL."""
    sc = deep_scene()
    ref = scaled_pruning(sc)
    assert ref.max() < np.log(1e-308)
    assert engine.ORIGINAL({"scene": sc, "requests": [{"property": "snnlogl"}]})["status"] == "infeasible"
    out = engine.process_json_in({"scene": sc, "requests": [{"property": "dnnlogl"}, {"property": "sdnderi"}]})
    assert out["status"] == "feasible"
    np.testing.assert_allclose(out["responses"][0], ref, rtol=1e-11)
    d = np.asarray(out["responses"][1])
    assert np.all(np.isfinite(d))
    t = np.array(sc["tree"]["edge_rate_scaling_factors"])
    for e in [0, len(t) // 2, len(t) - 1]:                 # top, middle, deepest edges
        a, b = copy.deepcopy(sc), copy.deepcopy(sc)
        h = 1e-5 * max(t[e], 1e-3)
        a["tree"]["edge_rate_scaling_factors"][e] += h
        b["tree"]["edge_rate_scaling_factors"][e] -= h
        la = engine.process_json_in({"scene": a, "requests": [{"property": "snnlogl"}]})["responses"][0]
        lb = engine.process_json_in({"scene": b, "requests": [{"property": "snnlogl"}]})["responses"][0]
        assert d[e] == pytest.approx(t[e] * (la - lb) / (2 * h), rel=1e-5, abs=1e-6)


def _igc(name):
    with tempfile.TemporaryDirectory() as d:
        with cases.quiet():
            if name == "psjs":
                obj = cases.psjs_ydr418w(d)
                obj.unpack_x(np.array(obj.x, float))
                return [obj.get_scene(n, pair) for n, pair in cases.PSJS_SCENES]
            obj = {"hky": cases.hky_ylr406c, "mg94": cases.mg94_ylr406c}[name](d)
            return [obj.get_scene()]


@pytest.mark.parametrize("name", ["hky", "mg94", "psjs"])
def test_igc_scenes(name):
    for sc in _igc(name):
        check(sc)


def test_jsonctmctree_grad_bug_observed_internal():
    """jsonctmctree's sdngrad disagrees with finite differences below an observed internal node;
    sdnderi and BEAGLE agree with them. (Leaf-only IGC scenes are not affected.)"""
    sc = naive_scene()
    J = engine.ORIGINAL

    def ll(s):
        return J({"scene": s, "requests": [{"property": "snnlogl"}]})["responses"][0]
    t = np.array(sc["tree"]["edge_rate_scaling_factors"])
    fd = []
    for e in range(len(t)):
        h = 1e-6 * t[e]
        a, b = copy.deepcopy(sc), copy.deepcopy(sc)
        a["tree"]["edge_rate_scaling_factors"][e] += h
        b["tree"]["edge_rate_scaling_factors"][e] -= h
        fd.append(t[e] * (ll(a) - ll(b)) / (2 * h))
    r = J({"scene": sc, "requests": [{"property": "sdnderi"}, {"property": "sdngrad"}]})["responses"]
    got = engine.process_json_in({"scene": sc, "requests": [{"property": "sdngrad"}]})["responses"][0]
    np.testing.assert_allclose(r[0], fd, atol=1e-6)
    np.testing.assert_allclose(got, fd, atol=1e-6)
    assert np.max(np.abs(np.asarray(r[1]) - fd)) > 0.1      # edges 2 and 3 (below observed node 2)


def _no_outgroup(d):
    """HKY with no outgroup: the root is the duplication node with a single child (unary root)."""
    nwk = os.path.join(d, "no_outgroup.newick")
    with open(nwk, "w") as f:
        f.write("(((((cerevisiae,paradoxus)N5,mikatae)N4,kudriavzevii)N3,bayanus)N2,castellii)N1;\n")
    aln = os.path.join(d, "no_outgroup.fasta")
    recs = open(os.path.join(cases.T, "YLR406C_YDL075W_test_input.fasta")).read().split(">")[1:]
    with open(aln, "w") as f:
        f.write("".join(">" + r for r in recs if not r.startswith("kluyveri")))
    obj = cases._recodon("HKY", nwk, aln, ["YLR406C", "YDL075W"], d)
    cases._set_x(obj, kappa=2.5, tau=1.5, blen_seed=7)
    return obj


@pytest.mark.parametrize("name", ["hky", "mg94", "no_outgroup", "psjs"])
def test_shadow_model_gradients(name):
    """IGCexpansion's own lnL + gradient methods (snnlogl + sdnderi) agree under shadow mode."""
    engine.install("shadow")              # rebinds jsonctmctree.interface.process_json_in
    engine.SHADOW_MAX.clear()
    calls0, bad0 = engine.STATS["beagle_calls"], engine.STATS["shadow_mismatch"]
    try:
        _run_model(name)
    finally:
        engine.uninstall()
    assert engine.STATS["beagle_calls"] > calls0
    assert engine.STATS["shadow_mismatch"] == bad0
    assert any("sdnderi" in k for k in engine.SHADOW_MAX)
    assert max(m for _, m in engine.SHADOW_MAX.values()) < 1e-9


def _run_model(name):
    with tempfile.TemporaryDirectory() as d, cases.quiet():
        if name == "psjs":
            obj = cases.psjs_ydr418w(d)
            obj.unpack_x(np.array(obj.x, float))
            for n in sorted(obj.data.space_list)[:3]:
                obj.loglikelihood_and_gradient_for_one_n(n)
        else:
            obj = {"hky": cases.hky_ylr406c, "mg94": cases.mg94_ylr406c, "no_outgroup": _no_outgroup}[name](d)
            if name == "no_outgroup":
                cs = compile_scene(obj.get_scene())
                assert any(op[5] == cs.ONES and op[0] == cs.idx[cs.root_b] for op in cs.post_ops)
            obj.loglikelihood_and_gradient()
