"""M4: IS-IGC (ReCodonGeneconv, JointAnalysis) on BEAGLE with exact gradients (plan L2-L4).

* IsIgcBuilder reproduces the object's own process definitions and root prior (1e-13) for HKY,
  MG94, MG94 + IGC_Omega, MG94 + Tau_Omega, and with Force.
* Exact gradients: GL with the analytic complex-step Jacobian against the block-generator
  oracle (1e-10), and against IGCexpansion's central differences (step ~1e-6).
* Golden lnL values (L2) through the patched class.
* RecodonPool (whole pairs and one pair split by edges) equals in-process evaluation.
* JointAnalysis objective equals the legacy one; L4: HKY fits on BEAGLE reach the legacy
  optimum.
"""

import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from beagle_util import lib_or_skip  # noqa: E402
import cases  # noqa: E402

T = cases.T


@pytest.fixture(scope="module")
def lib():
    return lib_or_skip()


@pytest.fixture
def recodon(lib):
    from IGCexpansion.beagle import recodon as rc
    yield rc
    rc.uninstall()


def _obj(model, aln="YLR406C_YDL075W_test_input.fasta", paralog=("YLR406C", "YDL075W"), **kw):
    from IGCexpansion.CodonGeneconv import ReCodonGeneconv
    d = tempfile.mkdtemp()
    with cases.quiet():
        o = ReCodonGeneconv(T + "/YeastTree.newick", T + "/" + aln, list(paralog), Model=model, clock=False,
                            save_path=d + "/", save_name=os.path.join(d, "s.txt"), **kw)
        cases._set_x(o, kappa=3.0, tau=2.0, omega=0.3 if model == "MG94" else None)
    return o


def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - b) / np.maximum(1.0, np.abs(b))))


@pytest.mark.parametrize("model,kw", [
    ("HKY", {}),
    ("MG94", {}),
    ("MG94", {"IGC_Omega": 0.5}),
    ("MG94", {"Tau_Omega": 0.7}),
    ("HKY", {"Force": {3: 2.5}}),
])
def test_builder_matches_scene(lib, model, kw):
    from IGCexpansion.beagle.recodon import IsIgcBuilder
    from IGCexpansion.beagle.scene import compile_scene
    o = _obj(model, **kw)
    cs = compile_scene(o.get_scene())
    assert IsIgcBuilder(o).check(o, cs) < 1e-13


@pytest.mark.parametrize("model,kw", [("HKY", {}), ("HKY", {"Force": {3: 2.5}}),
                                      pytest.param("MG94", {"IGC_Omega": 0.5}, marks=pytest.mark.slow)])
def test_gl_vs_block_and_fd(recodon, model, kw):
    from IGCexpansion.beagle._lib import get_lib
    from IGCexpansion.beagle.blockgen import block_gradient
    o = _obj(model, **kw)
    ll, g, _ = recodon.evaluate(o)
    st = recodon._STATE[o]
    D, dpi = st.jac.at(np.asarray(o.x_process, float))
    ref = block_gradient(get_lib(), st.bs, D, dpi)
    assert _rel(g[st.free], ref) < 1e-10
    for i in set(range(st.m)) - set(st.free):
        assert g[i] == 0.0
    recodon.detach(o)
    with cases.quiet():
        f0, g0 = o.loglikelihood_and_gradient()         # legacy class: jsonctmctree + central FD
    assert ll == pytest.approx(-f0, rel=1e-12)
    np.testing.assert_allclose(g, -np.asarray(g0), rtol=1e-6, atol=1e-6)


def test_per_site_scores_sum(recodon):
    o = _obj("HKY")
    ll, g, per = recodon.evaluate(o, per_site=True)
    np.testing.assert_allclose(per.sum(axis=0), g, rtol=1e-11, atol=1e-9)


def test_golden_through_patched_class(recodon):
    import json
    gold = json.load(open(os.path.join(HERE, "golden", "l2.json")))
    o = cases.hky_ylr406c(tempfile.mkdtemp())
    recodon.attach(o)
    with cases.quiet():
        f, g = o.loglikelihood_and_gradient()
        ll, ed = o._loglikelihood(edge_derivative=True)
    ref = gold["values"]["hky_ylr406c"]
    ref = ref["lnL"]
    assert -f == pytest.approx(ref, rel=1e-12)
    assert ll == -f
    np.testing.assert_allclose(ed, -g[-len(ed):], rtol=1e-14)


def test_pool_matches_in_process(recodon):
    objs = [_obj("MG94"), _obj("HKY")]
    ref = [recodon.evaluate(o) for o in objs]
    for k in (1, 3, 6):
        pool = recodon.RecodonPool(objs, processes=k)
        try:
            if k == 6:
                assert any(u[0][1] is not None for u in pool.units)       # the MG94 pair was split
            ll, gs = pool.evaluate([o.x for o in objs])
            ll2, _ = pool.evaluate([o.x for o in objs], grad=False)
        finally:
            pool.close()
        for (l0, g0, _), l, l2, g in zip(ref, ll, ll2, gs):
            assert l == pytest.approx(l0, rel=1e-12) and l2 == pytest.approx(l0, rel=1e-12)  # GPU: not bitwise
            assert _rel(g, g0) < 1e-11


def test_plan_units():
    from IGCexpansion.beagle.recodon import plan_units
    costs = [np.array([5.0, 1, 1, 3]), np.array([0.1, 0.1])]
    u = plan_units(costs, 3)
    assert sum(len(v) for v in u) == 3
    split = [x for v in u for x in v if x[0] == 0]
    assert sorted(e for x in split for e in x[1]) == [0, 1, 2, 3]
    assert sum(x[2] for x in split) == 1                                   # one primary unit
    u = plan_units(costs, 1)
    assert u == [[(0, None, True), (1, None, True)]] or u == [[(1, None, True), (0, None, True)]]


def _joint(Model="HKY", Shared=None):
    from IGCexpansion.JointAnalysis import JointAnalysis
    d = tempfile.mkdtemp() + "/"
    with cases.quiet():
        ja = JointAnalysis([T + "/YLR406C_YDL075W_test_input.fasta", T + "/YDR418W_YEL054C_old_input.fasta"],
                           T + "/YeastTree.newick", [["YLR406C", "YDL075W"], ["YDR418W", "YEL054C"]],
                           Model=Model, Shared=Shared, save_path=d)
    return ja


def test_joint_objective(recodon):
    ja = _joint(Shared=[4])                        # shared tau
    x = np.array(ja.x)
    with cases.quiet():
        f0, g0 = ja.objective_and_gradient(x)      # legacy (serial)
    recodon.attach(ja, processes=2)
    try:
        with cases.quiet():
            f1, g1 = ja.objective_and_gradient_multi_threaded(x)
            f2, g2 = ja.objective_and_gradient(x)
        assert f1 == pytest.approx(f0, rel=1e-12) and f2 == f1
        np.testing.assert_allclose(g1, g0, rtol=1e-6, atol=1e-6)
        np.testing.assert_array_equal(g1, g2)
        assert len(ja.ll) == 2
    finally:
        recodon.detach(ja)


@pytest.mark.slow
def test_l4_hky_fit(recodon):
    """L4: an HKY ReCodon fit on BEAGLE reaches the legacy (jsonctmctree, FD) optimum."""
    fits = {}
    for engine in ("legacy", "beagle"):
        o = _obj("HKY")
        if engine == "beagle":
            recodon.attach(o)
        with cases.quiet():
            r = o.get_mle(display=False)
        fits[engine] = (r.fun, np.array(r.x))
        recodon.detach(o)
    (f0, x0), (f1, x1) = fits["legacy"], fits["beagle"]
    assert f1 <= f0 + 1e-6
    assert abs(f1 - f0) < 1e-4
    np.testing.assert_allclose(x1, x0, atol=2e-3)
