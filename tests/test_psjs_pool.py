"""M2: PS-IGC objective on resident workers, engine selectable (jsonctmctree | beagle).

* QBuilder reproduces PSJSModel's process definitions (structure and values) for every codon-site
  pair, with IGC on and off.
* Both pool engines reproduce the legacy per-distance lnL and edge derivatives; their central-FD
  model gradients agree with each other and with the legacy forward differences.
* The merged-N0 gradient (unit weights) matches a central difference of the objective; the legacy
  weights (10, 10/9) do not.
* Installed methods: PSJSGeneconv.loglikelihood_and_gradient / objective_wo_gradient use the
  attached pool; detaching restores the legacy path.
* L2 (slow): the full 1,461-scene objective matches the golden value with both engines.
"""

import json
import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases  # noqa: E402

from IGCexpansion.beagle import psjs  # noqa: E402
from IGCexpansion.beagle.scene import generator_triplets  # noqa: E402

ENGINES = ["jsonctmctree", "beagle"]


def _need(engine):
    if engine == "beagle":
        from beagle_util import lib_or_skip
        lib_or_skip()


@pytest.fixture(scope="module")
def igc_obj():
    d = tempfile.mkdtemp()
    with cases.quiet():
        return cases.psjs_ydr418w_igc(d)


@pytest.fixture(scope="module")
def dist(igc_obj):
    ns = sorted(igc_obj.data.space_list)
    return [ns[0], ns[1], ns[len(ns) // 2], ns[-1]]


@pytest.fixture(scope="module")
def legacy(igc_obj, dist):
    f, g = 0.0, 0.0
    with cases.quiet():
        for n in dist:
            a, b = igc_obj.loglikelihood_and_gradient_for_one_n(n)
            f, g = f + a, g + b
    return -f, -np.asarray(g)


def _dense(t, S=256):
    M = np.zeros((S, S))
    np.add.at(M, (t[0], t[1]), t[2])
    return M


@pytest.mark.parametrize("igc_on", [False, True])
def test_qbuilder_matches_legacy(igc_on):
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w_igc(d) if igc_on else cases.psjs_ydr418w(d)
        obj.unpack_x(np.array(obj.x, float))
    sl = psjs.scene_list(obj)
    qb = psjs.QBuilder(obj.psjsmodel, sorted({cs for _, cs in sl}))
    coef = qb.coefficients()
    assert (coef[2] > 0) == igc_on
    rng = np.random.default_rng(1)
    seen = set()
    for j in list(rng.choice(len(sl), 25, replace=False)) + [0, len(sl) - 1]:
        n, cs = sl[j]
        seen.add(cs)
        with cases.quiet():
            sc = obj.get_scene(n, cs)
        vals = qb.values(cs, n, coef)
        for p, (s, v) in enumerate(zip(qb.struct[cs], vals)):
            pd = sc["process_definitions"][p]
            assert np.array_equal(np.asarray(pd["row_states"]), s["rows"])
            assert np.array_equal(np.asarray(pd["column_states"]), s["cols"])
            np.testing.assert_allclose(v, pd["transition_rates"], rtol=1e-13, atol=1e-15)
            np.testing.assert_allclose(_dense(psjs._generator(s["rows"], s["cols"], v, [4] * 4)),
                                       _dense(generator_triplets([4] * 4, pd)), atol=1e-13)
    assert len(seen) == 9


@pytest.mark.parametrize("engine", ENGINES)
def test_pool_matches_legacy(engine, igc_obj, dist, legacy):
    _need(engine)
    ll0, g0 = legacy
    x = np.array(igc_obj.x, float)
    p = psjs.PSJSPool(igc_obj, engine=engine, processes=2, fd="forward", merged_n0="legacy", distances=dist)
    try:
        ll, g = p.loglikelihood_and_gradient(x)
        assert p.loglikelihood(x) == pytest.approx(ll, rel=1e-13)
    finally:
        p.close()
    assert ll == pytest.approx(ll0, rel=1e-11)
    np.testing.assert_allclose(g[8:], g0[8:], atol=1e-9 * max(1, np.max(np.abs(g0[8:]))))   # edges
    np.testing.assert_allclose(g[:8], g0[:8], atol=1e-3)     # forward differences, step 1e-8


def test_engines_agree_central(igc_obj, dist):
    _need("beagle")
    x = np.array(igc_obj.x, float)
    out = {}
    for e in ENGINES:
        p = psjs.PSJSPool(igc_obj, engine=e, processes=2, distances=dist)
        try:
            out[e] = p.loglikelihood_and_gradient(x)
        finally:
            p.close()
    (la, ga), (lb, gb) = out["jsonctmctree"], out["beagle"]
    assert la == pytest.approx(lb, rel=1e-12)
    np.testing.assert_allclose(ga, gb, rtol=1e-8, atol=1e-6)


def test_merged_n0_fix(igc_obj, dist):
    """d lnL / d x_rate[0] equals d_N0D1 + d_N0out (unit weights), not the legacy 10, 10/9."""
    x = np.array(igc_obj.x, float)
    p = psjs.PSJSPool(igc_obj, engine="jsonctmctree", processes=2, distances=dist)
    q = psjs.PSJSPool(igc_obj, engine="jsonctmctree", processes=2, distances=dist, merged_n0="legacy")
    try:
        i = p.m                                   # first rate parameter: the merged N0 rate
        _, g = p.loglikelihood_and_gradient(x)
        _, gl = q.loglikelihood_and_gradient(x)
        h = 1e-5
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        fd = (p.loglikelihood(xp) - p.loglikelihood(xm)) / (2 * h)
    finally:
        p.close()
        q.close()
    assert g[i] == pytest.approx(fd, rel=1e-6)
    assert abs(gl[i] - fd) > 0.5 * abs(fd)
    np.testing.assert_allclose(np.delete(g, i), np.delete(gl, i), rtol=1e-12)


@pytest.mark.parametrize("engine", ENGINES)
def test_installed_methods(engine, dist):
    _need(engine)
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w_igc(d)
    psjs.install()
    try:
        psjs.attach(obj, engine=engine, processes=2, distances=dist)
        f, g = obj.loglikelihood_and_gradient()
        assert obj.ll == pytest.approx(-f)
        x = np.array(obj.x, float)
        assert obj.objective_wo_gradient(False, x) == pytest.approx(f, rel=1e-12)
        psjs.detach(obj)
        with cases.quiet():
            ll_legacy = sum(obj._loglikelihood_for_one_n(n)[0] for n in dist)
        assert -f == pytest.approx(ll_legacy, rel=1e-11)
    finally:
        psjs.detach(obj)
        psjs.uninstall()


@pytest.mark.slow
@pytest.mark.parametrize("engine", ENGINES)
def test_full_objective_golden(engine):
    _need(engine)
    golden = json.load(open(os.path.join(HERE, "golden", "l2.json")))["values"]["psjs_objective_at_x"]["lnL"]
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w(d)
        obj.unpack_x(np.array(obj.x, float))
    p = psjs.PSJSPool(obj, engine=engine, processes=4)
    try:
        assert len(p.scenes) == 1461
        assert p.loglikelihood(np.array(obj.x, float)) == pytest.approx(golden, rel=1e-11)
    finally:
        p.close()
