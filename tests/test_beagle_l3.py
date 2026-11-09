"""L3: exact parameter gradients (Gauss-Legendre on BEAGLE's edge-derivative kernel) against the
block-generator oracle (plan §6.2, §11).

* HKY: K sweep converges to the oracle (K=6 to ~1e-13); per-site scores sum to the gradient;
  plumbing check (D = Q_p reproduces the edge derivatives); agreement with IGCexpansion's FD.
* Random scenes (several processes, polytomies, observed internal nodes) with random D and dpi.
* Composite panels (J_e > 1) and the J_max clip counter.
* Long edges and large tau (log tau = 4 and 7, branch lengths x10).
* PS-IGC scenes with the analytic complex-step differential matrices, which also match the
  generic scene-FD differential matrices.
* MG94 (3721 states) log tau (slow).
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
from test_beagle_l1 import random_scene  # noqa: E402

from IGCexpansion.beagle.blockgen import block_gradient  # noqa: E402
from IGCexpansion.beagle.differentials import scene_fd_differentials  # noqa: E402
from IGCexpansion.beagle.scene import BeagleScene, compile_scene, scene_values  # noqa: E402

RTOL = 1e-10


def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - b) / np.maximum(1.0, np.abs(b))))


def _scene(lib, sc, m, **gl):
    cs = compile_scene(sc)
    g = dict(m=m)
    g.update(gl)
    bs = BeagleScene(lib, cs, gl=g, resource=os.environ.get("IGC_BEAGLE_RESOURCE", "cpu"))
    bs.set_values(*scene_values(sc, cs))
    bs.loglik()
    return cs, bs


def _model(name, d, **kw):
    with cases.quiet():
        obj = {"hky": cases.hky_ylr406c, "mg94": cases.mg94_ylr406c}[name](d)
        if kw:
            x = np.array(obj.x, float)
            k = len(obj.edge_to_blen)
            npr = len(x) - k
            if "log_tau" in kw:
                x[npr - 1] = kw["log_tau"]
            if "blen_scale" in kw:
                x[npr:] += np.log(kw["blen_scale"])
            obj.update_by_x(x)
    return obj


def _jac(obj, cs, free):
    def scene_at(xx):
        with cases.quiet():
            obj.update_by_x(xx)
        return obj.get_scene()
    return scene_fd_differentials(scene_at, np.array(obj.x, float), free, cs)


@pytest.fixture(scope="module")
def lib():
    return lib_or_skip()


@pytest.fixture(scope="module")
def hky(lib):
    d = tempfile.mkdtemp()
    obj = _model("hky", d)
    sc = obj.get_scene()
    m = len(obj.x) - len(obj.edge_to_blen)
    cs = compile_scene(sc)
    D, dpi = _jac(obj, cs, list(range(m)))
    return obj, sc, m, D, dpi


def test_hky_k_sweep(lib, hky):
    obj, sc, m, D, dpi = hky
    errs = []
    ref = None
    for K in (2, 4, 6, 8):
        cs, bs = _scene(lib, sc, m, K=K, K_graded=K, rho0=1e9)
        g, per = bs.gl_gradient(D, dpi, per_site=True)
        if ref is None:
            ref = block_gradient(lib, bs, D, dpi)
        errs.append(_rel(g, ref))
        np.testing.assert_allclose(per.sum(axis=0), g, rtol=1e-12, atol=1e-10)
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 1e-12 and errs[3] < 1e-12
    # against IGCexpansion's own central differences of the objective (step error ~1e-7)
    with cases.quiet():
        _, g_fd = obj.loglikelihood_and_gradient()
    np.testing.assert_allclose(g, -np.asarray(g_fd[:m]), rtol=1e-6, atol=1e-6)


def test_plumbing_constant_integrand(lib, hky):
    """D = Q_p makes the integrand constant: t * sum w f / L equals the edge derivative."""
    _, sc, _, _, _ = hky
    cs, bs = _scene(lib, sc, 1)
    ed, _ = bs.edge_derivatives()
    g, _ = bs.gl_gradient([[bs.trips[g] for g in range(cs.G)]])
    assert g[0] == pytest.approx(ed.sum(), rel=1e-12)


def _random_D(rng, cs, m):
    D = []
    for _ in range(m):
        Dk = []
        for g in range(cs.G):
            if rng.random() < 0.25:
                Dk.append(None)
                continue
            S = cs.S
            M = np.where(rng.random((S, S)) < 0.4, rng.normal(size=(S, S)), 0.0)
            np.fill_diagonal(M, 0.0)
            np.fill_diagonal(M, -M.sum(axis=1))
            r, c = np.nonzero(M)
            Dk.append((r.astype(np.int32), c.astype(np.int32), M[r, c]))
        D.append(Dk)
    dpi = []
    for _ in range(m):
        v = rng.normal(size=cs.S) * (rng.random() < 0.5)
        dpi.append(v - v.mean() if np.any(v) else None)
    return D, dpi


def test_random_scenes(lib):
    """GL vs block generator on random scenes with random D and dpi (at least 10 feasible scenes)."""
    done = 0
    for seed in range(30):
        sc = random_scene(seed)
        rng = np.random.default_rng(100 + seed)
        m = 3
        cs, bs = _scene(lib, sc, m, K=8)
        if not np.all(np.isfinite(bs.loglik()[1])):
            continue                                  # infeasible random scene
        D, dpi = _random_D(rng, cs, m)
        g, _ = bs.gl_gradient(D, dpi)
        ref = block_gradient(lib, bs, D, dpi)
        assert _rel(g, ref) < RTOL, seed
        done += 1
    assert done >= 10


def test_graded_panels(lib, hky):
    """Small rho0 forces graded multi-panel integration on every edge; it matches the oracle and
    beats one panel at the same K."""
    _, sc, m, D, dpi = hky
    cs, bs = _scene(lib, sc, m, K=4, K_graded=4, rho0=0.02)
    g, _ = bs.gl_gradient(D, dpi)
    panels = [len(v[1]) - 1 for v in bs._gl_nodes if v is not None]
    assert min(panels) > 2 and bs.gl_clipped == 0
    ref = block_gradient(lib, bs, D, dpi)
    cs1, bs1 = _scene(lib, sc, m, K=4, K_graded=4, rho0=1e9)
    g1, _ = bs1.gl_gradient(D, dpi)
    assert _rel(g, ref) < RTOL < _rel(g1, ref)


def test_panel_cap_is_counted(lib):
    d = tempfile.mkdtemp()
    obj = _model("hky", d, log_tau=7.0, blen_scale=10.0)
    sc = obj.get_scene()
    m = len(obj.x) - len(obj.edge_to_blen)
    cs = compile_scene(sc)
    D, dpi = _jac(obj, cs, [m - 1])
    cs, bs = _scene(lib, sc, 1, J_max=4)
    bs.gl_gradient(D, dpi)
    assert bs.gl_clipped > 0


@pytest.mark.parametrize("log_tau,scale", [(4.0, 10.0), (7.0, 10.0)])
def test_long_edges_large_tau(lib, log_tau, scale):
    d = tempfile.mkdtemp()
    obj = _model("hky", d, log_tau=log_tau, blen_scale=scale)
    sc = obj.get_scene()
    m = len(obj.x) - len(obj.edge_to_blen)
    cs = compile_scene(sc)
    D, dpi = _jac(obj, cs, list(range(m)))
    cs, bs = _scene(lib, sc, m)                               # default K=6, automatic panels
    g, _ = bs.gl_gradient(D, dpi)
    ref = block_gradient(lib, bs, D, dpi)
    assert max(len(v[1]) - 1 for v in bs._gl_nodes if v is not None) > 1     # graded panels were used
    assert bs.gl_clipped == 0
    assert _rel(g, ref) < RTOL


def test_psjs_analytic_differentials(lib):
    from IGCexpansion.beagle.psjs import PSJSDifferentials, QBuilder, scene_list
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w_igc(d)
    sl = scene_list(obj)
    qb = QBuilder(obj.psjsmodel, sorted({c for _, c in sl}))
    free = list(range(8))
    J = PSJSDifferentials(qb, [4] * 4, free)
    x = np.array(obj.x, float)
    J.at(obj.psjsmodel.x_js)
    rng = np.random.default_rng(3)
    for j in list(rng.choice(len(sl), 4, replace=False)) + [len(sl) - 1]:
        n, c = sl[j]
        with cases.quiet():
            sc = obj.get_scene(n, c)
        cs, bs = _scene(lib, sc, len(free))
        D, dpi = J.scene(c, n)
        g, _ = bs.gl_gradient(D, dpi)
        ref = block_gradient(lib, bs, D, dpi)
        assert _rel(g, ref) < RTOL

        def scene_at(xx):
            with cases.quiet():
                obj.unpack_x(xx)
                return obj.get_scene(n, c)
        D2, dpi2 = scene_fd_differentials(scene_at, x, free, cs)
        g2, _ = bs.gl_gradient(D2, dpi2)
        np.testing.assert_allclose(g2, g, rtol=1e-7, atol=1e-7)   # FD differential matrices (h = 1e-5) vs analytic


@pytest.mark.slow
def test_mg94_log_tau(lib):
    d = tempfile.mkdtemp()
    obj = _model("mg94", d)
    sc = obj.get_scene()
    npr = len(obj.x) - len(obj.edge_to_blen)
    cs = compile_scene(sc)
    D, dpi = _jac(obj, cs, [npr - 1])
    cs, bs = _scene(lib, sc, 1)
    g, _ = bs.gl_gradient(D, dpi)
    ref = block_gradient(lib, bs, D, dpi)
    assert _rel(g, ref) < RTOL
