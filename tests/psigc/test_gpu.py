"""Tests of psigc.gpu.GPUEngine against the CPU engine (runs on the CPU; CUDA-only tests skip without a GPU)."""
import os
import pathlib

for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pytest
import scipy.linalg as sla

torch = pytest.importorskip('torch')
from IGCexpansion import psigc                                                       # noqa: E402
from IGCexpansion.psigc import simulate as S                                    # noqa: E402
from IGCexpansion.psigc import gpu as G                                         # noqa: E402

YEAST = '((((((cerevisiae,paradoxus)N5,mikatae)N4,kudriavzevii)N3,bayanus)N2,castellii)N1,kluyveri)N0;'
X0 = np.array([-0.7421, -0.5077, -0.8372, 2.6722, -0.6257, 2.4455, 1.4663, -1.6503, -4.1702, -5.6241, -4.1406,
               -4.3407, -6.4641, -6.0789, -6.0463, -5.1153, -6.3955, -5.4676, -5.4562, -6.4584])   # legacy fit, rounded
DEV = 'cpu'


@pytest.fixture(scope='module')
def tree():
    return psigc.Tree(YEAST)


@pytest.fixture(scope='module')
def x0(tree):
    assert len(X0) == 8 + tree.n_rate
    return X0


@pytest.fixture(scope='module')
def ds(tree, x0):
    L = 45
    si = np.stack([np.arange(1, L + 1), np.arange(L) // 3 + 1, np.arange(L) % 3 + 1], 1)
    rng, _ = S.rng_for(11, 0, 0, 0)
    res = S.simulate(psigc.unpack_x(x0, tree), tree, si, rng)
    d = S.to_dataset(res, tree, si)
    cols = d.cols.copy()
    cols[::9, 2] = -1                                              # missing data -> mask tips
    return psigc.Dataset(cols, si[:, 0], si[:, 2], tree)


def points(x0):
    xs = {'x0': x0}
    for name, ip in (('p300', 300.0), ('p1.5', 1.5)):
        x = x0.copy(); x[7] = -np.log(ip); xs[name] = x
    return xs


@pytest.mark.parametrize('det', [False, True])
def test_matches_cpu(ds, x0, det):
    cpu = psigc.PSEngine(ds, K=16)
    eng = G.GPUEngine(ds, K=16, device=DEV, deterministic=det)
    for name, x in points(x0).items():
        fc, gc = cpu.lnl_grad(x)
        fg, gg = eng.lnl_grad(x)
        assert abs(fg - fc) <= 1e-9, name
        assert np.abs(gg - gc).max() <= 1e-10 * max(1.0, np.abs(gc).max()), name
        assert abs(eng.lnl(x) - fc) <= 1e-9
        _, pc = cpu.lnl(x, pairs=True)
        _, pg = eng.lnl(x, pairs=True)
        assert np.abs(pg - pc).max() <= 1e-12


def test_auto_k_matches_direct_and_ratchet(ds, x0):
    eng = G.GPUEngine(ds, device=DEV)
    x = x0.copy(); x[6] += np.log(16.0)
    eng.lnl_grad(x)
    assert eng.K_floor >= 16
    _, pp = eng.lnl(x, pairs=True)
    assert np.abs(pp - eng.pair_lnl_direct(x)).max() <= 1e-10
    k = eng.K_floor
    eng.lnl(x0)
    assert eng.K_floor == k                                        # ratchet: does not fall within a task
    eng.begin_task()
    assert eng.K_floor == G.K_LADDER[0]


def test_exact_fallback(ds, x0, monkeypatch):
    monkeypatch.setattr(G, 'ERR_CAP', 1e-300)
    eng = G.GPUEngine(ds, device=DEV, tail_tol=1e-300)
    f, g = eng.lnl_grad(x0)
    assert len(eng.last['fallback']) == len(eng.classes)
    pd = eng.pair_lnl_direct(x0)
    assert abs(f - pd.sum()) <= 1e-9
    fc, gc = psigc.PSEngine(ds, K=16).lnl_grad(x0)
    assert np.abs(g - gc).max() <= 1e-9 * max(1.0, np.abs(gc).max())
    _, pp = eng.lnl(x0, pairs=True)
    assert np.abs(pp - pd).max() <= 1e-12


def test_cpu_engine_auto_and_fallback(ds, x0, monkeypatch):
    x = x0.copy(); x[6] += np.log(16.0)
    ref = G.GPUEngine(ds, device=DEV).pair_lnl_direct(x).sum()
    eng = psigc.PSEngine(ds, K='auto')
    f, g = eng.lnl_grad(x)
    assert abs(f - ref) <= 1e-8
    # forced exact path: gradient against central differences of the exact lnL
    import IGCexpansion.psigc.cheb as C
    monkeypatch.setattr(C, 'tail_min', lambda F, **k: np.inf)
    e2 = psigc.PSEngine(ds, K='auto')
    f2, g2 = e2.lnl_grad(x)
    assert e2.last['fallback'] and abs(f2 - ref) <= 1e-8
    h = 1e-5
    for j in (6, 7, 9):
        xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
        fd = (e2.lnl(xp) - e2.lnl(xm)) / (2 * h)
        assert abs(fd - g2[j]) <= 1e-5 * max(1.0, abs(g2[j]))


def test_dup_root_tree():
    tr = psigc.Tree('((a,b)N1,(c,d)N2)D1;')
    rng = np.random.default_rng(3)
    cols = rng.integers(0, 4, size=(30, 8)); cols[rng.random(cols.shape) < 0.05] = -1
    d = psigc.Dataset(cols, np.arange(1, 31), (np.arange(30) % 3) + 1, tr)
    x = np.r_[X0[:8], np.log([0.1, 0.2, 0.05, 0.3, 0.15, 0.25])]
    fc, gc = psigc.PSEngine(d, K=16).lnl_grad(x)
    fg, gg = G.GPUEngine(d, K=16, device=DEV).lnl_grad(x)
    assert abs(fg - fc) <= 1e-9 and np.abs(gg - gc).max() <= 1e-9


def test_pair_scores(ds, x0):
    eng = G.GPUEngine(ds, device=DEV)
    lp, Sc = eng.pair_scores(x0)
    f, g = eng.lnl_grad(x0)
    assert abs(lp.sum() - f) <= 1e-9
    assert np.abs(Sc.sum(0) - g).max() <= 1e-10 * max(1.0, np.abs(g).max())
    j, h = 7, 1e-5
    xp, xm = x0.copy(), x0.copy(); xp[j] += h; xm[j] -= h
    fd = (eng.lnl(xp, pairs=True)[1] - eng.lnl(xm, pairs=True)[1]) / (2 * h)
    assert np.abs(fd - Sc[:, j]).max() <= 1e-6


def test_expm_unif_kernel():
    """V0 (small): uniformization vs SciPy, non-negativity, adjoint vs expm_frechet, exact d/dt."""
    rng = np.random.default_rng(0)
    for n, scale, s in ((16, 0.5, 0), (16, 3.0, 2), (64, 1.0, 0)):
        Q = rng.random((3, n, n)) * scale
        np.fill_diagonal(Q[0], 0); np.fill_diagonal(Q[1], 0); np.fill_diagonal(Q[2], 0)
        Q -= np.einsum('gij,ij->gi', Q, np.ones((n, n)))[:, :, None] * np.eye(n)
        t = np.array([0.3, 1.1])
        mu = np.abs(np.einsum('gii->gi', Q)).max() * t.max()
        m, s_ = G.unif_plan(mu)
        Qt = torch.tensor(Q, requires_grad=True); tt = torch.tensor(t, requires_grad=True)
        P = G.expm_unif(Qt, tt, m, s_)
        ref = np.stack([[sla.expm(tt_ * Qg) for tt_ in t] for Qg in Q])
        assert np.abs(P.detach().numpy() - ref).max() <= 1e-14
        assert (P.detach().numpy() >= 0).all()
        Gb = rng.standard_normal(P.shape)
        (P * torch.tensor(Gb)).sum().backward()
        for gi in range(3):
            adj = sum(t[e] * sla.expm_frechet(t[e] * Q[gi].T, Gb[gi, e], compute_expm=False) for e in range(2))
            assert np.abs(Qt.grad[gi].numpy() - adj).max() <= 1e-14 * np.abs(Gb).max() * 10
        dt = [sum(np.sum(Gb[gi, e] * (Q[gi] @ ref[gi, e])) for gi in range(3)) for e in range(2)]
        assert np.allclose(tt.grad.numpy(), dt, rtol=1e-13, atol=1e-13)


def test_no_matrix_exp_in_engine():
    src = (pathlib.Path(G.__file__)).read_text()
    assert 'linalg.matrix_exp' not in src.replace("'linalg.matrix_exp'", '')


@pytest.mark.skipif(not torch.cuda.is_available(), reason='needs CUDA')
def test_graph_capture(ds, x0):
    """V2: CUDA-graph replay equals eager evaluation (lnL bitwise, gradient to 1e-12); captures are reused."""
    e1 = G.GPUEngine(ds, K=16, device='cuda:0', graphs=True, graph_mem_gb=None)   # no budget clears here
    e2 = G.GPUEngine(ds, K=16, device='cuda:0', graphs=False)
    for x in points(x0).values():
        f1, g1 = e1.lnl_grad(x)
        f2, g2 = e2.lnl_grad(x)
        assert f1 == f2 and e1.lnl(x) == e2.lnl(x)                  # lnL bitwise (plan V2)
        assert np.abs(g1 - g2).max() <= 1e-12 * np.abs(g2).max()     # cuBLAS may pick other kernels in capture
        assert np.array_equal(e1.lnl_grad(x)[1], g1)                 # replays are bitwise repeatable
    n = e1.n_capture
    e1.lnl_grad(x0)
    assert e1.n_capture == n
