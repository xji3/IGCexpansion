"""Tests of psigc.uq (step 8) on the CPU with the GPU engine's code path."""
import os

for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pytest

torch = pytest.importorskip('torch')
from IGCexpansion import psigc                                                       # noqa: E402
from IGCexpansion.psigc import fit, uq                                          # noqa: E402
from IGCexpansion.psigc import simulate as S                                    # noqa: E402
from IGCexpansion.psigc.gpu import GPUEngine                                    # noqa: E402

YEAST = '((((((cerevisiae,paradoxus)N5,mikatae)N4,kudriavzevii)N3,bayanus)N2,castellii)N1,kluyveri)N0;'
X0 = np.array([-0.7421, -0.5077, -0.8372, 2.6722, -0.6257, 2.4455, 1.4663, -1.6503, -4.1702, -5.6241, -4.1406,
               -4.3407, -6.4641, -6.0789, -6.0463, -5.1153, -6.3955, -5.4676, -5.4562, -6.4584])
L = 45


@pytest.fixture(scope='module')
def setup():
    tree = psigc.Tree(YEAST)
    si = np.stack([np.arange(1, L + 1), np.arange(L) // 3 + 1, np.arange(L) % 3 + 1], 1)
    rng, _ = S.rng_for(11, 0, 0, 0)
    ds = S.to_dataset(S.simulate(psigc.unpack_x(X0, tree), tree, si, rng), tree, si)
    eng = GPUEngine(ds, device='cpu')
    return tree, si, ds, eng


def test_site_split_sums_to_score(setup):
    _, _, ds, eng = setup
    t = fit.theta_from_x(X0)
    v, u, U = uq.site_split(eng, X0, t)
    g = fit.theta_fg(eng)(t)[1]
    assert np.abs(U - g).max() <= 1e-9 * max(1, np.abs(g).max())
    assert np.abs(v.sum(0) - U).max() <= 1e-9 * max(1, np.abs(U).max())
    _, Ss = eng.site_scores(X0)
    assert np.abs(uq.theta_rows(Ss, t)[:, fit.I_PSI]).max() <= 1e-10      # no psi in the IS site score
    # far pairs factor exactly: d_ab = 0 when w_n underflows
    x = X0.copy(); x[7] = np.log(0.9)                                    # p = 0.9: w_n < 1e-15 beyond n = 15
    t2 = fit.theta_from_x(x)
    _, u2, _ = uq.site_split(eng, x, t2)
    s2 = uq.theta_rows(eng.site_scores(x)[1], t2)
    I, J = eng.pair_sites()
    far = ds.nt_index[J] - ds.nt_index[I] > 20
    assert np.abs(u2[far] - s2[I[far]] - s2[J[far]]).max() <= 1e-10


def test_hac_limits(setup):
    _, _, ds, eng = setup
    t = fit.theta_from_x(X0)
    v, u, U = uq.site_split(eng, X0, t)
    J, info = uq.j_hac(v, ds.nt_index, p_hat=1e-9)                       # B >> L: weights ~ 1 -> J ~ U U'
    assert not info['valid']
    w_min = 1 - L / (info['bandwidth'] + 1)
    assert np.allclose(J, np.outer(U, U), rtol=0, atol=(1 - w_min) * np.abs(v).sum(0).max() ** 2 + 1e-9)
    J2, info2 = uq.j_hac(v, ds.nt_index, p_hat=0.5)
    assert np.linalg.eigvalsh(J2).min() >= -1e-8 * np.abs(J2).max()        # PSD
    assert np.allclose(uq.j_legacy(u), u.T @ u)


def test_godambe_and_summary(setup):
    _, _, ds, eng = setup
    t = fit.theta_from_x(X0)
    fg = fit.theta_fg(eng)
    free = uq.free_mask(t)
    H = uq.hessian(fg, t, free)
    V, Hinv = uq.godambe(H, H, free)
    f = np.flatnonzero(free)
    assert np.allclose(V[np.ix_(f, f)], Hinv[np.ix_(f, f)], rtol=1e-8, atol=1e-12)
    s = uq.summarize(t, V, Hinv)
    assert abs(s['lambda_hat'] - 1) < 1e-8
    assert s['psi']['lo'] < t[fit.I_PSI] < s['psi']['hi']
    # a zero-length branch is held fixed; psi on its lower bound stays free (lambda_hat for p = 1 fits)
    # unless keep_psi=False
    t2 = t.copy(); t2[fit.I_PSI] = fit.PSI_MIN; t2[10] = -25
    m = uq.free_mask(t2)
    assert m[fit.I_PSI] and not m[10] and m[fit.I_LOGTAU]
    assert not uq.free_mask(t2, keep_psi=False)[fit.I_PSI]


def test_j_pb_and_driver(setup):
    tree, si, ds, eng = setup
    t = fit.theta_from_x(X0)
    out = uq.uq_dataset(eng, t, YEAST, si, R=4, entropy=5, key=(9, 0),
                        make_engine=lambda d: GPUEngine(d, device='cpu', graphs=False), workers=2)
    assert out['PB']['psi']['se'] > 0 and out['HAC']['psi']['se'] > 0 and out['legacy']['psi']['se'] > 0
    assert len(out['matrices']['J_PB']) == len(t)
    # same key -> same replicates
    J1, U1, _ = uq.j_pb(lambda d: GPUEngine(d, device='cpu', graphs=False), t, YEAST, si, 2, 5, (9, 0), workers=2)
    assert np.allclose(U1, uq.j_pb(lambda d: GPUEngine(d, device='cpu', graphs=False), t, YEAST, si, 2, 5,
                                   (9, 0), workers=2)[1])


def _quad(n=20, mu_psi=3.0, seed=0):
    rng = np.random.default_rng(seed)
    M = rng.standard_normal((n, n))
    A = M @ M.T / n + np.eye(n) * 0.5
    A *= 4.0
    mu = rng.standard_normal(n) * 0.3; mu[fit.I_PSI] = mu_psi
    return (lambda t: (-0.5 * (t - mu) @ A @ (t - mu), -A @ (t - mu))), A, mu


def test_profile_clr_quadratic():
    from IGCexpansion.psigc import clr
    fg, A, mu = _quad()
    r = fit.maximize(fg, np.zeros(20))
    t_hat, cl_hat = r['theta'], r['f']
    prof = clr.m3_profile(fg, t_hat)
    s2 = np.linalg.inv(A)[fit.I_PSI, fit.I_PSI]
    W = clr.W_values(cl_hat, prof)
    ps = np.array([p['psi'] for p in prof])
    assert np.allclose(W, (ps - mu[fit.I_PSI]) ** 2 / s2, atol=1e-3)
    for lam in (1.0, 2.0):
        ci = clr.clr_interval(fg, t_hat, cl_hat, lam, prof)
        h = np.sqrt(clr.CUT95 * lam * s2)
        assert abs(ci['lo'] - (mu[fit.I_PSI] - h)) < 0.03 and abs(ci['hi'] - (mu[fit.I_PSI] + h)) < 0.03
    bt = clr.boundary_test(cl_hat, prof, 1.0)
    assert abs(bt['W0'] - (fit.PSI_MIN - mu[fit.I_PSI]) ** 2 / s2) < 1e-3 and bt['reject_5pct']
    assert not clr.plateau(t_hat, prof)


def test_clr_open_upper_and_plateau():
    from IGCexpansion.psigc import clr
    fg, A, mu = _quad(mu_psi=20.0)                               # optimum beyond the psi bound
    r = fit.maximize(fg, np.zeros(20))
    assert r['theta'][fit.I_PSI] >= fit.PSI_MAX - 1e-6
    prof = clr.m3_profile(fg, r['theta'])
    assert clr.plateau(r['theta'], prof)
    ci = clr.clr_interval(fg, r['theta'], r['f'], 1.0, prof)
    assert ci['hi'] == np.inf and ci['lo'] > 0


def test_infer_dataset(setup):
    from IGCexpansion.psigc import infer
    tree, si, ds, eng = setup
    t = fit.theta_from_x(X0)
    r = fit.maximize(fit.theta_fg(eng), t)
    out = infer.infer_dataset(eng, r['theta'], r['f'], float(np.exp(t[fit.I_LOGTAU])), YEAST, si, R=3, entropy=1,
                              key=(4,), make_engine=lambda d: GPUEngine(d, device='cpu', graphs=False), workers=2)
    assert out['clr']['one'] is not None and out['boundary']['one'] is not None
    assert out['lambda']['PB'] > 0 and len(out['profile']) == 11


def test_chunking_multilocus_noncoding(setup):
    tree, si, ds, eng = setup
    # pattern chunking: same lnL / gradient; per-pair values equal after matching (I, J)
    ec = GPUEngine(ds, device='cpu', max_patterns=40)
    assert len(ec.classes) > len(eng.classes)
    f0, g0 = eng.lnl_grad(X0); f1, g1 = ec.lnl_grad(X0)
    assert abs(f1 - f0) <= 1e-9 and np.abs(g1 - g0).max() <= 1e-9 * max(1, np.abs(g0).max())
    p0 = eng.lnl(X0, pairs=True)[1]; p1 = ec.lnl(X0, pairs=True)[1]
    I0, J0 = eng.pair_sites(); I1, J1 = ec.pair_sites()
    o0, o1 = np.lexsort((J0, I0)), np.lexsort((J1, I1))
    assert np.abs(p0[o0] - p1[o1]).max() <= 1e-12
    # two loci: lnL is the sum over loci; v sums to U
    rng, _ = S.rng_for(11, 0, 0, 1)
    d2 = S.to_dataset(S.simulate(psigc.unpack_x(X0, tree), tree, si[:30], rng), tree, si[:30])
    dm = psigc.Dataset.concat([ds, d2])
    assert dm.n_pairs == ds.n_pairs + d2.n_pairs and dm.n_loci == 2
    em = GPUEngine(dm, device='cpu')
    fm, gm = em.lnl_grad(X0)
    f2, g2 = GPUEngine(d2, device='cpu').lnl_grad(X0)
    assert abs(fm - (f0 + f2)) <= 1e-9 and np.abs(gm - (g0 + g2)).max() <= 1e-9 * max(1, np.abs(gm).max())
    t = fit.theta_from_x(X0)
    v, u, U = uq.site_split(em, X0, t)
    assert np.abs(v.sum(0) - U).max() <= 1e-9 * max(1, np.abs(U).max())
    Jc, K = uq.j_cr(v, dm.locus)
    assert K == 2 and np.linalg.matrix_rank(Jc) <= 2
    J, info = uq.j_hac(v, dm.nt_index, 0.5, locus=dm.locus)
    assert np.linalg.eigvalsh(J).min() >= -1e-8 * np.abs(J).max()
    # non-coding: codon position 0 -> no dependence on r2, r3
    dn = psigc.Dataset(ds.cols, ds.nt_index, np.zeros(ds.L, int), tree)
    assert dn.noncoding and list(dn.blocks) == [(1, 1)]
    gn = GPUEngine(dn, device='cpu').lnl_grad(X0)[1]
    assert gn[4] == 0 and gn[5] == 0
    assert uq.fixed_coords(dn) == [4, 5] and not uq.free_mask(t, fixed=[4, 5])[4]
    # multi-locus simulation for J_PB
    dsim, st = uq._sim_one((X0, YEAST, [si, si[:30]], 3, (1, 2)))
    assert dsim.n_loci == 2 and dsim.L == L + 30


def test_boundary_psi_kept_free_and_forward_hessian():
    from IGCexpansion.psigc import uq, fit
    n = 12
    t = np.zeros(n); t[fit.I_PSI] = fit.PSI_MIN
    free = uq.free_mask(t)
    assert free[fit.I_PSI]
    assert not uq.free_mask(t, keep_psi=False)[fit.I_PSI]
    A = np.diag(np.arange(1, n + 1.0)); A[0, fit.I_PSI] = A[fit.I_PSI, 0] = 0.3

    def fg(th):                                   # cl = -1/2 th' A th (+ cubic in psi to test the forward step)
        g = -A @ th; g[fit.I_PSI] -= 0.5 * th[fit.I_PSI] ** 2
        return -0.5 * th @ A @ th - th[fit.I_PSI] ** 3 / 6, g
    H = uq.hessian(fg, t, free)
    assert np.allclose(H, A, atol=2e-4)


def test_jvp_columns_retries_on_oom():
    import torch
    from IGCexpansion.psigc import gpu
    calls = []
    real = gpu._vmap_jvp

    def flaky(fn, x, tans):
        calls.append(tans.shape[0])
        if tans.shape[0] > 2:
            raise torch.OutOfMemoryError('CUDA out of memory (test)')
        return real(fn, x, tans)
    gpu._vmap_jvp = flaky
    try:
        x = torch.arange(1.0, 8.0, dtype=gpu.DT)
        prim, S = gpu._jvp_columns(lambda v: torch.stack([v.sum(), (v ** 2).sum()]), x, 7, 5, torch.device('cpu'))
    finally:
        gpu._vmap_jvp = real
    assert calls[0] == 5 and max(calls[1:]) <= 2
    assert torch.allclose(S, torch.stack([torch.ones(7, dtype=gpu.DT), 2 * x]))


def test_expected_hessian_pb(setup):
    """H_E = mean_r sum_ab u_ab u_ab' (per-pair Bartlett identity): PSD, scores identical to lnl_grad,
    same scale as the observed Hessian, and a positive lambda_hat."""
    tree, si, ds, eng = setup
    t = fit.theta_from_x(X0)
    mk = lambda d: GPUEngine(d, device='cpu', graphs=False)
    J1, U1, i1 = uq.j_pb(mk, t, YEAST, si, 12, 5, (9, 3), workers=2)
    J2, U2, i2 = uq.j_pb(mk, t, YEAST, si, 12, 5, (9, 3), workers=2, expected_h=True)
    assert 'H_E' not in i1 and i2['H_E_R'] == 12
    assert np.abs(U1 - U2).max() <= 1e-8 * max(1, np.abs(U1).max())    # pair-score sum = gradient
    HE = i2['H_E']
    free = uq.free_mask(t, fixed=uq.fixed_coords(ds))
    f = np.flatnonzero(free)
    assert np.linalg.eigvalsh(HE[np.ix_(f, f)]).min() >= -1e-8 * np.abs(HE).max()
    H = uq.hessian(fit.theta_fg(eng), t, free)
    for i in (fit.I_PSI, fit.I_LOGTAU, 0, 3):
        assert 0.2 < HE[i, i] / H[i, i] < 5, (i, HE[i, i], H[i, i])
    out = uq.uq_dataset(eng, t, YEAST, si, R=12, entropy=5, key=(9, 3), make_engine=mk, workers=2,
                        expected_h=True)
    assert out['PB_HE']['lambda_hat'] > 0 and out['H_E']['R'] == 12 and 'H_E' in out['matrices']


def test_clr_interval_local_quadratic():
    from scipy.optimize import brentq
    from IGCexpansion.psigc import clr
    fg, A, mu = _quad()
    r = fit.maximize(fg, np.zeros(20))
    t_hat, cl_hat = r['theta'], r['f']
    prof = clr.m3_profile(fg, t_hat)
    s2 = np.linalg.inv(A)[fit.I_PSI, fit.I_PSI]
    m = mu[fit.I_PSI]
    # constant lambda: same as clr_interval
    loc = clr.clr_interval_local(fg, t_hat, cl_hat, prof, lambda i: 2.0, 2.0)
    ref = clr.clr_interval(fg, t_hat, cl_hat, 2.0, prof)
    assert abs(loc['lo'] - ref['lo']) < 0.03 and abs(loc['hi'] - ref['hi']) < 0.03
    # lambda varying with psi (log-linear between grid points is exact for exp-linear lambda)
    lam = lambda psi: np.exp(0.3 * psi)
    loc = clr.clr_interval_local(fg, t_hat, cl_hat, prof, lambda i: float(lam(prof[i]['psi'])),
                                 float(lam(t_hat[fit.I_PSI])))
    f = lambda p: (p - m) ** 2 / s2 - clr.CUT95 * lam(p)
    hi = brentq(f, m, 12.0)
    lo = brentq(f, -5, m) if f(fit.PSI_MIN) > 0 else 0.0
    assert abs(loc['hi'] - hi) < 0.05 and abs(loc['lo'] - max(lo, 0.0)) < 0.05, (loc, lo, hi)
    assert len(loc['probes']) < 11                                  # bisection, not every grid point
    # None (no PD H_E) counts as inside
    loc = clr.clr_interval_local(fg, t_hat, cl_hat, prof, lambda i: None, 1.0)
    assert loc['lo'] == 0.0 and loc['hi'] == np.inf
