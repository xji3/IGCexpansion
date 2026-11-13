import os
for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pytest
import scipy.linalg as sla

from IGCexpansion import psigc
from IGCexpansion.psigc import model as M
from IGCexpansion.psigc.prune import Pruner

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T = os.path.join(REPO, 'test')
YEAST_TREE = open(os.path.join(T, 'YeastTree.newick')).read()
FASTA = os.path.join(T, 'YDR418W_YEL054C_input.fasta')
SEQ_INDEX = os.path.join(T, 'YDR418W_YEL054C_seq_index.txt')
# legacy PS-IGC fit of the IGCTract study, outside this repository
IGCTRACT = os.environ.get('IGC_TRACT_DIR', os.path.join(os.path.dirname(REPO), 'IGCTract'))
LEGACY_SAVE = os.path.join(IGCTRACT, 'PSJSAnalyses', 'save', 'PSJS_HKY_YDR418W_YEL054C_One_rate_Guess_3_rv_SCOK_nonclock_save.txt')
needs_legacy_save = pytest.mark.skipif(not os.path.exists(LEGACY_SAVE), reason='legacy PS-IGC save file (IGCTract) not found')
LEGACY_NEG_LNL = 846874.7811633880483
PI = np.array([0.2869, 0.2265, 0.1894, 0.2972]); PI = PI / PI.sum()


@pytest.fixture(scope='module')
def tree():
    return psigc.Tree(YEAST_TREE)


@pytest.fixture(scope='module')
def real(tree):
    if not os.path.exists(FASTA):
        pytest.skip('real data not available')
    return psigc.Dataset.from_fasta(FASTA, SEQ_INDEX, tree, ('YDR418W', 'YEL054C'))


def synthetic(tree, L=24, seed=0, missing=()):
    rng = np.random.default_rng(seed)
    ncol = 2 * len(tree.post_tips) + len(tree.pre_tips)
    base = rng.integers(0, 4, L)
    cols = np.where(rng.random((L, ncol)) < 0.8, base[:, None], rng.integers(0, 4, (L, ncol)))
    for (i, j) in missing:
        cols[i, j] = -1
    nt = np.arange(1, L + 1)
    return psigc.Dataset(cols, nt, (nt - 1) % 3 + 1, tree)


def x_test(tree, inv_p=6.0):
    x = np.loadtxt(LEGACY_SAVE) if os.path.exists(LEGACY_SAVE) else None
    if x is None:
        par = psigc.Params(PI, 10.0, {1: 1.0, 2: 0.6, 3: 10.0}, 4.0, 1 / inv_p,
                           {e: 0.005 for e in tree.edges})
        x = psigc.pack_x(par, tree)
    x = x.copy(); x[7] = -np.log(inv_p)
    return x


# ------------------------------------------------------------------ 1. generator algebra
def test_generator_rules():
    Q1, Q3 = M.hky_Q(PI, 14.0, 1.0), M.hky_Q(PI, 14.0, 11.0)
    tau = 22.6
    A, B = M.PS256.build_AB([Q3, Q3], tau)
    A13, B13 = M.PS256.build_AB([Q1, Q3], tau)
    assert np.abs(A.sum(1)).max() < 1e-12 and np.abs(B.sum(1)).max() < 1e-12
    assert np.allclose(B, B13)                                # B is class independent
    assert int(((B != 0) & ~np.eye(256, dtype=bool)).sum()) == 864
    p, n = 1 / 30, 7
    w = (1 - p) ** n
    Qn = A + w * B
    i = M.idx256
    assert np.isclose(Qn[i(0, 1, 2, 1), i(2, 1, 2, 1)], Q3[0, 2] + tau)
    assert np.isclose(Qn[i(0, 1, 2, 3), i(2, 1, 2, 3)], Q3[0, 2] + tau * (1 - w))
    assert np.isclose(Qn[i(0, 1, 2, 3), i(2, 3, 2, 3)], tau * w)
    assert np.isclose(Qn[i(0, 1, 2, 3), i(1, 1, 2, 3)], Q3[0, 1])


def test_factorisation_and_swap():
    Q1, Q3 = M.hky_Q(PI, 14.0, 1.0), M.hky_Q(PI, 14.0, 11.0)
    tau, t = 22.6, 0.0159
    A13, B13 = M.PS256.build_AB([Q1, Q3], tau)
    perm = np.array([16 * (4 * (s // 64) + (s // 4) % 4) + (4 * ((s // 16) % 4) + s % 4) for s in range(256)])
    PK = np.kron(sla.expm(t * M.IS16.build([Q1], tau)), sla.expm(t * M.IS16.build([Q3], tau)))
    assert np.abs(sla.expm(t * A13) - PK[np.ix_(perm, perm)]).max() < 1e-13
    A31, B31 = M.PS256.build_AB([Q3, Q1], tau)
    sw = np.array([M.idx256((s // 16) % 4, s // 64, s % 4, (s // 4) % 4) for s in range(256)])
    w = 0.3
    assert np.abs((A31 + w * B31) - (A13 + w * B13)[np.ix_(sw, sw)]).max() < 1e-12


@needs_legacy_save
def test_legacy_x_roundtrip(tree):
    x = np.loadtxt(LEGACY_SAVE)
    par = psigc.unpack_x(x, tree)
    assert np.isclose(par.blen[('N0', 'kluyveri')] / par.blen[('N0', 'D1')], 9.0)
    assert np.allclose(psigc.pack_x(par, tree), x)


# ------------------------------------------------------------------ 2-3. interpolation and p = 1
def test_p_equals_one_identity(real, tree):
    x = x_test(tree); x[7] = 0.0                              # p = 1: w_n = 0 for every pair
    lnl_ps = psigc.PSEngine(real, K=16, w_tol=0.0).lnl(x)     # all pairs through the 256-state path
    lnl_is = psigc.ISEngine(real).lnl(x)
    assert abs(lnl_ps - (real.L - 1) * lnl_is) < 1e-6


def test_interpolation_vs_direct(real, tree):
    x = x_test(tree, inv_p=30.0)
    eng = psigc.PSEngine(real, K=16)
    _, lp = eng.lnl(x, pairs=True)
    par = psigc.unpack_x(x, tree)
    c = (3, 3); blk = real.blocks[c]
    Q3 = par.site_Q(3)
    A, B = M.PS256.build_AB([Q3, Q3], par.tau)
    Gpre = M.PRE16.build([Q3, Q3])
    pr = Pruner(tree, M.DIAG256)
    for n in (1, 50, 200, 480):
        sel = np.flatnonzero(blk.n == n)
        if not len(sel):
            continue
        w = (1 - par.p) ** n
        P = {e: sla.expm(par.blen[e] * (A + w * B)) for e in tree.post_edges}
        P.update({e: sla.expm(par.blen[e] * Gpre) for e in tree.pre_edges})
        L = pr.forward(P, blk.tips, np.outer(par.pi, par.pi).ravel())[0][blk.inv[sel]]
        assert np.abs(lp[eng._offsets[c] + sel] - np.log(L)).max() < 1e-10


@needs_legacy_save
def test_legacy_reproduction(real, tree):
    x = np.loadtxt(LEGACY_SAVE)
    assert abs(psigc.PSEngine(real, K=16).lnl(x) + LEGACY_NEG_LNL) < 1e-6


# ------------------------------------------------------------------ 4. missing data
def test_missing_data_marginalises(tree):
    j = 3
    ds_m = synthetic(tree, missing=[(5, j)])
    x = x_test(tree, inv_p=8.0)
    eng = psigc.PSEngine(ds_m, K=12)
    _, lp_m = eng.lnl(x, pairs=True)
    tot = 0.0
    for v in range(4):
        cols = ds_m.cols.copy(); cols[5, j] = v
        ds_v = psigc.Dataset(cols, ds_m.nt_index, ds_m.codon_pos, tree)
        tot = tot + np.exp(psigc.PSEngine(ds_v, K=12).lnl(x, pairs=True)[1])
    I, J = np.concatenate([b.pair_I for b in ds_m.blocks.values()]), np.concatenate([b.pair_J for b in ds_m.blocks.values()])
    touch = (I == 5) | (J == 5)
    assert np.abs(lp_m[touch] - np.log(tot[touch])).max() < 1e-10            # summed over the 4 completions
    assert np.abs(lp_m[~touch] - np.log(tot[~touch] / 4)).max() < 1e-10      # pairs not touching column 5


# ------------------------------------------------------------------ 5-6. exact gradient
@pytest.mark.parametrize('inv_p', [1.5, 6.0, 300.0])
def test_gradient_vs_richardson(tree, inv_p):
    ds = synthetic(tree, L=30, seed=1)
    x = x_test(tree, inv_p) + np.random.default_rng(3).normal(0, 0.1, 20)
    x[7] = -np.log(inv_p)
    eng = psigc.PSEngine(ds, K=16)
    f, g = eng.lnl_grad(x)
    assert np.isclose(f, eng.lnl(x))
    h = 1e-3
    for j in range(len(x)):
        def F(d):
            xx = x.copy(); xx[j] += d; return eng.lnl(xx)
        d1 = (F(h) - F(-h)) / (2 * h); d2 = (F(h / 2) - F(-h / 2)) / h
        fd = (4 * d2 - d1) / 3
        assert abs(g[j] - fd) <= 1e-6 * max(1.0, abs(fd)), (j, g[j], fd)


def test_tied_root_derivative_is_plain_sum(tree):
    """d lnL / d x_rate[0] must equal the FD along the tied parameter (fails on the legacy formula)."""
    ds = synthetic(tree, L=20, seed=2)
    x = x_test(tree, 6.0)
    eng = psigc.PSEngine(ds, K=12)
    _, g = eng.lnl_grad(x)
    j = 8
    h = 1e-4
    xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
    assert abs(g[j] - (eng.lnl(xp) - eng.lnl(xm)) / (2 * h)) < 1e-6 * max(1, abs(g[j]))


def test_root_by_dup_tree():
    t = psigc.Tree('((((cerevisiae,paradoxus)N5,mikatae)N4,kudriavzevii)N3,castellii)D1;')
    ds = synthetic(t, L=15, seed=4)
    par = psigc.Params(PI, 8.0, {1: 1.0, 2: 0.7, 3: 5.0}, 3.0, 0.1, {e: 0.004 for e in t.edges})
    x = psigc.pack_x(par, t)
    eng = psigc.PSEngine(ds, K=12)
    f, g = eng.lnl_grad(x)
    j, h = 7, 1e-4
    xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
    assert abs(g[j] - (eng.lnl(xp) - eng.lnl(xm)) / (2 * h)) < 1e-6 * max(1, abs(g[j]))


def test_is_gradient(tree):
    ds = synthetic(tree, L=30, seed=5)
    x = x_test(tree, 6.0) + np.random.default_rng(9).normal(0, 0.1, 20)
    eng = psigc.ISEngine(ds)
    f, g = eng.lnl_grad(x)
    h = 1e-4
    for j in range(len(x)):
        xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
        fd = (eng.lnl(xp) - eng.lnl(xm)) / (2 * h)
        assert abs(g[j] - fd) <= 1e-6 * max(1.0, abs(fd)), (j, g[j], fd)


def test_theta_transform_and_fit(tree):
    from IGCexpansion.psigc import fit
    ds = synthetic(tree, L=18, seed=6)
    x = x_test(tree, 6.0)
    t = fit.theta_from_x(x)
    assert np.allclose(fit.x_from_theta(t), x)
    eng = psigc.PSEngine(ds, K=12)
    fg = fit.theta_fg(eng)
    # psi-direction gradient via the theta chain rule vs FD
    h = 1e-4
    tp, tm = t.copy(), t.copy(); tp[7] += h; tm[7] -= h
    fd = (eng.lnl(fit.x_from_theta(tp)) - eng.lnl(fit.x_from_theta(tm))) / (2 * h)
    assert abs(fg(t)[1][7] - fd) < 1e-6 * max(1, abs(fd))
    # a 2-parameter fit (log tau, psi) increases the objective and stops with a small projected gradient
    r = fit.maximize(fg, t, free=[6, 7])
    assert r['f'] >= eng.lnl(x) - 1e-9 and r['max_pg'] < 1e-3


def test_theta_frequency_chain_rule(tree):
    from IGCexpansion.psigc import fit
    ds = synthetic(tree, L=12, seed=8)
    eng = psigc.ISEngine(ds)
    t = fit.theta_from_x(x_test(tree, 6.0)) + np.random.default_rng(1).normal(0, 0.2, 20)
    g = fit.theta_fg(eng)(t)[1]
    for j in range(3):
        tp, tm = t.copy(), t.copy(); tp[j] += 1e-5; tm[j] -= 1e-5
        fd = (eng.lnl(fit.x_from_theta(tp)) - eng.lnl(fit.x_from_theta(tm))) / 2e-5
        assert abs(g[j] - fd) < 1e-6 * max(1, abs(fd))


def test_newton_decrement_rule():
    """Decision 3: convergence is the Newton decrement. A strongly correlated quadratic, where the
    per-coordinate gradient x SE test is pessimistic, still converges, and the decrement matches the exact one."""
    from IGCexpansion.psigc import fit
    rho = 0.999
    H = np.array([[1.0, rho], [rho, 1.0]]) * 50
    mu = np.array([0.3, -0.2])
    fg = lambda t: (-0.5 * (t - mu) @ H @ (t - mu), -H @ (t - mu))
    r = fit.maximize(fg, np.zeros(2), bounds=[(None, None)] * 2)
    assert r['converged'] and r['newton_dec'] <= fit.ND_TOL
    g = r['grad']
    assert r['newton_dec'] <= 10 * np.sqrt(g @ np.linalg.solve(H, g)) + 1e-12
    assert fit.rule_converged(r) and fit.rule_converged(r, 'plan_unscaled') == r['converged_pg']
    r2 = fit.maximize(fg, np.zeros(2), bounds=[(None, None)] * 2, maxiter=1)
    assert not r2['converged'] and r2['newton_dec'] > fit.ND_TOL
