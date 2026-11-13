"""Uncertainty quantification for the M2 fit (SIMULATION_REDO_PLAN.md section 6; PSIGC_GPU_PLAN.md step 8).

All quantities are in theta coordinates (fit.theta_from_x): (3 frequency params, log kappa, log r2, log r3,
log tau, psi = log 1/p, edge log-rates). With cl the composite log-likelihood:

* U = sum_{a<b} grad l_ab (the score); H = -Hessian of cl by central differences of the exact gradient.
* V = H^-1 J H^-1 (no division by the number of pairs).
* J estimators: J_PB (primary: covariance of U(theta_hat; X*) over R datasets simulated at theta_hat, scores
  only), J_HAC (robustness: Bartlett-weighted sum over site contributions v_a), and the legacy Godambe J
  (sum of outer products of pair scores, i.e. pairs treated as independent; a negative control).
* lambda_hat = V_psipsi / (H^-1)_psipsi, the CLR adjustment; Wald intervals for psi, log tau, log eta.

Coordinates that sit on a bound (psi at PSI_MIN or PSI_MAX) or whose edge log-rate has drifted below
ZERO_BRANCH (a zero-length branch: the likelihood is flat there and the Hessian entry is ~1e-6) are held
fixed: H, J and V are computed on the remaining free coordinates; fixed rows of V are NaN.
"""
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np

from . import fit

ZERO_BRANCH = -15.0          # edge log-rate below this = zero-length branch, held fixed in V
I_PSI, I_LOGTAU = fit.I_PSI, fit.I_LOGTAU
Z95 = 1.959963984540054


# ------------------------------------------------------------------ coordinates
def theta_rows(Gx, t):
    """Row-wise chain rule from legacy-x gradients (n, nx) to theta gradients (n, nx)."""
    Gx = np.atleast_2d(np.asarray(Gx, float))
    G = Gx.copy()
    G[:, :3] = Gx[:, :3] @ fit._t3_jacobian(np.asarray(t[:3], float))
    G[:, I_LOGTAU] = Gx[:, 6]
    G[:, I_PSI] = -Gx[:, 6] - Gx[:, 7]
    return G


def fixed_coords(ds):
    """Coordinates the data cannot inform: r2, r3 (theta 4, 5) when every site is non-coding."""
    return [4, 5] if getattr(ds, 'noncoding', False) else []


def free_mask(t, bounds=None, zero_branch=ZERO_BRANCH, fixed=(), keep_psi=True):
    """Coordinates of the sandwich. With keep_psi, psi stays free on its lower bound (p = 1 fits): lambda_hat
    is still needed there for the upper CLR limit; its Hessian column then uses a forward difference."""
    t = np.asarray(t, float)
    bounds = fit.default_bounds(len(t)) if bounds is None else bounds
    free = np.ones(len(t), bool)
    free[list(fixed)] = False
    for i, (lo, hi) in enumerate(bounds):
        at_lo = lo is not None and t[i] <= lo + 1e-9
        if keep_psi and i == I_PSI and at_lo and not (hi is not None and t[i] >= hi - 1e-9):
            continue
        if at_lo or (hi is not None and t[i] >= hi - 1e-9):
            free[i] = False
    free[8:] &= t[8:] > zero_branch
    return free


def hessian(fg, t, free=None, h=1e-4, bounds=None):
    """H = -d2 cl / dtheta2 on the free coordinates (full-size matrix, zeros elsewhere). Central differences
    of the exact gradient; a coordinate within h of its lower bound gets a forward difference instead."""
    t = np.asarray(t, float)
    n = len(t)
    bounds = fit.default_bounds(n) if bounds is None else bounds
    idx = np.flatnonzero(np.ones(n, bool) if free is None else free)
    near = [j for j in idx if bounds[j][0] is not None and t[j] - h < bounds[j][0]]
    H = np.zeros((n, n))
    g0 = fg(t)[1] if near else None
    for j in idx:
        tp = t.copy(); tp[j] += h
        if j in near:
            H[:, j] = (fg(tp)[1] - g0) / h
        else:
            tm = t.copy(); tm[j] -= h
            H[:, j] = (fg(tp)[1] - fg(tm)[1]) / (2 * h)
    return -(H + H.T) / 2


# ------------------------------------------------------------------ sandwich
def godambe(H, J, free):
    """V = H^-1 J H^-1 on the free block; NaN on fixed rows/columns. Also returns H^-1 (same layout)."""
    n = H.shape[0]
    f = np.flatnonzero(free)
    Hf = H[np.ix_(f, f)]
    Hi = np.linalg.inv(Hf)
    Hi = (Hi + Hi.T) / 2
    Vf = Hi @ J[np.ix_(f, f)] @ Hi
    V = np.full((n, n), np.nan); V[np.ix_(f, f)] = (Vf + Vf.T) / 2
    Hinv = np.full((n, n), np.nan); Hinv[np.ix_(f, f)] = Hi
    return V, Hinv


def contrasts(n):
    """Linear functionals of theta reported with Wald intervals."""
    c_psi = np.zeros(n); c_psi[I_PSI] = 1
    c_tau = np.zeros(n); c_tau[I_LOGTAU] = 1
    return {'psi': c_psi, 'log_tau': c_tau, 'log_eta': c_tau - c_psi}


def summarize(t, V, Hinv):
    """SEs and 95% Wald intervals for psi, log tau, log eta, and lambda_hat."""
    t = np.asarray(t, float)
    out = {}
    for name, c in contrasts(len(t)).items():
        m = c != 0
        if np.isnan(V[np.ix_(m, m)]).any():
            out[name] = {'est': float(c @ t), 'se': None, 'lo': None, 'hi': None}
            continue
        cc = np.where(m, c, 0.0)
        Vz = np.nan_to_num(V)
        se = float(np.sqrt(max(cc @ Vz @ cc, 0.0)))
        est = float(cc @ t)
        out[name] = {'est': est, 'se': se, 'lo': est - Z95 * se, 'hi': est + Z95 * se}
    hp = Hinv[I_PSI, I_PSI]
    out['lambda_hat'] = None if not np.isfinite(hp) else float(V[I_PSI, I_PSI] / hp)
    out['se_naive_psi'] = None if not np.isfinite(hp) else float(np.sqrt(hp))
    return out


# ------------------------------------------------------------------ J estimators
def site_split(eng, x, t):
    """Exact split of the score: u_ab = s_a + s_b + d_ab, with s_a the single-site IS score. Returns
    v (L, n) with v_a = c_a s_a + 1/2 sum_b d_ab (so sum_a v_a = U), c_a the number of pairs of site a
    (L - 1 for one locus, L_locus - 1 with several), the pair scores u (P, n) and U."""
    _, Sp = eng.pair_scores(x)
    _, Ss = eng.site_scores(x)
    u, s = theta_rows(Sp, t), theta_rows(Ss, t)
    I, J = eng.pair_sites()
    c = eng.ds.pair_count_per_site()
    d = u - s[I] - s[J]
    half = np.zeros_like(s)
    np.add.at(half, I, 0.5 * d)
    np.add.at(half, J, 0.5 * d)
    v = c[:, None] * s + half
    return v, u, u.sum(0)


def andrews_bartlett(v):
    """Andrews (1991) AR(1) plug-in bandwidth for the Bartlett kernel, unit weights over the components
    with non-zero variance; v is (T, n) in site order."""
    T = v.shape[0]
    num = den = 0.0
    for c in range(v.shape[1]):
        y = v[:, c]
        if np.var(y) <= 1e-300:
            continue
        rho = float(np.dot(y[1:], y[:-1]) / np.dot(y[:-1], y[:-1]))
        rho = min(max(rho, -0.97), 0.97)
        s2 = float(np.mean((y[1:] - rho * y[:-1]) ** 2))
        num += 4 * rho ** 2 * s2 ** 2 / ((1 - rho) ** 6 * (1 + rho) ** 2)
        den += s2 ** 2 / (1 - rho) ** 4
    alpha = num / den if den > 0 else 0.0
    return 1.1447 * (alpha * T) ** (1 / 3)


def j_hac(v, pos, p_hat, L=None, locus=None):
    """Bartlett HAC estimate of J = Var(sum_a v_a) with weights max(0, 1 - |pos_a - pos_b| / (B + 1)) over
    nucleotide distance (a positive-definite function, so J is PSD for any spacing). Bandwidth
    B = max(2 n*, Andrews), n* = log 1e-3 / log(1 - p_hat), the distance at which w_n = 1e-3. Valid only
    when B <= L / 4 (plan section 6). With several loci the weights are zero across loci (independent
    loci), the Andrews bandwidth is computed on the concatenated loci, and L is the median locus length."""
    pos = np.asarray(pos, float)
    locus = np.zeros(len(pos), np.int64) if locus is None else np.asarray(locus)
    if L is None:
        L = float(np.median(np.bincount(locus)[np.unique(locus)]))
    order = np.lexsort((pos, locus))
    v, pos, locus = v[order], pos[order], locus[order]
    nstar = np.log(1e-3) / np.log1p(-p_hat) if p_hat < 1 else 0.0
    b_and = andrews_bartlett(v)
    B = max(2 * nstar, b_and)
    D = np.abs(pos[:, None] - pos[None, :])
    W = np.clip(1 - D / (B + 1), 0, None) * (locus[:, None] == locus[None, :])
    J = v.T @ W @ v
    return (J + J.T) / 2, {'bandwidth': float(B), 'n_star': float(nstar), 'andrews': float(b_and),
                           'valid': bool(B <= L / 4)}


def j_cr(v, locus):
    """Cluster J over loci: sum_l U_l U_l' with U_l = sum of v_a over locus l (rank <= K - 1 at theta_hat;
    use for scalar quantities only, plan section 6)."""
    K = int(np.max(locus)) + 1
    Ul = np.zeros((K, v.shape[1]))
    np.add.at(Ul, np.asarray(locus), v)
    return Ul.T @ Ul * K / max(K - 1, 1), K


def j_legacy(u):
    """Legacy Godambe J: pair scores treated as independent (negative control)."""
    return u.T @ u


# ------------------------------------------------------------------ J_PB
def _sim_one(args):
    """One simulated dataset at x_hat. seq_index: one (L, 3) array, or a list of them (one per locus, each
    simulated independently with key + (locus,), then concatenated)."""
    x_hat, newick, seq_index, entropy, key = args
    from . import simulate as S
    from .data import Dataset
    from .model import unpack_x
    from .tree import Tree
    tree = Tree(newick)
    par = unpack_x(np.asarray(x_hat), tree)
    if isinstance(seq_index, (list, tuple)):
        dss, N = [], 0
        for l, si in enumerate(seq_index):
            rng, _ = S.rng_for(entropy, *key, l)
            res = S.simulate(par, tree, si, rng)
            dss.append(S.to_dataset(res, tree, si)); N += res.summary()['N']
        return Dataset.concat(dss), {'N': N}
    rng, _ = S.rng_for(entropy, *key)
    res = S.simulate(par, tree, seq_index, rng)
    return S.to_dataset(res, tree, seq_index), res.summary()


def j_pb(make_engine, t_hat, newick, seq_index, R, entropy, key, workers=8, pool=None, expected_h=False, h_reps=None):
    """J_PB: covariance of U(theta_hat; X*) over R datasets simulated at theta_hat. Datasets are simulated
    and built by a CPU process pool (spawn) while the engine evaluates the previous ones; key + (r,) gives
    replicate r its own PCG64 stream. make_engine(ds) -> engine with lnl_grad (e.g. GPUEngine(ds,
    graphs=False): a new dataset would need new graph captures). Returns J, the score matrix (R, n) and timings.

    expected_h: also estimate the expected Hessian H_E = E[-d2 cl] = sum_ab E[u_ab u_ab'] (second Bartlett
    identity, which holds pair by pair because each pair term is a genuine likelihood of that pair's data),
    as the mean over replicates of sum_ab u_ab u_ab' from the per-pair scores. Positive semi-definite by
    construction. Replicates where pair_scores needs the exact fallback are skipped for H_E (counted).
    Returned in the info dict as 'H_E' (with 'H_E_R', the replicates used). h_reps: per-pair scores (and
    H_E) only on the first h_reps replicates, plain gradients on the rest; each H_E term already sums over
    all pairs, so H_E needs far fewer replicates than J (pair scores dominate the cost)."""
    t_hat = np.asarray(t_hat, float)
    x_hat = fit.x_from_theta(t_hat)
    si = list(map(np.asarray, seq_index)) if isinstance(seq_index, (list, tuple)) else np.asarray(seq_index)
    jobs = [(x_hat, newick, si, entropy, tuple(key) + (r,)) for r in range(R)]
    own = pool is None
    if own:
        pool = ProcessPoolExecutor(workers, mp_context=get_context('spawn'))
    Us, stats, t_eng, t0 = [], [], 0.0, time.perf_counter()
    HE, nHE = 0.0, 0
    try:
        for r, (ds, st) in enumerate(pool.map(_sim_one, jobs)):
            t1 = time.perf_counter()
            eng = make_engine(ds)
            u = None
            if expected_h and (h_reps is None or r < h_reps):
                try:
                    _, Sp = eng.pair_scores(x_hat)
                    u = theta_rows(Sp, t_hat)
                except RuntimeError as e:              # exact fallback needed: no forward-mode pair scores
                    if 'fallback' not in str(e):
                        raise
            if u is not None:
                Us.append(u.sum(0))
                HE = HE + u.T @ u; nHE += 1
            else:
                _, gx = eng.lnl_grad(x_hat)
                Us.append(theta_rows(gx, t_hat)[0])
            stats.append(st)
            t_eng += time.perf_counter() - t1
    finally:
        if own:
            pool.shutdown()
    Uarr = np.array(Us)
    J = np.cov(Uarr.T, ddof=1)
    info = {'wall_s': time.perf_counter() - t0, 'engine_s': t_eng, 'R': R,
            'mean_N': float(np.mean([s['N'] for s in stats]))}
    if expected_h:
        info['H_E'] = (HE / nHE) if nHE else None
        info['H_E_R'] = nHE
    return J, Uarr, info


def pb_refits(make_engine, t_hat, newick, seq_index, R, entropy, key, refit, workers=8, pool=None, log=None):
    """Parametric bootstrap with full refits: the R datasets J_PB simulates at theta_hat (same key, so replicate
    r is the same dataset), each fitted by refit(eng, ds) -> {'theta', 'f', 'converged', 'outcome'}, the
    dataset's own M2 procedure over all parameters. Returns one row per replicate."""
    t_hat = np.asarray(t_hat, float)
    x_hat = fit.x_from_theta(t_hat)
    si = list(map(np.asarray, seq_index)) if isinstance(seq_index, (list, tuple)) else np.asarray(seq_index)
    jobs = [(x_hat, newick, si, entropy, tuple(key) + (r,)) for r in range(R)]
    own = pool is None
    if own:
        pool = ProcessPoolExecutor(workers, mp_context=get_context('spawn'))
    rows = []
    try:
        for r, (ds, st) in enumerate(pool.map(_sim_one, jobs)):
            t1 = time.perf_counter()
            try:
                q = refit(make_engine(ds), ds)
                rows.append({'r': r, 'theta': [float(v) for v in q['theta']], 'psi': float(q['theta'][fit.I_PSI]),
                             'f': float(q['f']), 'converged': bool(q['converged']), 'outcome': q['outcome'],
                             'N': st.get('N'), 'wall_s': time.perf_counter() - t1})
            except Exception as e:                         # one failed replicate does not end the bootstrap
                if fit._is_oom(e):
                    fit._free_gpu_cache()
                rows.append({'r': r, 'error': repr(e)[:300], 'wall_s': time.perf_counter() - t1})
            if log:
                log(rows[-1])
    finally:
        if own:
            pool.shutdown()
    return rows


def pb_percentile(rows, level=0.95):
    """Percentile interval for psi (equivalently 1/p = exp(psi)) from the replicate estimates, as the legacy
    parametric bootstrap built it; replicates whose fit errored or failed (optimizer_failure) are left out and
    counted."""
    ok = [q for q in rows if 'psi' in q and q.get('outcome') != 'optimizer_failure']
    a = (1 - level) / 2
    out = {'level': level, 'R': len(rows), 'n_used': len(ok),
           'n_error': sum('error' in q for q in rows),
           'n_failed': sum(q.get('outcome') == 'optimizer_failure' for q in rows)}
    if len(ok) < 10:
        return {**out, 'lo': None, 'hi': None}
    ps = np.array([q['psi'] for q in ok])
    lo, med, hi = np.quantile(ps, [a, 0.5, 1 - a])
    return {**out, 'lo': float(lo), 'hi': float(hi), 'median': float(med),
            'inv_p': [float(np.exp(lo)), float(np.exp(med)), float(np.exp(hi))],
            'n_at_lower_bound': int(np.sum(ps <= fit.PSI_MIN + 1e-6)),
            'n_at_upper_box': int(np.sum(ps >= fit.PSI_MAX - 1e-3))}


def lambda_pb_he(make_engine, t, newick, seq_index, R, entropy, key, workers=8, pool=None, fixed=(), h_reps=None):
    """Local lambda_hat at theta t (a restricted fit theta~(psi) or theta_hat) from J_PB and the expected
    Hessian H_E (per-pair Bartlett identity), both from the same R replicates simulated at t; no observed
    Hessian, so lambda_hat > 0 whenever H_E is positive definite. Returns dict(lambda_hat, V_psipsi,
    Hinv_psipsi, H_E_min_eig, R_used, wall_s)."""
    t = np.asarray(t, float)
    t0 = time.perf_counter()
    free = free_mask(t, fixed=fixed)
    J, _, info = j_pb(make_engine, t, newick, seq_index, R, entropy, key, workers, pool, expected_h=True,
                      h_reps=h_reps)
    HE = info['H_E']
    out = {'psi': float(t[I_PSI]), 'R_used': info['H_E_R'], 'lambda_hat': None}
    if HE is None:
        return {**out, 'wall_s': time.perf_counter() - t0}
    ev = np.linalg.eigvalsh(HE[np.ix_(free, free)])
    out['H_E_min_eig'] = float(ev[0])
    if ev[0] > 0:
        V, Hinv = godambe(HE, J, free)
        out.update(lambda_hat=float(V[I_PSI, I_PSI] / Hinv[I_PSI, I_PSI]), V_psipsi=float(V[I_PSI, I_PSI]),
                   Hinv_psipsi=float(Hinv[I_PSI, I_PSI]))
    out['wall_s'] = time.perf_counter() - t0
    return out


# ------------------------------------------------------------------ driver
def uq_dataset(eng, t_hat, newick=None, seq_index=None, R=100, entropy=None, key=(), make_engine=None,
               workers=8, pool=None, do_pb=True, expected_h=False):
    """All step-8 quantities for one fitted dataset. eng: engine on the observed data (pair_scores,
    site_scores, lnl_grad). Returns a JSON-ready dict."""
    t_hat = np.asarray(t_hat, float)
    x_hat = fit.x_from_theta(t_hat)
    fg = fit.theta_fg(eng)
    free = free_mask(t_hat, fixed=fixed_coords(eng.ds))
    T0 = time.perf_counter()
    H = hessian(fg, t_hat, free)
    t_H = time.perf_counter() - T0
    v, u, U = site_split(eng, x_hat, t_hat)
    t_sc = time.perf_counter() - T0 - t_H
    p_hat = float(np.exp(-t_hat[I_PSI]))
    Jh, hac = j_hac(v, eng.ds.nt_index, p_hat, locus=eng.ds.locus)
    Js = {'HAC': Jh, 'legacy': j_legacy(u)}
    if eng.ds.n_loci >= 20:
        Js['CR'], _ = j_cr(v, eng.ds.locus)
    out = {'free': free.tolist(), 'U_norm': float(np.linalg.norm(U[free])), 'hac': hac,
           'timing': {'hessian_s': t_H, 'scores_s': t_sc}}
    if do_pb:
        Jpb, Upb, info = j_pb(make_engine, t_hat, newick, seq_index, R, entropy, key, workers, pool,
                              expected_h=expected_h)
        HE = info.pop('H_E', None)
        Js['PB'] = Jpb
        out['pb'] = {**info, 'mean_U_psi': float(Upb[:, I_PSI].mean()),
                     'mcse_mean_U_psi': float(Upb[:, I_PSI].std(ddof=1) / np.sqrt(R))}
    ev = np.linalg.eigvalsh(H[np.ix_(free, free)])
    out['H'] = {'cond': float(ev[-1] / ev[0]) if ev[0] > 0 else None, 'min_eig': float(ev[0])}
    out['naive'] = summarize(t_hat, *godambe(H, H, free))          # V = H^-1 (J = H)
    for name, J in Js.items():
        V, Hinv = godambe(H, J, free)
        out[name] = summarize(t_hat, V, Hinv)
    out['matrices'] = {'H': H.tolist(), **{f'J_{k}': J.tolist() for k, J in Js.items()}}
    if do_pb and expected_h and HE is not None:
        # expected-information sandwich: H_E in place of the observed Hessian (PB = J_PB, HE_only = J = H_E)
        evE = np.linalg.eigvalsh(HE[np.ix_(free, free)])
        out['H_E'] = {'cond': float(evE[-1] / evE[0]) if evE[0] > 0 else None, 'min_eig': float(evE[0]),
                      'R': out['pb']['H_E_R']}
        out['PB_HE'] = summarize(t_hat, *godambe(HE, Js['PB'], free))
        out['naive_HE'] = summarize(t_hat, *godambe(HE, HE, free))
        out['matrices']['H_E'] = HE.tolist()
    return out
