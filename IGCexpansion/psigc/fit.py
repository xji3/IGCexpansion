"""Fitting: theta parameterisation, L-BFGS-B with exact gradients, psi scans, profiles, Hessian,
convergence rules and outcome classes of SIMULATION_REDO_PLAN.md section 5.

theta = x with
  x[0:3] (legacy nested logs of nucleotide fractions) replaced by log(pi_A/pi_T, pi_C/pi_T, pi_G/pi_T),
         which is unconstrained and keeps every frequency positive at any optimizer step;
  x[6], x[7] (log eta, log p) replaced by (log tau, psi), psi = log(1/p) >= 0.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .model import freq_from_x, x_from_freq

I_LOGTAU, I_PSI = 6, 7
PSI_MAX = np.log(1e6)
PSI_MIN = 1e-4     # 1/p >= 1.0001: below this d/dpsi is ill-conditioned (PSIGC_GPU_PLAN.md section 3.3);
                   # at psi = 0 exactly every pair factorises and the psi-derivative through w_1 is lost


def _x3_from_t3(t3):
    e = np.exp(np.concatenate([t3, [0.0 * t3[0]]]))
    pi = e / e.sum()
    ag = pi[0] + pi[2]
    return np.log(np.stack([ag, pi[0] / ag, pi[1] / (1 - ag)]))


def _t3_jacobian(t3):
    """d x[0:3] / d theta[0:3] by complex step."""
    J = np.zeros((3, 3))
    for j in range(3):
        tc = np.array(t3, complex); tc[j] += 1e-30j
        J[:, j] = _x3_from_t3(tc).imag / 1e-30
    return J


def theta_from_x(x):
    t = np.array(x, float)
    pi = freq_from_x(np.asarray(x[:3], float))
    t[:3] = np.log(pi[:3] / pi[3])
    t[I_LOGTAU] = x[6] - x[7]
    t[I_PSI] = -x[7]
    return t


def x_from_theta(t):
    x = np.array(t, float)
    x[:3] = _x3_from_t3(np.asarray(t[:3], float))
    x[6] = t[I_LOGTAU] - t[I_PSI]
    x[7] = -t[I_PSI]
    return x


def grad_theta(gx, t=None, x=None):
    """Chain rule from the legacy-x gradient to theta. Needs theta (or x) for the frequency block."""
    g = np.array(gx, float)
    if t is None:
        t = theta_from_x(x)
    g[:3] = _t3_jacobian(np.asarray(t[:3], float)).T @ np.asarray(gx[:3], float)
    g[I_LOGTAU] = gx[6]
    g[I_PSI] = -gx[6] - gx[7]
    return g


def theta_fg(engine):
    """fg(theta) -> (lnL, d lnL / d theta) for any engine with lnl_grad(x)."""
    def fg(t):
        f, gx = engine.lnl_grad(x_from_theta(t))
        return f, grad_theta(gx, t=t)
    return fg


def default_bounds(n, psi=(PSI_MIN, PSI_MAX)):
    """Only psi is bounded. Do not box every variable: SciPy's L-BFGS-B then takes the raw gradient
    as its first step (no 1/||g|| scaling), which on this objective jumps tens of log units."""
    b = [(None, None)] * n
    b[I_PSI] = psi
    return b


def _projected_grad(t, g, bounds):
    pg = g.copy()
    for i, (lo, hi) in enumerate(bounds):
        if lo is not None and t[i] <= lo + 1e-12 and g[i] < 0:
            pg[i] = 0.0
        if hi is not None and t[i] >= hi - 1e-12 and g[i] > 0:
            pg[i] = 0.0
    return pg


ND_TOL = 1e-2     # Newton decrement threshold: <= 5e-5 lnL left to gain, <= 0.01 SE from the optimum


def _is_oom(e):
    return type(e).__name__ == 'OutOfMemoryError' or 'out of memory' in str(e)[:200].lower()


def _free_gpu_cache():
    import gc, sys
    gc.collect()
    if 'torch' in sys.modules:
        sys.modules['torch'].cuda.empty_cache()


def maximize(fg, t0, free=None, bounds=None, maxiter=500, pg_tol=1e-3, df_tol=1e-6, nd_tol=ND_TOL, maxcor=20):
    """Maximise fg(theta) -> (f, grad) over the `free` coordinates (default all) with L-BFGS-B.

    Converged (decision 3, PSIGC_GPU_PLAN.md section 10): Newton decrement sqrt(pg' H^-1 pg) <= nd_tol, with
    H^-1 L-BFGS-B's inverse-Hessian approximation of -cl, and |delta f| < df_tol over the last 3 iterations
    or L-BFGS-B's own convergence. If H^-1 is unavailable the unscaled rule (max |pg| <= pg_tol) is used.
    pg_tol also sets L-BFGS-B's gtol (pg_tol / 100). 'converged_pg' records the old unscaled rule."""
    t0 = np.array(t0, float)
    n = len(t0)
    free = np.arange(n) if free is None else np.asarray(free)
    bounds = default_bounds(n) if bounds is None else bounds
    fb = [bounds[i] for i in free]
    hist, cache = [], {}

    def full(z):
        t = t0.copy(); t[free] = z; return t

    def obj(z):
        key = z.tobytes()
        if key not in cache:
            oom = False
            try:
                with np.errstate(all='ignore'):
                    f, g = fg(full(z))
            except Exception as e:                      # rejected trial step (extreme rates on a line search)
                if not (isinstance(e, (ArithmeticError, ValueError)) or _is_oom(e)):
                    raise
                oom = _is_oom(e)
                f, g = -np.inf, np.zeros(n)
            if oom:                                     # after the traceback (and its tensors) is released
                _free_gpu_cache()
            if not np.isfinite(f) or not np.all(np.isfinite(g)):
                f, g = -1e300, np.zeros(n)              # forces the line search to backtrack
            cache.clear(); cache[key] = (f, g)
        f, g = cache[key]
        return -f, -g[free]

    def cb(z):
        hist.append(-obj(z)[0])

    res = minimize(obj, t0[free], jac=True, method='L-BFGS-B', bounds=fb, callback=cb,
                   options={'maxiter': maxiter, 'ftol': 1e-15, 'gtol': pg_tol * 1e-2, 'maxcor': maxcor})
    t = full(res.x)
    f, g = fg(t)
    pg = _projected_grad(t[free], g[free], fb)
    dfs = np.abs(np.diff(hist[-4:])) if len(hist) >= 4 else np.array([np.inf])
    flat = bool(dfs.max() < df_tol or res.success)
    conv_pg = bool(np.abs(pg).max() <= pg_tol and flat)
    sm = scaled_metrics(pg, res)
    converged = bool(sm['newton_dec'] <= nd_tol and flat) if sm['newton_dec'] is not None else conv_pg
    return {'theta': t, 'f': f, 'grad': g, 'max_pg': float(np.abs(pg).max()), 'converged': converged,
            'converged_pg': conv_pg, 'max_df3': float(dfs.max()), **sm,
            'nit': int(res.nit), 'nfev': int(res.nfev), 'message': str(res.message), 'history': hist}


def scaled_metrics(pg, res):
    """Scaled convergence measures (decision 3 chose the Newton decrement, PSIGC_GPU_PLAN.md section 10), from L-BFGS-B's inverse-Hessian approximation H^-1 of -cl:
    max_i |pg_i| sqrt(H^-1_ii) (gradient times an approximate standard error, in lnL units) and the Newton
    decrement sqrt(pg' H^-1 pg). Recorded for every fit so the rule can be applied after the fact; an exact
    variant uses the FD Hessian (hessian()) at the optimum."""
    try:
        Hi = res.hess_inv.todense() if hasattr(res.hess_inv, 'todense') else np.asarray(res.hess_inv)
        Hi = (Hi + Hi.T) / 2
        se = np.sqrt(np.clip(np.diag(Hi), 0, None))
        return {'max_pg_se': float(np.max(np.abs(pg) * se)), 'newton_dec': float(np.sqrt(max(pg @ Hi @ pg, 0.0)))}
    except Exception:
        return {'max_pg_se': None, 'newton_dec': None}


def rule_converged(r, rule='newton', tol=ND_TOL, df_tol=1e-6):
    """Re-apply a convergence rule to a maximize() result: 'newton' (the adopted rule: Newton decrement <= tol
    and |dcl| < df_tol over 3 iterations, or L-BFGS-B success), 'pg_se' (max |pg_i| SE_i <= tol, same dcl
    condition), 'plan_unscaled' (the old rule as maximize() evaluated it, with its pg_tol)."""
    if rule == 'plan_unscaled':
        return r.get('converged_pg', r['converged'])
    v = r['max_pg_se'] if rule == 'pg_se' else r['newton_dec']
    if v is None:
        return False
    return bool(v <= tol and (r['max_df3'] < df_tol or r['message'].startswith('CONVERGENCE')))


def psi_scan(fg, t0, psis, free=(I_LOGTAU,)):
    """For each psi, re-optimise the `free` coordinates (default log tau only), warm-started."""
    out, t = [], np.array(t0, float)
    for psi in psis:
        t = t.copy(); t[I_PSI] = psi
        r = maximize(fg, t, free=list(free))
        t = r['theta']
        out.append({'psi': float(psi), 'f': r['f'], 'theta': t, 'converged': r['converged']})
    return out


def profile(fg, t0, psis):
    """M3 profile: all parameters except psi re-optimised at each psi."""
    n = len(t0)
    return psi_scan(fg, t0, psis, free=[i for i in range(n) if i != I_PSI])


def hessian(fg, t, h=1e-4, idx=None):
    """Central differences of the exact gradient, symmetrised (2 gradient calls per parameter)."""
    n = len(t)
    idx = range(n) if idx is None else idx
    H = np.zeros((n, n))
    for j in idx:
        tp, tm = t.copy(), t.copy(); tp[j] += h; tm[j] -= h
        H[:, j] = (fg(tp)[1] - fg(tm)[1]) / (2 * h)
    return (H + H.T) / 2


def classify(theta_hat, tau0, scan=None):
    psi, tau = theta_hat[I_PSI], np.exp(theta_hat[I_LOGTAU])
    if tau < 1e-3 * tau0:
        return 'tau_to_zero'
    if psi < 1e-3:
        return 'at_bound'
    if psi >= np.log(1e4) and scan is not None:
        pts = sorted((s['psi'], s['f']) for s in scan if s['psi'] <= np.log(1e5) + 1e-9)
        fs = [f for _, f in pts]
        if all(b >= a - 1e-6 for a, b in zip(fs, fs[1:])):
            return 'plateau'
    return 'interior'


def fit_m1(is_eng, ps_eng, t0):
    """M1: IS-IGC MLE (psi irrelevant), then psi only on the PS likelihood. Returns (stage1, stage2)."""
    xt = x_from_theta
    n = len(t0)
    fg_is = theta_fg(is_eng)
    s1 = maximize(fg_is, t0, free=[i for i in range(n) if i != I_PSI])
    f_ps = lambda psi: -ps_eng.lnl(xt(np.r_[s1['theta'][:I_PSI], psi, s1['theta'][I_PSI + 1:]]))
    r = minimize_scalar(f_ps, bounds=(PSI_MIN, PSI_MAX), method='bounded', options={'xatol': 1e-3})
    t = s1['theta'].copy(); t[I_PSI] = r.x
    return s1, {'theta': t, 'f': -r.fun, 'nfev': int(r.nfev)}


def fit_m2(ps_eng, t_start, tau0=None, n_scan=15, alt_starts=(), theta0=None, fixed=()):
    """M2: 15-point psi scan on [PSI_MIN, log 1e6] re-optimising tau, then L-BFGS-B on all parameters from
    the best scan point and any alternative starts. Refit rule (SIMULATION_REDO_PLAN.md section 5): a fit
    fails if it is not converged or cl(theta_hat) < max(cl(theta0), best scan value) - 0.01; up to 2 refits,
    from theta0 (when given) and from the best scan point. Calls ps_eng.begin_task() first if it has one
    (the K ratchet of the adaptive engines restarts per task). `fixed`: coordinates held at their start values
    (uq.fixed_coords: r2, r3 for non-coding data, which cl does not depend on; left free they have zero gradient
    and zero curvature, which corrupts L-BFGS-B's inverse Hessian and so the Newton-decrement test)."""
    if hasattr(ps_eng, 'begin_task'):
        ps_eng.begin_task()
    fg = theta_fg(ps_eng)
    psis = np.linspace(PSI_MIN, PSI_MAX, n_scan)
    scan = psi_scan(fg, t_start, psis)
    best = max(scan, key=lambda s: s['f'])
    starts = [best['theta']] + [np.asarray(a, float) for a in alt_starts]
    free = [i for i in range(len(t_start)) if i not in set(fixed)]
    fits = [maximize(fg, s, free=free) for s in starts]
    top = max(fits, key=lambda r: r['f'])
    f0 = ps_eng.lnl(x_from_theta(theta0)) if theta0 is not None else -np.inf
    target = max(f0, best['f'])
    refit_from = ([np.asarray(theta0, float)] if theta0 is not None else [top['theta']]) + [best['theta']]
    refit = 0
    while (not top['converged'] or top['f'] < target - 0.01) and refit < 2:
        r = maximize(fg, refit_from[refit], free=free)
        refit += 1
        if r['f'] > top['f'] or (r['converged'] and not top['converged']):
            top = r
    tau0 = tau0 if tau0 is not None else np.exp(t_start[I_LOGTAU])
    failed = (not top['converged']) or top['f'] < target - 0.01
    return {'fit': top, 'scan': scan, 'all_fits': fits, 'refits': refit, 'cl_theta0': f0,
            'outcome': 'optimizer_failure' if failed else classify(top['theta'], tau0, scan)}
