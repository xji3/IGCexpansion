"""Profile-based inference (SIMULATION_REDO_PLAN.md sections 5-6; PSIGC_GPU_PLAN.md step 9).

* M3 profile: psi on log 10^{0, 0.5, ..., 5} (the first point is PSI_MIN), every other parameter re-maximised,
  warm-started outward from the fit with a linear predictor of the nuisance parameters (evaluation reduction).
* W(psi) = 2 [cl(theta_hat) - cl_prof(psi)].
* Adjusted profile CLR interval {psi >= 0 : W(psi) <= 3.84 lambda_hat}: endpoints by Brent root-finding on
  W(psi) - 3.84 lambda_hat (xtol 0.02 in psi) inside the bracket the profile grid gives, the nuisance
  parameters re-maximised at each trial from the interpolated bracket ends. The upper limit is +inf when W
  stays below the cutoff up to the last profile point; the lower limit is 0 when W(PSI_MIN) is below it.
* Test of p = 1: W(0) / lambda_hat against 1/2 chi2_0 + 1/2 chi2_1 (reject above 2.71). psi = 0 is evaluated
  at PSI_MIN (1/p = 1.0001), the bound of the fit.
* Plateau rule: psi_hat >= log 1e4 and the profile non-decreasing (within tol) up to log 1e5; sensitivity
  variants with other tolerances and thresholds. For plateau fits lambda_hat is computed at the restricted
  fit at psi = log L (sensitivity: lambda_hat = 1, psi = log 1e3, a bootstrap cutoff).
* Negative controls: the p-only parametric bootstrap (psi re-estimated with everything else fixed at
  theta_hat, as the legacy analyses did) and bootstrap refits for the warp-speed and nested checks.
"""
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

from . import fit

I_PSI = fit.I_PSI
CUT95 = 3.841458820694124
BOUNDARY_CRIT = 2.705543454095404          # 1/2 chi2_0 + 1/2 chi2_1 at 5%
PROFILE_LOG10 = np.arange(0, 5.01, 0.5)


def profile_grid():
    g = np.log(10.0 ** PROFILE_LOG10)
    g[0] = fit.PSI_MIN
    return g


def _nuisance(n):
    return [i for i in range(n) if i != I_PSI]


def _prof_point(fg, t_start, psi):
    t = np.array(t_start, float); t[I_PSI] = psi
    r = fit.maximize(fg, t, free=_nuisance(len(t)))
    return {'psi': float(psi), 'f': float(r['f']), 'theta': r['theta'], 'converged': r['converged'],
            'nfev': r['nfev']}


def m3_profile(fg, t_hat, psis=None, predictor=True):
    """Profile at psis (default profile_grid()), walking outward from the grid point nearest psi_hat in both
    directions. With predictor, each start extrapolates the nuisance parameters linearly in psi from the
    previous two points of the walk."""
    t_hat = np.asarray(t_hat, float)
    psis = profile_grid() if psis is None else np.asarray(psis, float)
    i0 = int(np.argmin(np.abs(psis - t_hat[I_PSI])))
    pts = {}
    for direction in (range(i0, len(psis)), range(i0 - 1, -1, -1)):
        prev = [(t_hat[I_PSI], t_hat)]
        for i in direction:
            if len(prev) >= 2 and predictor:
                (p1, t1), (p2, t2) = prev[-2], prev[-1]
                start = t2 + (t2 - t1) * (psis[i] - p2) / (p2 - p1) if abs(p2 - p1) > 1e-12 else t2
                start = np.clip(start, -30, 30)
            else:
                start = prev[-1][1]
            q = _prof_point(fg, start, psis[i])
            if predictor and len(prev) >= 2:                      # guard: predictor start worse than plain warm start
                q2 = _prof_point(fg, prev[-1][1], psis[i]) if not q['converged'] else None
                if q2 is not None and q2['f'] > q['f']:
                    q = q2
            pts[i] = q
            prev.append((psis[i], q['theta']))
    return [pts[i] for i in range(len(psis))]


def W_values(cl_hat, prof):
    return np.array([2 * (cl_hat - p['f']) for p in prof])


def _root(fg, cl_hat, cut, a, b, nfev):
    """Root of W(psi) - cut between profile points a and b (W(a) - cut and W(b) - cut of opposite sign)."""
    pa, pb = a['psi'], b['psi']

    def g(psi):
        lam = (psi - pa) / (pb - pa)
        q = _prof_point(fg, (1 - lam) * a['theta'] + lam * b['theta'], psi)
        nfev.append(q['nfev'])
        return 2 * (cl_hat - q['f']) - cut

    return brentq(g, pa, pb, xtol=0.02, rtol=1e-6)


def clr_interval(fg, t_hat, cl_hat, lam, prof, cut=CUT95):
    """Adjusted profile CLR interval for psi. Returns dict with lo, hi (hi = inf allowed), the cutoff used,
    and the number of optimizer evaluations spent on the endpoints."""
    c = cut * lam
    psi_hat = float(np.asarray(t_hat)[I_PSI])
    W = W_values(cl_hat, prof)
    ps = np.array([p['psi'] for p in prof])
    nfev = []
    hit = {'psi': psi_hat, 'f': cl_hat, 'theta': np.asarray(t_hat, float)}
    # upper limit: first grid point above psi_hat with W > c
    above = [i for i in range(len(ps)) if ps[i] > psi_hat]
    hi = np.inf
    prev = hit
    for i in above:
        if W[i] > c:
            hi = _root(fg, cl_hat, c, prev, prof[i], nfev)
            break
        prev = prof[i]
    below = [i for i in range(len(ps) - 1, -1, -1) if ps[i] < psi_hat]
    lo = 0.0
    prev = hit
    for i in below:
        if W[i] > c:
            lo = _root(fg, cl_hat, c, prof[i], prev, nfev)
            break
        prev = prof[i]
    return {'lo': float(lo), 'hi': float(hi), 'cutoff': float(c), 'endpoint_nfev': int(sum(nfev)),
            'W_min': float(W.min()), 'refit_needed': bool(W.min() < -0.02)}


def clr_interval_local(fg, t_hat, cl_hat, prof, lam_at, lam_hat, cut=CUT95):
    """Test-inversion CLR with a local lambda: {psi : W(psi) <= cut * lambda(psi)}, lambda(psi) at the
    restricted fit theta~(psi). lam_at(i) gives lambda at profile point i (cached by the caller), lam_hat
    lambda at theta_hat. On each side the first grid point outside the region is found by bisection over
    the grid indices (assumes one crossing per side; every probed point is recorded so this can be
    checked), then the endpoint is the root of W(psi) - cut * lambda(psi) with W from re-maximised profile
    points and log lambda interpolated linearly in psi between the bracket ends. lam_at may return None
    (no positive-definite H_E): that point is treated as inside (conservative)."""
    psi_hat = float(np.asarray(t_hat)[I_PSI])
    W = W_values(cl_hat, prof)
    ps = np.array([p['psi'] for p in prof])
    nfev, probes = [], {}

    def outside(i):
        lam = lam_at(i)
        probes[i] = lam
        return lam is not None and W[i] > cut * lam

    hat = {'psi': psi_hat, 'f': cl_hat, 'theta': np.asarray(t_hat, float)}

    def side(idx):
        """idx: grid indices ordered outward from psi_hat. Returns (bracket inner, outer) or None."""
        if not idx or not outside(idx[-1]):
            return None
        lo, hi = -1, len(idx) - 1               # inside at lo (psi_hat when -1), outside at hi
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if outside(idx[mid]):
                hi = mid
            else:
                lo = mid
        return (None if lo < 0 else idx[lo]), idx[hi]

    def root(i_in, i_out):
        a = hat if i_in is None else prof[i_in]
        b = prof[i_out]
        la = lam_hat if i_in is None else (probes.get(i_in) or lam_hat)
        lb = probes[i_out]
        pa, pb = a['psi'], b['psi']

        def g(psi):
            w = (psi - pa) / (pb - pa)
            q = _prof_point(fg, (1 - w) * a['theta'] + w * b['theta'], psi)
            nfev.append(q['nfev'])
            lam = float(np.exp((1 - w) * np.log(la) + w * np.log(lb)))
            return 2 * (cl_hat - q['f']) - cut * lam

        ga, gb = -cut * la if i_in is None else W[i_in] - cut * la, W[i_out] - cut * lb
        if ga >= 0:                               # interpolation inconsistent at the inner end: use the grid point
            return pa
        return brentq(g, min(pa, pb), max(pa, pb), xtol=0.02, rtol=1e-6) if gb > 0 else pb

    up = [i for i in range(len(ps)) if ps[i] > psi_hat]
    down = [i for i in range(len(ps) - 1, -1, -1) if ps[i] < psi_hat]
    bu, bd = side(up), side(down)
    hi = np.inf if bu is None else root(*bu)
    lo = 0.0 if bd is None else root(*bd)
    return {'lo': float(lo), 'hi': float(hi), 'endpoint_nfev': int(sum(nfev)),
            'probes': {int(k): v for k, v in probes.items()}, 'W': W.tolist(),
            'refit_needed': bool(W.min() < -0.02)}


def boundary_test(cl_hat, prof, lam):
    """W(0) / lambda_hat against 1/2 chi2_0 + 1/2 chi2_1; prof[0] must be the PSI_MIN point."""
    assert abs(prof[0]['psi'] - fit.PSI_MIN) < 1e-12
    W0 = max(2 * (cl_hat - prof[0]['f']), 0.0)
    stat = W0 / lam
    return {'W0': float(W0), 'stat': float(stat), 'p_value': float(0.5 * chi2.sf(stat, 1)) if stat > 0 else 0.5,
            'reject_5pct': bool(stat > BOUNDARY_CRIT)}


def plateau(t_hat, prof, psi_min=np.log(1e4), psi_end=np.log(1e5), tol=1e-3):
    """Plateau rule: psi_hat >= psi_min and the profile non-decreasing (within tol lnL) up to psi_end."""
    if np.asarray(t_hat)[I_PSI] < psi_min:
        return False
    fs = [p['f'] for p in prof if p['psi'] <= psi_end + 1e-9]
    return bool(all(b >= a - tol for a, b in zip(fs, fs[1:])))


def plateau_sensitivity(t_hat, prof):
    return {'base': plateau(t_hat, prof),
            'tol_1e-6': plateau(t_hat, prof, tol=1e-6), 'tol_1e-2': plateau(t_hat, prof, tol=1e-2),
            'psi_min_1e3': plateau(t_hat, prof, psi_min=np.log(1e3)),
            'from_1e3': plateau(t_hat, [p for p in prof if p['psi'] >= np.log(1e3) - 1e-9])}


def restricted_fit(fg, t_hat, psi):
    """Fit with psi fixed (used for lambda_hat of plateau fits at psi = log L, sensitivity at log 1e3)."""
    return _prof_point(fg, t_hat, psi)


# ------------------------------------------------------------------ bootstrap-based checks
def p_only_fit(eng, t_fixed, n_grid=15):
    """psi only, every other parameter at t_fixed (legacy-style): grid on [PSI_MIN, PSI_MAX] then bounded
    Brent around the best grid point."""
    t_fixed = np.asarray(t_fixed, float)

    def f(psi):
        t = t_fixed.copy(); t[I_PSI] = psi
        return eng.lnl(fit.x_from_theta(t))

    psis = np.linspace(fit.PSI_MIN, fit.PSI_MAX, n_grid)
    vals = [f(p) for p in psis]
    i = int(np.argmax(vals))
    lo, hi = psis[max(i - 1, 0)], psis[min(i + 1, n_grid - 1)]
    r = minimize_scalar(lambda p: -f(p), bounds=(lo, hi), method='bounded', options={'xatol': 1e-3})
    return (float(r.x), float(-r.fun)) if -r.fun >= vals[i] else (float(psis[i]), float(vals[i]))


def bootstrap(make_engine, t_hat, newick, seq_index, B, entropy, key, kind, workers=8, pool=None):
    """Parametric bootstrap at theta_hat. kind: 'p_only' (psi re-estimated, rest fixed at theta_hat; the
    legacy negative control) or 'full' (M2 refit on all parameters, started at theta_hat; warp-speed and
    nested checks). Returns psi* (B,) and per-replicate details."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    from .uq import _sim_one
    t_hat = np.asarray(t_hat, float)
    x_hat = fit.x_from_theta(t_hat)
    si = list(map(np.asarray, seq_index)) if isinstance(seq_index, (list, tuple)) else np.asarray(seq_index)
    jobs = [(x_hat, newick, si, entropy, tuple(key) + (b,)) for b in range(B)]
    own = pool is None
    if own:
        pool = ProcessPoolExecutor(workers, mp_context=get_context('spawn'))
    rows, t0 = [], time.perf_counter()
    try:
        for ds, st in pool.map(_sim_one, jobs):
            eng = make_engine(ds)
            if kind == 'p_only':
                psi, f = p_only_fit(eng, t_hat)
                rows.append({'psi': psi, 'f': f, 'N': st['N']})
            elif kind == 'full':
                r = fit.fit_m2(eng, t_hat, tau0=float(np.exp(t_hat[fit.I_LOGTAU])), theta0=t_hat)
                rows.append({'psi': float(r['fit']['theta'][I_PSI]), 'f': float(r['fit']['f']),
                             'outcome': r['outcome'], 'theta': r['fit']['theta'].tolist(), 'N': st['N']})
            else:
                raise ValueError(kind)
    finally:
        if own:
            pool.shutdown()
    psi = np.array([r['psi'] for r in rows])
    return psi, {'rows': rows, 'wall_s': time.perf_counter() - t0, 'sd_psi': float(psi.std(ddof=1)) if B > 1 else None}


def percentile_interval(psi_star, level=0.95):
    a = (1 - level) / 2
    return float(np.quantile(psi_star, a)), float(np.quantile(psi_star, 1 - a))
