"""Per-dataset inference after M2: M3 profile, outcome class, sandwich UQ, CLR intervals, p = 1 test
(PSIGC_GPU_PLAN.md steps 8-9). One call per fitted dataset; returns a JSON-ready dict."""
import time

import numpy as np

from . import clr, fit, uq

I_PSI = fit.I_PSI


def _lam(s):
    return s['lambda_hat'] if s and s.get('lambda_hat') is not None else None


def infer_dataset(eng, t_hat, cl_hat, tau0, newick, seq_index, R=100, entropy=None, key=(), make_engine=None,
                  workers=8, pool=None, L=None, do_pb=True, local_he=False, h_reps=None):
    """local_he: primary interval and p = 1 test by test inversion with a local lambda_hat(psi) from J_PB and
    the expected Hessian H_E at the restricted fits (clr.clr_interval_local; G5 diagnostic 2026-11-08: the
    observed-Hessian lambda at psi_hat under-covers and turns negative off the optimum). The observed-H
    outputs ('clr', 'boundary', 'lambda') are kept as sensitivity analyses."""
    T0 = time.perf_counter()
    t_hat = np.asarray(t_hat, float)
    fg = fit.theta_fg(eng)
    L = L if L is not None else int(np.bincount(eng.ds.locus).max())     # longest locus
    out = {'theta_m2': t_hat.tolist(), 'cl_m2': float(cl_hat)}
    # M3 profile; the plan's failure rule also compares with the best profile value
    prof = clr.m3_profile(fg, t_hat)
    best = max(prof, key=lambda p: p['f'])
    out['profile_refit'] = False
    if best['f'] > cl_hat + 0.01:
        r = fit.maximize(fg, best['theta'])
        if r['f'] > cl_hat:
            t_hat, cl_hat = r['theta'], float(r['f'])
            out.update(profile_refit=True, theta=t_hat.tolist(), cl=cl_hat, refit_converged=r['converged'])
            prof = clr.m3_profile(fg, t_hat)
    out['profile'] = [{'psi': p['psi'], 'f': p['f'], 'converged': p['converged'], 'nfev': p['nfev']} for p in prof]
    out['profile_nfev'] = int(sum(p['nfev'] for p in prof))
    out['outcome'] = fit.classify(t_hat, tau0, scan=prof)
    out['plateau'] = clr.plateau_sensitivity(t_hat, prof)
    t_prof = time.perf_counter() - T0
    # sandwich at theta_hat, or for plateau fits at the restricted fit psi = log L (sensitivity: log 1e3)
    if out['plateau']['base']:
        rf = clr.restricted_fit(fg, t_hat, np.log(L))
        t_uq = rf['theta']
        out['lambda_point'] = {'psi': float(np.log(L)), 'f': rf['f']}
    else:
        t_uq = t_hat
        out['lambda_point'] = {'psi': float(t_hat[I_PSI]), 'f': cl_hat}
    U = uq.uq_dataset(eng, t_uq, newick, seq_index, R, entropy, key, make_engine, workers, pool, do_pb=do_pb)
    if out['plateau']['base']:
        U['note'] = 'Wald summaries at the restricted fit psi = log L'
    out['uq'] = {k: v for k, v in U.items() if k != 'matrices'}
    out['matrices'] = U['matrices']
    t_uqs = time.perf_counter() - T0 - t_prof
    # CLR intervals: primary lambda from J_PB; sensitivity lambda from J_HAC (when valid) and lambda = 1
    lams = {'PB': _lam(U.get('PB')), 'HAC': _lam(U['HAC']) if U['hac']['valid'] else None, 'one': 1.0}
    if out['plateau']['base']:
        rf3 = clr.restricted_fit(fg, t_hat, np.log(1e3))
        U3 = uq.uq_dataset(eng, rf3['theta'], newick, seq_index, R, entropy, tuple(key) + (1,), make_engine,
                           workers, pool, do_pb=do_pb)
        lams['PB_at_1e3'] = _lam(U3.get('PB'))
    out['lambda'] = lams
    out['clr'] = {k: (clr.clr_interval(fg, t_hat, cl_hat, lam, prof) if lam is not None and lam > 0 else None)
                  for k, lam in lams.items()}
    out['boundary'] = {k: (clr.boundary_test(cl_hat, prof, lam) if lam is not None and lam > 0 else None)
                       for k, lam in lams.items()}
    t_clr = time.perf_counter() - T0 - t_prof - t_uqs
    if local_he and do_pb:
        fixed = uq.fixed_coords(eng.ds)
        cache = {}

        def lam_at(i):
            if i not in cache:
                cache[i] = uq.lambda_pb_he(make_engine, prof[i]['theta'], newick, seq_index, R, entropy,
                                           tuple(key) + (2, i), workers, pool, fixed=fixed, h_reps=h_reps)
            return cache[i]['lambda_hat']

        lh = uq.lambda_pb_he(make_engine, t_uq, newick, seq_index, R, entropy, tuple(key) + (3,), workers, pool,
                             fixed=fixed, h_reps=h_reps)
        lam_hat = lh['lambda_hat'] if lh['lambda_hat'] and lh['lambda_hat'] > 0 else 1.0
        loc = clr.clr_interval_local(fg, t_hat, cl_hat, prof, lam_at, lam_hat)
        W0 = max(2 * (cl_hat - prof[0]['f']), 0.0)
        lam0 = (lam_at(0) if W0 > 0 else None) or lam_hat
        out['local_he'] = {'h_reps': h_reps, 'lambda_hat': lh, 'lambda_hat_used': lam_hat, 'clr': loc,
                           'boundary': clr.boundary_test(cl_hat, prof, lam0), 'lambda_psi0_bound': lam0,
                           'points': {int(k): v for k, v in cache.items()},
                           'n_lambda_evals': 1 + len(cache)}
    out['timing'] = {'profile_s': t_prof, 'uq_s': t_uqs, 'clr_s': t_clr,
                     'local_he_s': time.perf_counter() - T0 - t_prof - t_uqs - t_clr,
                     'total_s': time.perf_counter() - T0}
    return out
