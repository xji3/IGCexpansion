"""Collect every number of the results report into one JSON.
  python scripts/report_data.py <results dir> <out.json>
<results dir> holds armL_UQ, armL_PB, realdata_fits and realdata_pb as written by the run commands in README.md."""
import json, glob, math, os, sys
import numpy as np
R = sys.argv[1]
cells = [3, 10, 50, 100, 200, 300, 400, 500]
UB = math.log(1e4)
A = 0.05


def get(d, path):
    for k in path:
        d = d.get(k) if isinstance(d, dict) else None
    return d


M = {
    'naive': ['uq', 'naive', 'psi'],
    'sandwich': ['uq', 'PB', 'psi'],
    'clr_global': ['clr', 'PB'],
    'clr_local': ['local_he', 'clr'],
}


def ivl(r, m):
    d = get(r, M[m])
    if not d:
        return None
    lo, hi = d.get('lo'), d.get('hi')
    if lo is None or not np.isfinite(lo):
        return None
    return float(lo), (float(hi) if hi is not None and np.isfinite(hi) else None)


def score(lo, hi, x):
    hi_ = UB if hi is None else hi
    return (hi_ - lo) + (2 / A) * max(lo - x, 0) + (2 / A) * max(x - hi_, 0), hi_ - lo, (lo <= x <= hi_), x < lo, x > hi_


rows = {}
for c in cells:
    for f in sorted(glob.glob(f'{R}/armL_UQ/UQ/tract_{c}/sim_*.json')):
        d = json.load(open(f))
        # sim_88 at 1/p = 400 and 500 is included: the verified 2017 simulator gives the same alignment for both
        # cells (shared seed, tracts covering the same differing sites), so it is a genuine draw in each
        r = d['result']; rep = d['meta']['rep']
        e = {m: ivl(r, m) for m in M}
        e['psi_hat'] = r['lambda_point']['psi']
        e['boundary_p'] = get(r, ['local_he', 'boundary', 'p_value'])
        pf = f'{R}/armL_PB/PB/tract_{c}/sim_{rep}.json'
        if os.path.exists(pf):
            iv = json.load(open(pf))['result']['interval']
            e['pb'] = (iv['lo'], iv['hi'])
        rows[(c, rep)] = e


def summarize(keys, methods):
    out = {}
    for m in methods:
        S = [score(*rows[k][m], math.log(k[0])) for k in keys if rows[k].get(m)]
        if not S:
            continue
        s = np.array([q[0] for q in S]); w = np.array([q[1] for q in S])
        out[m] = {'n': len(S), 'n_missing': len(keys) - len(S), 'cover': float(np.mean([q[2] for q in S])),
                  'width_mean': float(w.mean()), 'width_med': float(np.median(w)),
                  'is_mean': float(s.mean()), 'is_med': float(np.median(s)),
                  'penalty_mean': float((s - w).mean()),
                  'miss_low': int(sum(q[3] for q in S)), 'miss_high': int(sum(q[4] for q in S))}
    return out


meths = list(M)
D = {'cells': {}, 'pb_cells': {}, 'pb_pooled': None, 'pb_cases': [], 'realdata': []}
for c in cells:
    ks = [k for k in rows if k[0] == c]
    ph = np.array([rows[k]['psi_hat'] for k in ks])
    bp = [rows[k]['boundary_p'] for k in ks if rows[k]['boundary_p'] is not None]
    D['cells'][c] = {'n': len(ks), 'methods': summarize(ks, meths),
                     'psi_hat_median_invp': float(math.exp(np.median(ph))),
                     'psi_hat_q25_invp': float(math.exp(np.quantile(ph, .25))),
                     'psi_hat_q75_invp': float(math.exp(np.quantile(ph, .75))),
                     'frac_psi_hat_below': float(np.mean(ph < math.log(c))),
                     'boundary_reject': float(np.mean(np.array(bp) < 0.05)) if bp else None}
pk = sorted(k for k in rows if 'pb' in rows[k])
for c in cells:
    D['pb_cells'][c] = summarize([k for k in pk if k[0] == c], meths + ['pb'])
D['pb_pooled'] = summarize(pk, meths + ['pb'])
for k in pk:
    e = rows[k]
    D['pb_cases'].append({'cell': k[0], 'rep': k[1], 'psi_hat': e['psi_hat'], **{m: e[m] for m in meths + ['pb']}})

for f in sorted(glob.glob(f'{R}/realdata_fits/*.json')):
    d = json.load(open(f)); n = d['name']; inf = d['infer']
    pf = f'{R}/realdata_pb/{n}.pb.json'
    pb = None
    if os.path.exists(pf):
        P = json.load(open(pf)); iv = P.get('interval_recheck') or P['interval']
        pb = {k: iv.get(k) for k in ('lo', 'hi', 'median', 'n_used', 'n_failed', 'n_at_lower_bound')}
        pb['rechecked'] = 'interval_recheck' in P
    D['realdata'].append({
        'name': n, 'noncoding': d['noncoding'], 'n_loci': d['n_loci'], 'summary': d['summary'],
        'm2_outcome': d['m2']['outcome'], 'psi_hat': inf['lambda_point']['psi'],
        'sandwich': ivl(inf, 'sandwich'), 'clr_global': ivl(inf, 'clr_global'), 'clr_local': ivl(inf, 'clr_local'),
        'boundary_p': get(inf, ['local_he', 'boundary', 'p_value']), 'pb': pb})
json.dump(D, open(sys.argv[2], 'w'), indent=1, default=str)
print('ok', len(rows), len(pk), len(D['realdata']))
