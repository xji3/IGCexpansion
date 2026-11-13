"""Arm L of SIMULATION_REDO_PLAN.md: M0, M1/M4 and M2 on the 800 archived replicates.

Subcommands
  bundle  --src <IGCTract/SimulationStudy> --out <dir>   copy FASTAs, event-log stats, tree, seq_index, truth
  run     --bundle <dir> --out <dir> --workers N [--estimators M0,M1,M2] [--cells 3,10] [--reps 1-5]
          [--backend cpu|gpu --devices 0,1 --workers-per-gpu 1]
  status  --bundle <dir> --out <dir>
  collect --out <dir>                                    all rows -> rows.csv (and rows.parquet if pyarrow)

Each task is one (dataset, estimator); it writes <out>/<estimator>/tract_<g>/sim_<k>.json atomically and
is skipped if that file exists, so the runner can be killed and restarted at any time. Each JSON holds one
row per start (plan section 7, orchestration) plus the scan/grid used.
"""
import argparse, hashlib, json, os, shutil, sys, time, traceback

for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from IGCexpansion import psigc                                           # noqa: E402
from IGCexpansion.psigc import fit                                  # noqa: E402

CELLS = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
REPS = range(1, 101)
PARALOGS = ('YDR418W', 'YEL054C')
K_NODES = 'auto'       # adaptive Chebyshev K (ladder 16/24/32/48, ratchet per M2 task, exact fallback)
PG_TOL = 1e-2          # L-BFGS-B gtol = PG_TOL/100; convergence itself is the Newton decrement <= fit.ND_TOL
N_GRID = 15


def code_hash():
    """Hash of every engine source (psigc/*.py, including the GPU engine psigc/gpu.py) and this runner."""
    h = hashlib.sha256()
    from IGCexpansion import psigc as _psigc
    pkg = os.path.dirname(os.path.abspath(_psigc.__file__))
    for p in [os.path.join(pkg, f) for f in sorted(os.listdir(pkg))] + [os.path.abspath(__file__)]:
        if p.endswith('.py') and os.path.exists(p):
            h.update(os.path.basename(p).encode())
            h.update(open(p, 'rb').read())
    return h.hexdigest()[:16]


def provenance(backend, device):
    """Software and hardware versions recorded in every result file."""
    import platform, scipy
    v = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
         'psigc': psigc.__version__}
    if backend == 'gpu':
        import torch
        v.update(torch=torch.__version__, cuda=torch.version.cuda, gpu=torch.cuda.get_device_name(device),
                 cudnn=torch.backends.cudnn.version())
        try:
            import subprocess
            v['driver'] = subprocess.run(['nvidia-smi', '--query-gpu=driver_version', '--format=csv,noheader'],
                                         capture_output=True, text=True, timeout=10).stdout.split()[0]
        except Exception:
            pass
    return v


# ------------------------------------------------------------------ bundle
def event_stats(log):
    rows = [l.rstrip('\n').split('\t') for l in open(log)][1:]
    rows = [r for r in rows if len(r) >= 8]
    lens = [int(r[6]) - int(r[5]) + 1 for r in rows]
    return {'N': len(rows), 'N_inf': int(sum(int(r[7]) > 0 for r in rows)),
            'm_real': float(np.mean(lens)) if lens else None}


def cmd_bundle(a):
    src, out = a.src, a.out
    os.makedirs(os.path.join(out, 'data'), exist_ok=True)
    shutil.copy(os.path.join(src, 'YeastTree.newick'), os.path.join(out, 'YeastTree.newick'))
    shutil.copy(os.path.join(src, '..', 'MafftAlignment', 'YDR418W_YEL054C', 'YDR418W_YEL054C_seq_index.txt'),
                os.path.join(out, 'seq_index.txt'))
    is_save = np.loadtxt(os.path.join(src, 'save', 'YDR418W_YEL054C_HKY_rv_IGC_nonclock_save.txt'))
    assert len(is_save) == 19
    man, hashes = [], {}
    for g in CELLS:
        for k in REPS:
            d = os.path.join(src, f'Tract_{g}_HKY', f'sim_{k}')
            fa = os.path.join(d, f'YDR418W_YEL054C_sim_{k}.fasta')
            dst = os.path.join('data', f'tract_{g:g}', f'sim_{k}.fasta')
            os.makedirs(os.path.dirname(os.path.join(out, dst)), exist_ok=True)
            shutil.copy(fa, os.path.join(out, dst))
            hsh = hashlib.sha256(open(fa, 'rb').read()).hexdigest()
            dup = hashes.get(hsh)
            hashes.setdefault(hsh, (g, k))
            # truth in the legacy PS x layout: x_pm, log eta = log(tau / g), log p = -log g, edge rates
            x0 = np.concatenate([is_save[:6], [is_save[6] - np.log(g), -np.log(g)], is_save[7:]])
            man.append({'cell': g, 'rep': k, 'fasta': dst, 'sha256': hsh, 'duplicate_of': dup,
                        # sim_88 is identical at 1/p = 400 and 500 but a genuine draw in both (the 2017 simulator gives
                        # it in each cell: shared seed, tracts covering the same differing sites), so nothing is excluded
                        'excluded': False,
                        'truth_x': x0.tolist(), **event_stats(os.path.join(d, f'YDR418W_YEL054C_sim_{k}_IGC.log'))})
    with open(os.path.join(out, 'manifest.jsonl'), 'w') as f:
        for m in man:
            f.write(json.dumps(m) + '\n')
    print(f'{len(man)} datasets; duplicates: {[(m["cell"], m["rep"], m["duplicate_of"]) for m in man if m["duplicate_of"]]}')


# ------------------------------------------------------------------ estimators
def load(bundle, m):
    tree = psigc.Tree(open(os.path.join(bundle, 'YeastTree.newick')).read())
    ds = psigc.Dataset.from_fasta(os.path.join(bundle, m['fasta']), os.path.join(bundle, 'seq_index.txt'), tree, PARALOGS)
    return tree, ds


def grid_brent(f, n=N_GRID):
    """Maximise a scalar function of psi: grid on [PSI_MIN, PSI_MAX], then bounded Brent around the best."""
    from scipy.optimize import minimize_scalar
    psis = np.linspace(fit.PSI_MIN, fit.PSI_MAX, n)
    vals = [f(p) for p in psis]
    i = int(np.argmax(vals))
    lo, hi = psis[max(i - 1, 0)], psis[min(i + 1, n - 1)]
    r = minimize_scalar(lambda p: -f(p), bounds=(lo, hi), method='bounded', options={'xatol': 1e-3})
    best = (float(r.x), float(-r.fun)) if -r.fun >= vals[i] else (float(psis[i]), float(vals[i]))
    return best, [{'psi': float(p), 'f': float(v)} for p, v in zip(psis, vals)], int(r.nfev) + n


def theta_row(t):
    return {'theta': [float(v) for v in t], 'psi': float(t[fit.I_PSI]), 'inv_p': float(np.exp(t[fit.I_PSI])),
            'tau': float(np.exp(t[fit.I_LOGTAU])), 'eta': float(np.exp(t[fit.I_LOGTAU] - t[fit.I_PSI]))}


def est_M0(eng, is_eng, t0, tau0):
    def f(psi):
        t = t0.copy(); t[fit.I_PSI] = psi
        t[fit.I_LOGTAU] = t0[fit.I_LOGTAU]                       # tau fixed at truth (eta = tau p)
        return eng.lnl(fit.x_from_theta(t))
    (psi, lnl), grid, nfev = grid_brent(f)
    t = t0.copy(); t[fit.I_PSI] = psi
    return [{'start': 'grid', 'lnL': lnl, **theta_row(t), 'nfev': nfev,
             'outcome': fit.classify(t, tau0, grid), 'grid': grid}]


def is_fit(is_eng, t0):
    fg = fit.theta_fg(is_eng)
    free = [i for i in range(len(t0)) if i != fit.I_PSI]
    return fit.maximize(fg, t0, free=free, pg_tol=PG_TOL)


def est_M1(eng, is_eng, t0, tau0):
    s1 = is_fit(is_eng, t0)
    ts = s1['theta']

    def f(psi):
        t = ts.copy(); t[fit.I_PSI] = psi
        return eng.lnl(fit.x_from_theta(t))
    (psi, lnl), grid, nfev = grid_brent(f)
    t = ts.copy(); t[fit.I_PSI] = psi
    m4 = {'start': 'truth', 'estimator': 'M4', 'lnL_IS': s1['f'], **theta_row(ts), 'converged': s1['converged'],
          'max_pg': s1['max_pg'], 'nit': s1['nit'], 'nfev': s1['nfev']}
    m4['psi'] = m4['inv_p'] = m4['eta'] = None                 # IS fit: tau only
    m1 = {'start': 'IS_fit', 'estimator': 'M1', 'lnL': lnl, **theta_row(t), 'nfev': nfev,
          'outcome': fit.classify(t, tau0, grid), 'grid': grid}
    return [m4, m1]


def fit_row(name, r, t1):
    return {'start': name, 'lnL': r['f'], **theta_row(r['theta']), 'converged': r['converged'],
            'max_pg': r['max_pg'], 'max_pg_se': r.get('max_pg_se'), 'newton_dec': r.get('newton_dec'),
            'max_df3': r.get('max_df3'), 'nit': r['nit'], 'nfev': r['nfev'], 'message': r['message'],
            'wall_s': time.perf_counter() - t1}


def est_M2(eng, is_eng, t0, tau0):
    fg = fit.theta_fg(eng)
    s1 = is_fit(is_eng, t0)
    ts = s1['theta']
    if hasattr(eng, 'begin_task'):
        eng.begin_task()                 # K ratchet: IS warm start excluded; scan and all starts included
    t_sc = time.perf_counter()
    scan = fit.psi_scan(fg, ts, np.linspace(fit.PSI_MIN, fit.PSI_MAX, N_GRID))
    t_sc = time.perf_counter() - t_sc
    order = sorted(range(len(scan)), key=lambda i: -scan[i]['f'])
    best = scan[order[0]]
    alt = next((scan[i] for i in order[1:] if abs(scan[i]['psi'] - best['psi']) >= np.log(10.0)), None)
    starts = [('scan_best', best['theta'])] + ([('scan_alt', alt['theta'])] if alt is not None else [('truth', t0)])
    rows, fits = [], []
    for name, st in starts:
        t1 = time.perf_counter()
        r = fit.maximize(fg, st, pg_tol=PG_TOL)
        fits.append(r)
        rows.append(fit_row(name, r, t1))
    top = max(fits, key=lambda r: r['f'])
    refits = 0
    cl0 = eng.lnl(fit.x_from_theta(t0))                                 # plan: cl(theta0) enters the refit rule
    best_ref = max(cl0, max(s['f'] for s in scan))
    while (not top['converged'] or top['f'] < best_ref - 0.01) and refits < 2:
        refits += 1
        st = (t0, 'refit_truth') if refits == 1 else (best['theta'], 'refit_scan_best')
        t1 = time.perf_counter()
        r = fit.maximize(fg, st[0], pg_tol=PG_TOL)
        rows.append(fit_row(st[1], r, t1))
        if r['f'] > top['f'] + 1e-9 or (r['converged'] and not top['converged'] and r['f'] >= top['f'] - 0.01):
            top = r
    failed = (not top['converged']) or top['f'] < best_ref - 0.01
    outcome = 'optimizer_failure' if failed else fit.classify(top['theta'], tau0, scan)
    pilot = pilot_extras(eng, fg, top, ts, t0) if os.environ.get('PSIGC_PILOT') == '1' else None
    sel = {'start': 'selected', 'lnL': top['f'], **theta_row(top['theta']), 'converged': top['converged'],
           'max_pg': top['max_pg'], 'max_pg_se': top.get('max_pg_se'), 'newton_dec': top.get('newton_dec'),
           'max_df3': top.get('max_df3'), 'cl_theta0': cl0,
           'refits': refits, 'outcome': outcome, 'scan_s': t_sc,
           'scan': [{'psi': s['psi'], 'f': s['f'], 'tau': float(np.exp(s['theta'][fit.I_LOGTAU]))} for s in scan],
           'is_fit': {'lnL_IS': s1['f'], 'converged': s1['converged'], **theta_row(ts)}}
    if pilot is not None:
        rows += pilot.pop('rows')
        sel['pilot'] = pilot
    return rows + [sel]


COLD_STARTS = (30.0, 200.0, 500.0)      # Guess.py initial tract lengths


def _phi_fg(fg):
    """(log eta, log p) parameterisation of coordinates 6, 7 (G5 sensitivity refit), exact gradients."""
    def to_t(z):
        t = np.array(z, float); t[fit.I_LOGTAU] = z[6] - z[7]; t[fit.I_PSI] = -z[7]; return t

    def f(z):
        v, g = fg(to_t(z))
        gz = g.copy(); gz[6] = g[fit.I_LOGTAU]; gz[7] = -g[fit.I_LOGTAU] - g[fit.I_PSI]
        return v, gz
    return f, to_t


def legacy_emulation(eng, ts, L0):
    """Guess.py settings: legacy x layout (log eta, log p), start = IS fit with r2, r3 = 0.7, 3.0 and
    1/p = L0, L-BFGS-B with the legacy bounds and ftol = 1e-9 (stringent 'high'), SciPy defaults otherwise."""
    from scipy.optimize import minimize
    x = fit.x_from_theta(ts)
    x[4], x[5] = np.log(0.7), np.log(3.0)
    x[7] = -np.log(L0); x[6] = ts[fit.I_LOGTAU] + x[7]
    bnds = [(None, -0.001)] * 3 + [(None, None)] * 4 + [(None, 0.0)] + [(None, None)] * (len(x) - 8)

    def obj(z):
        f, g = eng.lnl_grad(z)
        return (-f, -g) if np.isfinite(f) else (1e300, np.zeros_like(g))
    r = minimize(obj, x, jac=True, method='L-BFGS-B', bounds=bnds, options={'ftol': 1e-9})
    t = fit.theta_from_x(r.x)
    return {'L0': L0, 'lnL': float(-r.fun), 'psi': float(t[fit.I_PSI]), 'nit': int(r.nit), 'nfev': int(r.nfev),
            'message': str(r.message)}


def pilot_extras(eng, fg, top, ts, t0):
    """G5 pilot (plan section 7): 4 extra starts (theta0 and the Guess.py cold starts), optimizer-sensitivity
    refits from the M2 optimum, and the M2-legacy emulation (reps 1-5 of each cell)."""
    rows, fits = [], []
    starts = [('pilot_truth', t0)]
    for L0 in COLD_STARTS:
        st = ts.copy(); st[fit.I_PSI] = np.log(L0); starts.append((f'pilot_cold_{L0:g}', st))
    for name, st in starts:
        t1 = time.perf_counter()
        r = fit.maximize(fg, st, pg_tol=PG_TOL)
        fits.append(r); rows.append(fit_row(name, r, t1))
    best6 = max([top['f']] + [r['f'] for r in fits])
    th = top['theta']
    sens = {}
    t1 = time.perf_counter()
    sens['tol_x100'] = fit.maximize(fg, th, pg_tol=PG_TOL / 100, df_tol=1e-8, nd_tol=fit.ND_TOL / 100)
    sens['psi_bound_1e8'] = fit.maximize(fg, th, bounds=fit.default_bounds(len(th), psi=(fit.PSI_MIN, np.log(1e8))),
                                         pg_tol=PG_TOL)
    sens['maxcor_30'] = fit.maximize(fg, th, pg_tol=PG_TOL, maxcor=30)
    fphi, to_t = _phi_fg(fg)
    z = th.copy(); z[6] = th[fit.I_LOGTAU] - th[fit.I_PSI]; z[7] = -th[fit.I_PSI]
    bph = fit.default_bounds(len(th)); bph[7] = (-fit.PSI_MAX, -fit.PSI_MIN); bph[6] = (None, None)
    rph = fit.maximize(fphi, z, bounds=bph, pg_tol=PG_TOL)
    rph['theta'] = to_t(rph['theta'])
    sens['log_eta_log_p'] = rph
    sens_rows = {k: {'lnL': float(r['f']), 'd_lnL': float(r['f'] - top['f']), 'psi': float(r['theta'][fit.I_PSI]),
                     'converged': r['converged'], 'nfev': r['nfev']} for k, r in sens.items()}
    out = {'rows': rows, 'best6_lnL': float(best6), 'm2_minus_best6': float(top['f'] - best6),
           'sensitivity': sens_rows, 'sensitivity_s': time.perf_counter() - t1}
    if os.environ.get('PSIGC_LEGACY_EMULATION') == '1':
        t1 = time.perf_counter()
        em = [legacy_emulation(eng, ts, L0) for L0 in COLD_STARTS]
        b = max(em, key=lambda e: e['lnL'])
        out['legacy'] = {'starts': em, 'best_lnL': b['lnL'], 'best_psi': b['psi'],
                         'd_lnL': float(b['lnL'] - top['f']), 'd_psi': float(b['psi'] - th[fit.I_PSI]),
                         'wall_s': time.perf_counter() - t1}
    return out


ESTIMATORS = {'M0': est_M0, 'M1': est_M1, 'M2': est_M2}


def make_engine(ds, backend, device):
    """'cpu': psigc.PSEngine; 'gpu': psigc.gpu.GPUEngine (CUDA graphs, deterministic); both adaptive K."""
    if backend == 'cpu':
        return psigc.PSEngine(ds, K=K_NODES)
    from IGCexpansion.psigc.gpu import GPUEngine
    return GPUEngine(ds, K=K_NODES, device=device)


def eval_summary(eng):
    """Per-task summary of the engine's per-evaluation diagnostics (GPU engine: eng.eval_log)."""
    log = getattr(eng, 'eval_log', None)
    if not log:
        return None
    Ks = [e['K'] for e in log]
    return {'n_eval': len(log), 'K_hist': {str(k): Ks.count(k) for k in sorted(set(Ks))},
            'n_replay': int(sum(e.get('n_try', 1) - 1 for e in log)),
            'n_fallback': int(sum(bool(e['fallback']) for e in log)), 'tail_max': float(max(e['tail_max'] for e in log)),
            'n_nonpos': int(sum(e['nonpos'] for e in log)), 'engine_s': float(sum(e['time'] for e in log)),
            'n_capture': getattr(eng, 'n_capture', None), 'n_graph_clears': getattr(eng, 'n_graph_clears', None),
            # compact per-evaluation trace: K, tail/threshold (max over classes), attempts, exact-fallback flag
            'trace': {'K': Ks, 'tail': [float(f"{e['tail_max']:.3g}") for e in log],
                      'n_try': [e.get('n_try', 1) for e in log], 'fallback': [int(bool(e['fallback'])) for e in log],
                      'kind': ''.join('g' if e['kind'] == 'lnl_grad' else 'f' for e in log)}}


def task_path(out, est, m):
    return os.path.join(out, est, f'tract_{m["cell"]:g}', f'sim_{m["rep"]}.json')


_WORKER = {}


def init_worker(backend, devq, numa_cpus):
    """Pool initializer: GPU workers take a device from the queue and pin themselves to its NUMA node."""
    if backend == 'gpu':
        dev = devq.get()
        _WORKER['device'] = f'cuda:{dev}'
        cpus = numa_cpus.get(str(dev))
        if cpus and hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(0, cpus)
        import torch
        torch.cuda.set_device(_WORKER['device'])


def run_task(args):
    bundle, out, est, m, chash, backend, device = args
    device = _WORKER.get('device', device)
    path = task_path(out, est, m)
    if os.path.exists(path):
        return est, m['cell'], m['rep'], 'skip', 0.0
    t_start = time.perf_counter()
    try:
        tree, ds = load(bundle, m)
        eng, is_eng = make_engine(ds, backend, device), psigc.ISEngine(ds)
        t0 = fit.theta_from_x(np.array(m['truth_x']))
        tau0 = float(np.exp(t0[fit.I_LOGTAU]))
        if os.environ.get('PSIGC_PILOT') == '1':
            os.environ['PSIGC_LEGACY_EMULATION'] = '1' if m['rep'] <= 5 else '0'
            t0[fit.I_PSI] = max(t0[fit.I_PSI], fit.PSI_MIN)          # 1/p = 1 cell: truth on the bound
        with np.errstate(all='ignore'):
            rows = ESTIMATORS[est](eng, is_eng, t0, tau0)
        wall = time.perf_counter() - t_start
        meta = {'estimator': est, 'cell': m['cell'], 'rep': m['rep'], 'sha256': m['sha256'],
                'excluded': m['excluded'], 'N': m['N'], 'N_inf': m['N_inf'], 'm_real': m['m_real'],
                'truth_psi': float(np.log(m['cell'])), 'truth_tau': tau0, 'K': K_NODES,
                'backend': backend, 'device': device if backend == 'gpu' else None, 'pg_tol': PG_TOL, 'conv_rule': 'newton', 'nd_tol': fit.ND_TOL,
                'psi_bounds': [fit.PSI_MIN, fit.PSI_MAX], 'evals': eval_summary(eng),
                'versions': provenance(backend, device),
                'code_hash': chash, 'host': os.uname().nodename, 'task_wall_s': wall,
                'finished': time.strftime('%Y-%m-%dT%H:%M:%S')}
        for r in rows:
            r.setdefault('estimator', est)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + f'.tmp{os.getpid()}'
        with open(tmp, 'w') as f:
            json.dump({'meta': meta, 'rows': rows}, f)
        os.replace(tmp, path)
        return est, m['cell'], m['rep'], 'ok', wall
    except Exception:
        os.makedirs(os.path.join(out, 'errors'), exist_ok=True)
        with open(os.path.join(out, 'errors', f'{est}_tract_{m["cell"]:g}_sim_{m["rep"]}.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return est, m['cell'], m['rep'], 'error', time.perf_counter() - t_start


def read_manifest(bundle):
    return [json.loads(l) for l in open(os.path.join(bundle, 'manifest.jsonl'))]


def select(a):
    man = read_manifest(a.bundle)
    if a.cells:
        cells = {float(c) for c in a.cells.split(',')}
        man = [m for m in man if m['cell'] in cells]
    if a.reps:
        lo, hi = (int(v) for v in a.reps.split('-'))
        man = [m for m in man if lo <= m['rep'] <= hi]
    ests = a.estimators.split(',')
    # cheap estimators first, then M2; within an estimator interleave cells so partial results span the grid
    return [(e, m) for e in ests for m in sorted(man, key=lambda m: (m['rep'], m['cell']))]


def cmd_run(a):
    from multiprocessing import get_context
    if getattr(a, 'pilot', False):
        os.environ['PSIGC_PILOT'] = '1'                  # inherited by the spawned workers
    tasks = [(e, m) for e, m in select(a) if not os.path.exists(task_path(a.out, e, m))]
    chash = code_hash()
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, 'runner.log'), 'a')
    log.write(f'{time.strftime("%F %T")} start {len(tasks)} tasks, workers {a.workers}, code {chash}\n'); log.flush()
    done = 0
    if a.backend == 'gpu':
        ctx = get_context('spawn')
        devs = [d for d in a.devices.split(',')]
        devq = ctx.Queue()
        for d in devs:
            for _ in range(a.workers_per_gpu):
                devq.put(d)
        nw = len(devs) * a.workers_per_gpu
        init = (init_worker, (a.backend, devq, numa_cpus(devs)))
        maxt = None                      # keep graph caches and the CUDA context for the life of a worker
    else:
        ctx, nw, init, maxt = get_context('fork'), a.workers, (None, ()), 20
    with ctx.Pool(nw, initializer=init[0], initargs=init[1], maxtasksperchild=maxt) as pool:
        for est, g, k, status, wall in pool.imap_unordered(
                run_task, [(a.bundle, a.out, e, m, chash, a.backend, 'cuda:0') for e, m in tasks]):
            done += 1
            log.write(f'{time.strftime("%F %T")} {done}/{len(tasks)} {est} tract {g:g} sim {k} {status} {wall:.0f}s\n')
            log.flush()
    log.write(f'{time.strftime("%F %T")} finished\n')


def numa_cpus(devs):
    """CPU list of each GPU's NUMA node (from nvidia-smi topo); empty if unavailable."""
    out = {}
    try:
        import subprocess
        topo = subprocess.run(['nvidia-smi', 'topo', '-m'], capture_output=True, text=True, timeout=20).stdout
        for line in topo.splitlines():
            parts = line.split()
            if parts and parts[0].startswith('GPU') and parts[0][3:].isdigit() and parts[0][3:] in devs:
                aff = next((p for p in parts if '-' in p and p.replace('-', '').replace(',', '').isdigit()), None)
                if aff:
                    cpus = set()
                    for rng in aff.split(','):
                        lo, hi = (int(v) for v in rng.split('-'))
                        cpus.update(range(lo, hi + 1))
                    out[parts[0][3:]] = cpus
    except Exception:
        pass
    return out


def cmd_status(a):
    a.cells = a.reps = None
    a.estimators = 'M0,M1,M2'
    tot, have = {}, {}
    for e, m in select(a):
        tot[e] = tot.get(e, 0) + 1
        have[e] = have.get(e, 0) + os.path.exists(task_path(a.out, e, m))
    err = len(os.listdir(os.path.join(a.out, 'errors'))) if os.path.isdir(os.path.join(a.out, 'errors')) else 0
    print({e: f'{have[e]}/{tot[e]}' for e in tot}, 'errors:', err)


def cmd_collect(a):
    import csv
    rows = []
    for root, _, files in os.walk(a.out):
        for fn in files:
            if fn.endswith('.json'):
                d = json.load(open(os.path.join(root, fn)))
                for r in d['rows']:
                    rr = {**d['meta'], **{k: v for k, v in r.items() if k not in ('grid', 'scan', 'theta', 'is_fit')}}
                    rr['theta'] = json.dumps(r.get('theta'))
                    rr['estimator'] = r.get('estimator', d['meta']['estimator'])
                    rows.append(rr)
    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(a.out, 'rows.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    try:
        import pandas as pd
        pd.DataFrame(rows).to_parquet(os.path.join(a.out, 'rows.parquet'))
    except Exception as e:
        print('parquet skipped:', e)
    print(len(rows), 'rows')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('bundle'); b.add_argument('--src', required=True); b.add_argument('--out', required=True)
    r = sub.add_parser('run')
    for s in (r, sub.add_parser('status')):
        s.add_argument('--bundle', required=True); s.add_argument('--out', required=True)
    r.add_argument('--workers', type=int, default=8); r.add_argument('--estimators', default='M0,M1,M2')
    r.add_argument('--cells'); r.add_argument('--reps')
    r.add_argument('--backend', choices=['cpu', 'gpu'], default='cpu')
    r.add_argument('--devices', default='0,1'); r.add_argument('--workers-per-gpu', type=int, default=1)
    r.add_argument('--pilot', action='store_true', help='G5 pilot extras (6 starts, sensitivity refits, legacy emulation)')
    c = sub.add_parser('collect'); c.add_argument('--out', required=True)
    a = p.parse_args()
    {'bundle': cmd_bundle, 'run': cmd_run, 'status': cmd_status, 'collect': cmd_collect}[a.cmd](a)
