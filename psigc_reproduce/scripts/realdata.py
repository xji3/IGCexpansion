"""Real-data configurations for the M2 refits and arm X3 (PSIGC_GPU_PLAN.md step 10; SIMULATION_REDO_PLAN.md
section 4, X3): the 14 filtered yeast pairs, EDN/ECP and HBG on primate trees, and the 20 Harpak introns
(joint, within-locus pairs only, non-coding).

  list   --root <IGCTract project dir>
  sizes  --root <dir> [--names a,b]                     columns, pairs, unique patterns (CPU)
  bench  --root <dir> --names a,b --device cuda:0       IS fit, then lnL + gradient time and peak GPU memory
  fit    --root <dir> --names a,b --device cuda:0 --out <dir> [--R 100]   M2 (Guess.py cold starts as
         alternatives) + profile/UQ/CLR (psigc.infer); one JSON per dataset

Large datasets (> MAX_PATTERNS unique pair patterns in a class) use pattern chunking (GPUEngine max_patterns)
and eager evaluation. The device is the one given; nothing selects or switches hardware.
"""
import argparse, json, os, sys, time

for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from IGCexpansion import psigc                                           # noqa: E402
from IGCexpansion.psigc import fit, uq                              # noqa: E402
from IGCexpansion.psigc.seeds import STAGE, STUDY_ENTROPY           # noqa: E402

MAX_PATTERNS = 20000          # per chunk; one chunk at K = 48 stays within a few GB
CHUNK_ABOVE = 100000          # datasets with more patterns than this in a class are chunked


def configs(root):
    T = os.path.join(root, 'IGCTract')
    out = {}
    for pair in [l.strip() for l in open(os.path.join(T, 'Filtered_pairs.txt')) if l.strip()]:
        d = os.path.join(T, 'MafftAlignment', pair)
        out[pair] = dict(tree=os.path.join(T, 'YeastTree.newick'), paralogs=tuple(pair.split('_')), outgroup_gene=None,
                         loci=[(os.path.join(d, f'{pair}_input.fasta'), os.path.join(d, f'{pair}_seq_index.txt'))])
    P = os.path.join(T, 'ParametricBootstrap')
    out['EDN_ECP'] = dict(tree=os.path.join(P, 'input_tree.newick'), paralogs=('ECP', 'EDN'), outgroup_gene='EDN',
                          loci=[(os.path.join(P, 'EDN_ECP_Cleaned_input.fasta'), os.path.join(P, 'EDN_ECP_seq_index.txt'))])
    H = os.path.join(T, 'Harpak')
    names = [l.strip().replace('.pos.seq.formatted', '').replace('.', '_')
             for l in open(os.path.join(H, 'missing_0_species_list.txt')) if l.strip()]
    out['Harpak_introns'] = dict(tree=os.path.join(H, 'HarpakTree.newick'), paralogs=('Paralog1', 'Paralog2'),
                                 outgroup_gene='Paralog1',
                                 loci=[(os.path.join(H, 'prepared_input', n + '.fasta'),
                                        os.path.join(H, 'prepared_input', n + '_seq_index.txt')) for n in names])
    G = os.path.join(root, 'Datasets', 'Human1_HBG', 'model_input')
    for kind, stem in (('genomic', 'HBG_genomic'), ('cds', 'HBG_CDS')):
        out[f'HBG_{kind}'] = dict(tree=os.path.join(G, kind, 'HBG_tree.newick'), paralogs=('HBG2', 'HBG1'),
                                  outgroup_gene='HBG2',
                                  loci=[(os.path.join(G, kind, f'{stem}.fasta'), os.path.join(G, kind, f'{stem}_seq_index.txt'))])
    return out


def load(cfg):
    tree = psigc.Tree(open(cfg['tree']).read())
    dss = [psigc.Dataset.from_fasta(fa, si, tree, cfg['paralogs'], cfg['outgroup_gene']) for fa, si in cfg['loci']]
    ds = dss[0] if len(dss) == 1 else psigc.Dataset.concat(dss)
    sis = [np.loadtxt(si, dtype=np.int64) for _, si in cfg['loci']]
    return tree, ds, (sis[0] if len(sis) == 1 else sis)


LONG_N = 1000                 # loci with near-pair distances beyond this get a longer K ladder
LONG_LADDER = (16, 24, 32, 48, 64, 96)
LONG_MAX_PATTERNS = 10000
PARTS_EDGES = 12              # trees with more post-duplication edges than this are evaluated class by class


def engine_kwargs(ds):
    """Engine options by dataset size. Long loci (HBG genomic: 4,277 sites) can need more Chebyshev nodes than
    48 at moderate 1/p; the exact fallback (one node per distinct distance, ~1,000) is far costlier than
    K = 64-96, so long loci get a longer ladder and smaller chunks (memory ~ K x chunk)."""
    from IGCexpansion.psigc.gpu import K_LADDER
    big = max(b.U for b in ds.blocks.values()) > CHUNK_ABOVE
    long_ = ds.n_max > LONG_N
    kw = dict(K='auto', K_ladder=LONG_LADDER if long_ else K_LADDER)
    if long_ or len(ds.tree.post_edges) > PARTS_EDGES:
        # one part per class: the saved squarings of all classes at once grow with tau (HBG_cds, 23 post edges,
        # K = 32: 22 / 34 / 41 GB at tau x 1 / 30 / 300 vs 3.7 / 5.7 / 6.9 GB in parts, same time;
        # gpu2/mem_probe.py), which ran out of memory on an 80 GB A100 in the HBG_cds profile
        kw['parts'] = True
    if big:
        kw.update(max_patterns=LONG_MAX_PATTERNS if long_ else MAX_PATTERNS, graphs=False)
    return kw, big


def engine(ds, device):
    from IGCexpansion.psigc.gpu import GPUEngine
    kw, big = engine_kwargs(ds)
    return GPUEngine(ds, device=device, **kw), big


def start_theta(ds, tree):
    """IS-IGC fit from a generic start (legacy Guess.py style values), psi then set per start."""
    x = np.r_[np.log([0.3, 0.4, 0.5]), np.log(4.0), 0.0, 0.0, np.log(5.0) - np.log(30.0), -np.log(30.0),
              np.log(0.1) * np.ones(tree.n_rate)]
    if not ds.noncoding:
        x[4], x[5] = np.log(0.7), np.log(3.0)
    t = fit.theta_from_x(x)
    fixed = [4, 5] if ds.noncoding else []
    is_eng = psigc.ISEngine(ds)
    free = [i for i in range(len(t)) if i != fit.I_PSI and i not in fixed]
    r = fit.maximize(fit.theta_fg(is_eng), t, free=free)
    return r['theta'], r


def cmd_list(a):
    for k, c in configs(a.root).items():
        print(k, len(c['loci']), 'loci', c['tree'])


def _names(a):
    C = configs(a.root)
    return C, (a.names.split(',') if a.names else list(C))


def cmd_sizes(a):
    C, names = _names(a)
    for n in names:
        t0 = time.time()
        tree, ds, _ = load(C[n])
        print(n, ds.summary(), 'loci', ds.n_loci, 'noncoding', ds.noncoding, 'post edges', len(tree.post_edges),
              f'{time.time() - t0:.0f}s', flush=True)


def cmd_bench(a):
    import torch
    C, names = _names(a)
    torch.cuda.set_device(torch.device(a.device))
    for n in names:
        tree, ds, _ = load(C[n])
        t_is, r_is = start_theta(ds, tree)
        torch.cuda.reset_peak_memory_stats()
        t0 = time.time()
        eng, big = engine(ds, a.device)
        t_build = time.time() - t0
        out = {'name': n, 'chunked': big, 'n_classes': len(eng.classes), 'build_s': t_build, 'IS_lnL': r_is['f']}
        for L0 in (30.0, 300.0):
            t = t_is.copy(); t[fit.I_PSI] = np.log(L0)
            x = fit.x_from_theta(t)
            eng.lnl_grad(x)                                  # warm-up / K settling
            torch.cuda.synchronize()
            t1 = time.time(); f, g = eng.lnl_grad(x); torch.cuda.synchronize()
            out[f'L0_{L0:g}'] = {'lnL': f, 'grad_s': time.time() - t1, 'K': eng.last['K'],
                                 'n_fallback': len(eng.last['fallback']), 'n_classes': len(eng.classes),
                                 'tail_max': float(max(eng.last.get('tails', {0: 0}).values()))}
            print(json.dumps({'name': n, 'L0': L0, **out[f'L0_{L0:g}']}), flush=True)
        out['peak_gb'] = torch.cuda.max_memory_allocated() / 1e9
        print(json.dumps(out), flush=True)


def cmd_fit(a):
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    C, names = _names(a)
    os.makedirs(a.out, exist_ok=True)
    pool = ProcessPoolExecutor(a.sim_workers, mp_context=get_context('spawn'))
    for ni, n in enumerate(sorted(C)):
        if n not in names:
            continue
        path = os.path.join(a.out, f'{n}.json')
        if os.path.exists(path):
            continue
        T0 = time.time()
        try:
            _fit_one(a, C, n, ni, path, pool, T0)
        except Exception as e:                               # one failed dataset does not stop the batch
            import traceback
            import torch
            json.dump({'name': n, 'error': repr(e)[:2000], 'traceback': traceback.format_exc()[-4000:],
                       'wall_s': time.time() - T0, 'device': a.device},
                      open(os.path.join(a.out, f'{n}.error.json'), 'w'))
            print(n, 'FAILED', repr(e)[:200], flush=True)
            torch.cuda.empty_cache()
    pool.shutdown()


def _fit_one(a, C, n, ni, path, pool, T0):
    from IGCexpansion.psigc.gpu import GPUEngine
    from IGCexpansion.psigc.infer import infer_dataset
    tree, ds, si = load(C[n])
    t_is, r_is = start_theta(ds, tree)
    eng, big = engine(ds, a.device)
    alts = []
    for L0 in (30.0, 200.0, 500.0):                      # Guess.py cold starts
        t = t_is.copy(); t[fit.I_PSI] = np.log(L0); alts.append(t)
    r = fit.fit_m2(eng, t_is, tau0=float(np.exp(t_is[fit.I_LOGTAU])), alt_starts=alts,
                   fixed=uq.fixed_coords(ds))
    top = r['fit']
    kw, _ = engine_kwargs(ds)
    res = infer_dataset(eng, top['theta'], top['f'], float(np.exp(t_is[fit.I_LOGTAU])), open(C[n]['tree']).read(), si,
                        R=a.R, entropy=STUDY_ENTROPY, key=(STAGE['j_pb'], 9, ni, 0),
                        make_engine=lambda d: GPUEngine(d, device=a.device, **{**kw, 'graphs': False}),
                        pool=pool, local_he=a.local_he, h_reps=a.h_reps)
    json.dump({'name': n, 'summary': ds.summary(), 'n_loci': ds.n_loci, 'noncoding': ds.noncoding,
               'is_fit': {'lnL': r_is['f'], 'theta': t_is.tolist()},
               'm2': {'lnL': top['f'], 'theta': top['theta'].tolist(), 'outcome': r['outcome'], 'refits': r['refits'],
                      'converged': top['converged'], 'newton_dec': top.get('newton_dec'),
                      'scan': [{'psi': s['psi'], 'f': s['f']} for s in r['scan']],
                      'all_fits': [{'f': q['f'], 'psi': float(q['theta'][fit.I_PSI])} for q in r['all_fits']]},
               'infer': res, 'wall_s': time.time() - T0, 'device': a.device, 'local_he': a.local_he},
              open(path, 'w'))
    print(n, 'done', f'{time.time() - T0:.0f}s', flush=True)


def m2_refit(eng, ds):
    """The real-data M2 procedure of cmd_fit on any dataset: IS fit from the generic start, then fit_m2 with
    the Guess.py cold starts."""
    t_is, _ = start_theta(ds, ds.tree)
    alts = []
    for L0 in (30.0, 200.0, 500.0):
        t = t_is.copy(); t[fit.I_PSI] = np.log(L0); alts.append(t)
    r = fit.fit_m2(eng, t_is, tau0=float(np.exp(t_is[fit.I_LOGTAU])), alt_starts=alts,
                   fixed=uq.fixed_coords(ds))
    return {'theta': r['fit']['theta'], 'f': r['fit']['f'], 'converged': r['fit']['converged'], 'outcome': r['outcome']}


def cmd_pb(a):
    """Parametric bootstrap with full refits: the J_PB replicates of cmd_fit (same seed key) at the M2
    estimate in <fits>/<name>.json, each refitted by m2_refit; percentile interval for psi (uq.pb_percentile)."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    import torch
    from IGCexpansion.psigc import uq
    from IGCexpansion.psigc.gpu import GPUEngine
    C, names = _names(a)
    os.makedirs(a.out, exist_ok=True)
    torch.cuda.set_device(torch.device(a.device))
    pool = ProcessPoolExecutor(a.sim_workers, mp_context=get_context('spawn'))
    for ni, n in enumerate(sorted(C)):
        if n not in names:
            continue
        path = os.path.join(a.out, f'{n}.pb.json')
        src = os.path.join(a.fits, f'{n}.json')
        if os.path.exists(path) or not os.path.exists(src):
            continue
        T0 = time.time()
        try:
            tree, ds, si = load(C[n])
            t_hat = np.array(json.load(open(src))['m2']['theta'])
            kw, _ = engine_kwargs(ds)
            key = (STAGE['j_pb'], 9, ni, 0)
            log = lambda q: print(n, json.dumps({k: q.get(k) for k in ('r', 'psi', 'f', 'outcome', 'error', 'wall_s')}),
                                  flush=True)
            rows = uq.pb_refits(lambda d: GPUEngine(d, device=a.device, **{**kw, 'graphs': False}), t_hat,
                                open(C[n]['tree']).read(), si, a.R, STUDY_ENTROPY, key, m2_refit, pool=pool, log=log)
            json.dump({'name': n, 'theta_hat': t_hat.tolist(), 'interval': uq.pb_percentile(rows), 'rows': rows,
                       'R': a.R, 'seed_key': list(key), 'wall_s': time.time() - T0, 'device': a.device},
                      open(path, 'w'))
            print(n, 'pb done', f'{time.time() - T0:.0f}s', flush=True)
        except Exception as e:
            import traceback
            json.dump({'name': n, 'error': repr(e)[:2000], 'traceback': traceback.format_exc()[-4000:]},
                      open(os.path.join(a.out, f'{n}.pb.error.json'), 'w'))
            print(n, 'pb FAILED', repr(e)[:200], flush=True)
            torch.cuda.empty_cache()
    pool.shutdown()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    for nm in ('list', 'sizes', 'bench', 'fit', 'pb'):
        q = sub.add_parser(nm); q.add_argument('--root', required=True); q.add_argument('--names')
        if nm in ('bench', 'fit', 'pb'):
            q.add_argument('--device', required=True)
        if nm == 'pb':
            q.add_argument('--fits', required=True, help='dir with <name>.json from fit (theta_hat)')
            q.add_argument('--out', required=True); q.add_argument('--R', type=int, default=100)
            q.add_argument('--sim-workers', type=int, default=24)
        if nm == 'fit':
            q.add_argument('--out', required=True); q.add_argument('--R', type=int, default=100)
            q.add_argument('--sim-workers', type=int, default=24)
            q.add_argument('--local-he', action='store_true', help='primary CLR by local lambda with H_E')
            q.add_argument('--h-reps', type=int, help='replicates with per-pair scores for H_E (default: all R)')
    a = p.parse_args()
    {'list': cmd_list, 'sizes': cmd_sizes, 'bench': cmd_bench, 'fit': cmd_fit, 'pb': cmd_pb}[a.cmd](a)
