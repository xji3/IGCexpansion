"""Steps 8-9 for fitted datasets: M3 profile, outcome class, sandwich UQ (J_PB, J_HAC, legacy Godambe), CLR
intervals and the p = 1 test, one JSON per dataset (psigc.infer.infer_dataset).

  run    --bundle <armL bundle> --m2 <dir with M2/tract_g/sim_k.json> --out <dir> [--cells 3,10] [--reps 1-5]
         --devices 0,1 [--workers-per-gpu 1] [--sim-workers 24] [--R 100] [--arm 0]
  pb     same options as run (no --local-he/--h-reps): parametric bootstrap with full refits of the J_PB
         replicates (armL.est_M2), percentile interval for psi; one JSON per dataset in <out>/PB/
  status --bundle <dir> --out <dir>

GPU only; the devices are the ones given (no selection or switching). Each GPU worker keeps a CPU process
pool that simulates and builds the J_PB replicates while the GPU evaluates. Idempotent like armL.py.
"""
import argparse, json, os, sys, time, traceback

for v in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import armL                                            # noqa: E402
from IGCexpansion.psigc import fit                                  # noqa: E402
from IGCexpansion.psigc.seeds import STAGE, STUDY_ENTROPY           # noqa: E402

_W = {}


def worker_main(dev, numa, sim_workers, task_q, res_q):
    """One GPU worker: a plain process (not an executor worker) that owns its simulation pool and shuts it
    down itself; nested executors hang at exit because multiprocessing joins children before cleanup."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    _W['device'] = f'cuda:{dev}'
    cpus = numa.get(str(dev))
    if cpus and hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, cpus)
    import torch
    torch.cuda.set_device(_W['device'])
    _W['pool'] = ProcessPoolExecutor(sim_workers, mp_context=get_context('spawn'))
    try:
        while True:
            args = task_q.get()
            if args is None:
                break
            res_q.put(run_task(args))
    finally:
        _W['pool'].shutdown(wait=True, cancel_futures=True)


def out_path(out, m):
    return os.path.join(out, 'UQ', f'tract_{m["cell"]:g}', f'sim_{m["rep"]}.json')


def pb_path(out, m):
    return os.path.join(out, 'PB', f'tract_{m["cell"]:g}', f'sim_{m["rep"]}.json')


def run_pb_task(args):
    """Parametric bootstrap with full refits for one arm-L dataset: the J_PB replicates at the M2 estimate
    (same seed key as uq_run's J_PB, so the same R datasets), each fitted by armL.est_M2 with theta_hat in the
    place of the truth; percentile interval for psi (uq.pb_percentile)."""
    _, bundle, m2dir, out, m, R, arm, chash, ci = args
    path = pb_path(out, m)
    if os.path.exists(path):
        return m['cell'], m['rep'], 'skip', 0.0
    t0 = time.perf_counter()
    fit._free_gpu_cache()
    try:
        from IGCexpansion import psigc
        from IGCexpansion.psigc import uq
        from IGCexpansion.psigc.gpu import GPUEngine
        dev = _W['device']
        J = json.load(open(armL.task_path(m2dir, 'M2', m)))
        sel = [r for r in J['rows'] if r.get('start') == 'selected'][0]
        t_hat = np.array(sel['theta'])
        tau_hat = float(np.exp(t_hat[fit.I_LOGTAU]))
        si = np.loadtxt(os.path.join(bundle, 'seq_index.txt'), dtype=np.int64)
        newick = open(os.path.join(bundle, 'YeastTree.newick')).read()

        def refit(eng, ds):
            rows = armL.est_M2(eng, psigc.ISEngine(ds), t_hat.copy(), tau_hat)
            s = rows[-1]
            return {'theta': s['theta'], 'f': s['lnL'], 'converged': s['converged'], 'outcome': s['outcome']}

        key = (STAGE['j_pb'], arm, ci, m['rep'])
        rows = uq.pb_refits(lambda d: GPUEngine(d, K=armL.K_NODES, device=dev, graphs=False), t_hat, newick, si, R,
                            STUDY_ENTROPY, key, refit, pool=_W['pool'])
        res = {'theta_hat': t_hat.tolist(), 'interval': uq.pb_percentile(rows), 'rows': rows}
        meta = {'cell': m['cell'], 'rep': m['rep'], 'truth_theta': fit.theta_from_x(np.array(m['truth_x'])).tolist(),
                'R': R, 'seed_key': list(key), 'entropy': str(STUDY_ENTROPY), 'device': dev, 'code_hash': chash,
                'versions': armL.provenance('gpu', dev), 'task_wall_s': time.perf_counter() - t0,
                'finished': time.strftime('%Y-%m-%dT%H:%M:%S')}
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + f'.tmp{os.getpid()}'
        with open(tmp, 'w') as f:
            json.dump({'meta': meta, 'result': res}, f)
        os.replace(tmp, path)
        return m['cell'], m['rep'], 'ok', time.perf_counter() - t0
    except Exception as e:
        os.makedirs(os.path.join(out, 'errors'), exist_ok=True)
        with open(os.path.join(out, 'errors', f'PB_tract_{m["cell"]:g}_sim_{m["rep"]}.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return m['cell'], m['rep'], ('oom' if fit._is_oom(e) else 'error'), time.perf_counter() - t0


def run_task(args):
    """One dataset; after an out-of-memory error the GPU cache is freed and the task is retried once (the
    caching allocator can fragment over a long run: 69 GB reserved but unallocated, armL_UQ2 tract 300 sim 14)."""
    one = run_pb_task if args[0] == 'pb' else _run_task
    r = one(args)
    if r[2] == 'oom':
        fit._free_gpu_cache()
        r = one(args)
        r = r[:2] + (('retried_' + r[2]) if r[2] in ('ok', 'error', 'oom') else r[2],) + r[3:]
    return r


def _run_task(args):
    bundle, m2dir, out, m, R, arm, chash, ci, local_he, h_reps = args
    path = out_path(out, m)
    if os.path.exists(path):
        return m['cell'], m['rep'], 'skip', 0.0
    t0 = time.perf_counter()
    fit._free_gpu_cache()
    try:
        from IGCexpansion.psigc.gpu import GPUEngine
        from IGCexpansion.psigc.infer import infer_dataset
        dev = _W['device']
        J = json.load(open(armL.task_path(m2dir, 'M2', m)))
        sel = [r for r in J['rows'] if r.get('start') == 'selected'][0]
        tree, ds = armL.load(bundle, m)
        si = np.loadtxt(os.path.join(bundle, 'seq_index.txt'), dtype=np.int64)
        newick = open(os.path.join(bundle, 'YeastTree.newick')).read()
        eng = GPUEngine(ds, K=armL.K_NODES, device=dev)
        eng.begin_task()
        t_truth = fit.theta_from_x(np.array(m['truth_x']))
        key = (STAGE['j_pb'], arm, ci, m['rep'])
        res = infer_dataset(eng, np.array(sel['theta']), sel['lnL'], float(np.exp(t_truth[fit.I_LOGTAU])),
                            newick, si, R=R, entropy=STUDY_ENTROPY, key=key,
                            make_engine=lambda d: GPUEngine(d, K=armL.K_NODES, device=dev, graphs=False),
                            pool=_W['pool'], local_he=local_he, h_reps=h_reps)
        res['m2_outcome'] = sel['outcome']
        meta = {'cell': m['cell'], 'rep': m['rep'], 'sha256': m['sha256'], 'excluded': m['excluded'],
                'truth_theta': t_truth.tolist(), 'R': R, 'seed_key': list(key), 'entropy': str(STUDY_ENTROPY),
                'device': dev, 'versions': armL.provenance('gpu', dev), 'code_hash': chash,
                'm2_code_hash': J['meta'].get('code_hash'), 'local_he': local_he, 'h_reps': h_reps, 'evals': armL.eval_summary(eng),
                'host': os.uname().nodename, 'task_wall_s': time.perf_counter() - t0,
                'finished': time.strftime('%Y-%m-%dT%H:%M:%S')}
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + f'.tmp{os.getpid()}'
        with open(tmp, 'w') as f:
            json.dump({'meta': meta, 'result': res}, f)
        os.replace(tmp, path)
        return m['cell'], m['rep'], 'ok', time.perf_counter() - t0
    except Exception as e:
        os.makedirs(os.path.join(out, 'errors'), exist_ok=True)
        with open(os.path.join(out, 'errors', f'UQ_tract_{m["cell"]:g}_sim_{m["rep"]}.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return m['cell'], m['rep'], ('oom' if fit._is_oom(e) else 'error'), time.perf_counter() - t0


def select(a):
    man = armL.read_manifest(a.bundle)
    if a.cells:
        cells = {float(c) for c in a.cells.split(',')}
        man = [m for m in man if m['cell'] in cells]
    if a.reps:
        lo, hi = (int(v) for v in a.reps.split('-'))
        man = [m for m in man if lo <= m['rep'] <= hi]
    return sorted(man, key=lambda m: (m['rep'], m['cell']))


def cmd_run(a):
    from multiprocessing import get_context
    pb = getattr(a, 'cmd', None) == 'pb'
    done_path = pb_path if pb else out_path
    tasks = [m for m in select(a) if not os.path.exists(done_path(a.out, m))
             and os.path.exists(armL.task_path(a.m2, 'M2', m))]
    chash = armL.code_hash()
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, 'runner.log'), 'a')
    log.write(f'{time.strftime("%F %T")} start {len(tasks)} tasks, devices {a.devices}, code {chash}\n'); log.flush()
    devs = a.devices.split(',')
    ctx = get_context('spawn')
    task_q, res_q = ctx.Queue(), ctx.Queue()
    numa = armL.numa_cpus(devs)
    procs = [ctx.Process(target=worker_main, args=(d, numa, a.sim_workers, task_q, res_q))
             for d in devs for _ in range(a.workers_per_gpu)]
    for p in procs:
        p.start()
    cells_all = sorted({m['cell'] for m in armL.read_manifest(a.bundle)})      # cell index for the seed key
    for m in tasks:
        ci = cells_all.index(m['cell'])
        task_q.put(('pb', a.bundle, a.m2, a.out, m, a.R, a.arm, chash, ci) if pb else
                   (a.bundle, a.m2, a.out, m, a.R, a.arm, chash, ci, a.local_he, a.h_reps))
    for _ in procs:
        task_q.put(None)
    done = 0
    import queue
    while done < len(tasks):
        try:
            g, k, status, wall = res_q.get(timeout=60)
        except queue.Empty:
            if not any(p.is_alive() for p in procs):
                log.write(f'{time.strftime("%F %T")} all workers exited with {len(tasks) - done} tasks left\n')
                break
            continue
        done += 1
        log.write(f'{time.strftime("%F %T")} {done}/{len(tasks)} {"PB" if pb else "UQ"} tract {g:g} sim {k} {status} {wall:.0f}s\n')
        log.flush()
    for p in procs:
        p.join()
    log.write(f'{time.strftime("%F %T")} finished\n')


def cmd_status(a):
    a.cells = a.reps = None
    man = select(a)
    have = sum(os.path.exists(out_path(a.out, m)) for m in man)
    err = len(os.listdir(os.path.join(a.out, 'errors'))) if os.path.isdir(os.path.join(a.out, 'errors')) else 0
    print(f'UQ {have}/{len(man)}', 'errors:', err)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run'); s = sub.add_parser('status')
    b = sub.add_parser('pb', help='parametric bootstrap with full refits (percentile interval for psi)')
    for q in (r, s, b):
        q.add_argument('--bundle', required=True); q.add_argument('--out', required=True)
    for q in (r, b):
        q.add_argument('--m2', required=True); q.add_argument('--cells'); q.add_argument('--reps')
        q.add_argument('--devices', required=True); q.add_argument('--workers-per-gpu', type=int, default=1)
        q.add_argument('--sim-workers', type=int, default=24); q.add_argument('--R', type=int, default=100)
        q.add_argument('--arm', type=int, default=0, help='0 = arm L (part of the J_PB seed key)')
    r.add_argument('--local-he', action='store_true',
                   help='primary CLR / p = 1 test by local lambda(psi) with the expected Hessian (psigc.infer)')
    r.add_argument('--h-reps', type=int, help='replicates with per-pair scores for H_E (default: all R)')
    a = p.parse_args()
    {'run': cmd_run, 'pb': cmd_run, 'status': cmd_status}[a.cmd](a)
