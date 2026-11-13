"""Re-check the convergence of bootstrap refits labelled optimizer_failure with an exact Newton decrement:
regenerate replicate r (same seed key as cmd_pb), FD Hessian of the exact gradient at the saved refit theta, on
the parameters cl depends on (non-coding data: r2, r3 drop out). Writes one row per replicate.
  python gpu2/recheck_pb_nd.py <realdata root> <pb dir> <name> <device> <out.json> [workers]"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

if __name__ == '__main__':
    import torch
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    import realdata as rd
    from IGCexpansion.psigc import fit, uq
    from IGCexpansion.psigc.gpu import GPUEngine
    from IGCexpansion.psigc.seeds import STUDY_ENTROPY
    root, pbdir, name, dev, out = sys.argv[1:6]
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 16
    torch.cuda.set_device(dev)
    C = rd.configs(root)
    tree, ds0, si = rd.load(C[name])
    P = json.load(open(os.path.join(pbdir, f'{name}.pb.json')))
    kw, _ = rd.engine_kwargs(ds0)
    x_hat = fit.x_from_theta(np.array(P['theta_hat']))
    newick = open(C[name]['tree']).read()
    sil = list(map(np.asarray, si)) if isinstance(si, (list, tuple)) else np.asarray(si)
    rows = [q for q in P['rows'] if 'theta' in q]
    jobs = [(x_hat, newick, sil, STUDY_ENTROPY, tuple(P['seed_key']) + (q['r'],)) for q in rows]
    res = []
    with ProcessPoolExecutor(workers, mp_context=get_context('spawn')) as pool:
        for q, (ds, _) in zip(rows, pool.map(uq._sim_one, jobs)):
            t0 = time.time()
            eng = GPUEngine(ds, device=dev, **{**kw, 'graphs': False})
            fg = fit.theta_fg(eng)
            t = np.array(q['theta'])
            f, g = fg(t)
            pg = fit._projected_grad(t, g, fit.default_bounds(len(t)))
            H = -fit.hessian(fg, t)
            act = [i for i in range(len(t)) if np.abs(H[i]).max() > 1e-6]
            Ha = H[np.ix_(act, act)]
            w = np.linalg.eigvalsh(Ha)
            nd = float(np.sqrt(max(pg[act] @ np.linalg.solve(Ha, pg[act]), 0))) if w.min() > 0 else None
            row = {'r': q['r'], 'psi': q['psi'], 'outcome': q['outcome'], 'f_saved': q['f'], 'f_diff': float(f - q['f']),
                   'max_pg': float(np.abs(pg).max()), 'nd_fd': nd, 'min_eig': float(w.min()), 'n_active': len(act),
                   'wall_s': time.time() - t0}
            res.append(row)
            print(json.dumps(row), flush=True)
            del eng
            torch.cuda.empty_cache()
    json.dump({'name': name, 'rows': res}, open(out, 'w'))
    print('done', flush=True)
