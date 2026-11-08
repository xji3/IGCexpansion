"""M2/M3 benchmark: the full PS-IGC objective on resident workers (plan M2, D1; M3 GL gradients).

    .venv/bin/python benchmarks/m2_psigc.py --procs 10,32,64 --engines beagle,jsonctmctree
    .venv/bin/python benchmarks/m2_psigc.py --procs 10,32,64 --fd beagle:gl,beagle:central --tag m3

Case: PSJSGeneconv HKY 'One rate' with rate variation on YDR418W_YEL054C (1,461 scenes), at
unpack_x(self.x), IGC on (force=None) so every a(n), b(n) term is live. For each engine and
process count: pool start-up, lnL (best of --repeats), lnL + gradient (best of --grad-repeats;
central differences, 7-8 free model parameters), worker time balance and peak RSS per worker.
The gradient is checked between engines. Writes benchmarks/results/m2_psigc_<host>.json.
"""

import argparse
import json
import os
import platform
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(os.path.dirname(HERE), "tests")]

import numpy as np  # noqa: E402

import cases  # noqa: E402
from IGCexpansion.beagle import psjs  # noqa: E402


def best(fn, k):
    out, t = None, float("inf")
    for _ in range(k):
        t0 = time.perf_counter()
        out = fn()
        t = min(t, time.perf_counter() - t0)
    return t, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", default="10")
    ap.add_argument("--engines", default="beagle,jsonctmctree")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--grad-repeats", type=int, default=2)
    ap.add_argument("--fd", default="", help="model-parameter gradient per engine, e.g. beagle:gl,beagle:central")
    ap.add_argument("--tag", default="m2")
    a = ap.parse_args()
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w_igc(d)
    x = np.array(obj.x, float)
    R = {"case": "PSJSGeneconv HKY One rate + rv, YDR418W_YEL054C, IGC on, at x", "x": x.tolist(),
         "env": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform(),
                 "cpu_count": os.cpu_count(), "date": time.strftime("%Y-%m-%d %H:%M"),
                 "beagle_prefix": os.environ.get("IGC_BEAGLE_PREFIX", "default")},
         "runs": []}
    grads = {}
    combos = [tuple(c.split(":")) for c in a.fd.split(",")] if a.fd else [(e, None) for e in a.engines.split(",")]
    for k in [int(v) for v in a.procs.split(",")]:
        for eng, fd in combos:
            p = psjs.PSJSPool(obj, engine=eng, processes=k, fd=fd)
            try:
                t_ll, ll = best(lambda: p.loglikelihood(x), a.repeats)
                bal_ll = p.last["worker_sec"]
                t_g, (llg, g) = best(lambda: p.loglikelihood_and_gradient(x), a.grad_repeats)
                run = {"engine": eng, "fd": p.fd, "procs": k, "gl_clipped": p.last.get("gl_clipped", 0), "scenes": len(p.scenes), "free_params": len(p.free),
                       "sec_start": p.sec_start, "sec_lnL": t_ll, "sec_lnL_grad": t_g, "lnL": ll,
                       "worker_sec_lnL_max_over_mean": float(np.max(bal_ll) / np.mean(bal_ll)),
                       "worker_sec_grad_max_over_mean": float(np.max(p.last["worker_sec"]) / np.mean(p.last["worker_sec"])),
                       "worker_rss_mb_start": p.worker_rss_mb, "worker_rss_mb_peak": p.last["worker_rss_mb"],
                       "lnL_calls_per_grad": len(p.scenes) * (1 + 2 * len(p.free))}
                grads[(eng, p.fd, k)] = g
                R["runs"].append(run)
                print(json.dumps({kk: run[kk] for kk in ("engine", "fd", "procs", "sec_start", "sec_lnL",
                                                         "sec_lnL_grad", "lnL")}), flush=True)
            finally:
                p.close()
    keys = list(grads)
    R["grad_max_abs_diff_between_runs"] = float(max(np.max(np.abs(grads[i] - grads[j])) for i in keys for j in keys))
    R["grad_max_abs"] = float(np.max(np.abs(grads[keys[0]])))
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    R["grads"] = {"%s:%s:%d" % key: g.tolist() for key, g in grads.items()}
    path = os.path.join(HERE, "results", "%s_psigc_%s.json" % (a.tag, platform.node().split(".")[0]))
    with open(path, "w") as f:
        json.dump(R, f, indent=1)
        f.write("\n")
    print(path, "grad diff", R["grad_max_abs_diff_between_runs"], "of", R["grad_max_abs"])


if __name__ == "__main__":
    main()
