"""M4 benchmark: IS-IGC MG94 lnL + exact gradient on BEAGLE (plan M4, D3).

    .venv/bin/python benchmarks/m4_isigc.py --configs 1x1,1x8,8x1,8x4,16x4 [--gpu]

Case: ReCodonGeneconv MG94 on YGL062W_YBR218C (1129 codons, 1013 patterns, 3721 states,
6 process parameters + 12 edges), the M0 baseline case. A config ``PxT`` is P processes with
T BEAGLE threads each: P = 1 evaluates in this process, P > 1 uses a RecodonPool that splits the
GL edges over the workers. ``--gpu`` adds one in-process run on IGC_BEAGLE_RESOURCE=gpu.
For each config: start-up, lnL alone and lnL + full gradient (best of --repeats), and the
gradient's largest relative difference from the first config. The M0 baseline (jsonctmctree,
central FD over 8 processes) is benchmarks/results/m0_isigc_<host>.json.
Writes benchmarks/results/m4_isigc_<host>.json.
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
from IGCexpansion.beagle import recodon  # noqa: E402


def best(fn, k):
    out, t = None, float("inf")
    for _ in range(k):
        t0 = time.perf_counter()
        out = fn()
        t = min(t, time.perf_counter() - t0)
    return t, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="1x1,1x8,8x1,8x4")
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--tag", default="m4")
    a = ap.parse_args()
    with cases.quiet():
        obj = cases.mg94_ygl062w(tempfile.mkdtemp())
    x = np.array(obj.x, float)
    R = {"case": "ReCodonGeneconv MG94 YGL062W_YBR218C, 3721 states, exact GL gradient", "x": x.tolist(),
         "env": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform(),
                 "cpu_count": os.cpu_count(), "date": time.strftime("%Y-%m-%d %H:%M"),
                 "beagle_prefix": os.environ.get("IGC_BEAGLE_PREFIX", "default")},
         "runs": []}
    ref = None
    configs = [tuple(int(v) for v in c.split("x")) for c in a.configs.split(",")]
    if a.gpu:
        configs.append(("gpu", 1))
    for P, T in configs:
        row = {"processes": P, "threads": T, "loadavg_before": os.getloadavg()[0]}
        t0 = time.perf_counter()
        if P == "gpu" or P == 1:
            res = "gpu" if P == "gpu" else "cpu"
            recodon.detach(obj)
            obj.update_by_x(x)
            recodon.evaluate(obj, grad=False, threads=T, resource=res)
            row["sec_start"] = time.perf_counter() - t0
            row["sec_lnl"], (ll, _, _) = best(lambda: recodon.evaluate(obj, grad=False), a.repeats)
            row["sec_grad"], (ll, g, _) = best(lambda: recodon.evaluate(obj), a.repeats)
            st = recodon._STATE[obj]
            row["gl_nodes"] = int(sum(v[0] for v in st.bs._gl_nodes if v))
            row["gl_clipped"] = st.bs.gl_clipped
            recodon.detach(obj)
        else:
            pool = recodon.RecodonPool([obj], processes=P, threads=T)
            row["sec_start"] = time.perf_counter() - t0
            row["units"] = [[(u[0], u[1], u[2]) for u in us] for us in pool.units]
            row["sec_lnl"], (lls, _) = best(lambda: pool.evaluate([x], grad=False), a.repeats)
            row["sec_grad"], (lls, gs) = best(lambda: pool.evaluate([x]), a.repeats)
            row["worker_sec"] = pool.last["worker_sec"]
            row["worker_rss_mb"] = pool.worker_rss_mb
            ll, g = float(lls[0]), gs[0]
            pool.close()
        row["lnL"] = ll
        if ref is None:
            ref = g
            R["gradient"] = g.tolist()
        row["grad_rel_diff"] = float(np.max(np.abs(g - ref) / np.maximum(1.0, np.abs(ref))))
        print(json.dumps({k: v for k, v in row.items() if k != "units"}), flush=True)
        R["runs"].append(row)
    host = platform.node().split(".")[0]
    out = os.path.join(HERE, "results", "%s_isigc_%s.json" % (a.tag, host))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(R, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
