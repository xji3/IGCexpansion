"""L4: IS-IGC MG94 fit on BEAGLE (exact GL gradient) vs the legacy fit (jsonctmctree, central FD).

    .venv/bin/python benchmarks/l4_mg94_fit.py [--procs 4] [--engines legacy,beagle]

Case: ReCodonGeneconv MG94 on the YLR406C_YDL075W test alignment from the golden-case start
point; ReCodonGeneconv.get_mle (L-BFGS-B, legacy bounds) for each engine. Records the optimum,
x, function evaluations and wall time. Writes benchmarks/results/l4_mg94_fit_<host>.json.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--engines", default="legacy,beagle")
    a = ap.parse_args()
    R = {"case": "ReCodonGeneconv MG94 YLR406C_YDL075W test alignment, get_mle (L-BFGS-B)",
         "env": {"python": platform.python_version(), "platform": platform.platform(),
                 "cpu_count": os.cpu_count(), "date": time.strftime("%Y-%m-%d %H:%M"),
                 "procs": a.procs, "threads": a.threads}, "fits": {}}
    for eng in a.engines.split(","):
        with cases.quiet():
            o = cases.mg94_ylr406c(tempfile.mkdtemp())
        R["x0"] = list(map(float, o.x))
        if eng == "beagle":
            recodon.attach(o, processes=a.procs, threads=a.threads)
        t0 = time.perf_counter()
        with cases.quiet():
            r = o.get_mle(display=False)
        sec = time.perf_counter() - t0
        recodon.detach(o)
        R["fits"][eng] = {"fun": float(r.fun), "x": list(map(float, r.x)), "nfev": int(r.nfev), "nit": int(r.nit),
                          "message": str(r.message), "sec": sec, "grad_at_opt": list(map(float, r.jac))}
        print(eng, json.dumps(R["fits"][eng]), flush=True)
    if "legacy" in R["fits"] and "beagle" in R["fits"]:
        a0, b0 = R["fits"]["legacy"], R["fits"]["beagle"]
        R["dfun"] = b0["fun"] - a0["fun"]
        R["max_abs_dx"] = float(np.max(np.abs(np.array(b0["x"]) - a0["x"])))
        print("dfun", R["dfun"], "max|dx|", R["max_abs_dx"])
    out = os.path.join(HERE, "results", "l4_mg94_fit_%s.json" % platform.node().split(".")[0])
    with open(out, "w") as f:
        json.dump(R, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
