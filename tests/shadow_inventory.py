"""Run the L2b legacy-fixture inventory with IGC_ENGINE=shadow and report BEAGLE vs jsonctmctree.

    .venv/bin/python tests/shadow_inventory.py

Every jsonctmctree call made by the ~50 legacy fixtures (ReCodon, IndCodon, JS, PSJS, HMM...)
is answered by both engines; jsonctmctree's answer is used, and the largest scaled difference
per call is recorded. Writes tests/golden/shadow_inventory.json (the inventory's own outputs
are written to a temporary directory, not over the checked-in report).
"""

import json
import os
import sys
import tempfile
import time

os.environ["IGC_ENGINE"] = "shadow"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

import IGCexpansion  # noqa: E402,F401
from IGCexpansion.beagle import engine  # noqa: E402
import inventory_l2b as inv  # noqa: E402


def main():
    assert engine.get_mode() == "shadow"
    tmp = tempfile.mkdtemp()
    inv.OUT_JSON = os.path.join(tmp, "l2b.json")
    inv.OUT_MD = os.path.join(tmp, "l2b.md")
    sys.argv = [sys.argv[0]]
    t0 = time.perf_counter()
    inv.main()
    wall = time.perf_counter() - t0
    log = list(engine.SHADOW_LOG)
    out = {
        "wall_sec": wall,
        "stats": dict(engine.STATS),
        "shadow_calls_logged": len(log),
        "max_scaled_diff": max((m for _, m in engine.SHADOW_MAX.values()), default=None),
        "mismatches": [r for r in log if not r["ok"]][:20],
        "by_request": {k: {"calls": n, "max_scaled_diff": m} for k, (n, m) in engine.SHADOW_MAX.items()},
        "note": "times are summed over all shadow calls (engine.STATS)",
    }
    path = os.path.join(HERE, "golden", "shadow_inventory.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, default=str)
        f.write("\n")
    print(json.dumps({k: out[k] for k in ("stats", "max_scaled_diff", "by_request")}, indent=1, default=str))


if __name__ == "__main__":
    main()
