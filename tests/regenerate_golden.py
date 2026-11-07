"""Regenerate tests/golden/l2.json from the jsonctmctree engine.

    .venv/bin/python tests/regenerate_golden.py [case ...]

Run it only when the oracle itself is meant to change; the result is checked in.
"""

import json
import os
import platform
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
import scipy  # noqa: E402

import cases  # noqa: E402

OUT = os.path.join(HERE, "golden", "l2.json")


def main(names):
    old = json.load(open(OUT)) if os.path.exists(OUT) else {"values": {}}
    values = dict(old["values"])
    with tempfile.TemporaryDirectory() as d:
        for name in names:
            try:
                values[name] = cases.evaluate(name, d)
                print(name, values[name], file=sys.stderr)
            except FileNotFoundError as e:
                print("skip %s: missing %s" % (name, e), file=sys.stderr)
    import jsonctmctree

    doc = {
        "about": "jsonctmctree golden log-likelihoods at fixed points (plan L2); see tests/cases.py",
        "env": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                "jsonctmctree": getattr(jsonctmctree, "__version__", "0.2.0 (xji3 d9d26cc)"),
                "machine": platform.machine()},
        "values": values,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(doc, f, indent=1, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main(sys.argv[1:] or cases.CASES)
