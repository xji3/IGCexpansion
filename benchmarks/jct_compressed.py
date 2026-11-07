"""Pattern-compressed jsonctmctree calls: the fair single-engine baseline (plan §9, M0).

jsonctmctree evaluates every row of ``iid_observations``; its ``w`` observation reduction
only weights rows afterwards. Compression therefore means passing the unique rows and
their counts: ``wnnlogl`` / ``wdnderi`` / ``wdngrad`` with ``observation_reduction``.
"""

import numpy as np

import IGCexpansion  # noqa: F401  (installs the jsonctmctree shims; must precede jsonctmctree)
import jsonctmctree.interface


def compress_scene(scene):
    """Return (scene with unique observation rows, counts)."""
    od = scene["observed_data"]
    iid = np.asarray(od["iid_observations"])
    uniq, counts = np.unique(iid, axis=0, return_counts=True)
    out = dict(scene)
    out["observed_data"] = dict(od, iid_observations=uniq.tolist())
    return out, counts


def run(scene, edge=None, compressed=True):
    """lnL and optionally edge derivatives (edge in {None, 'deri', 'grad'}).

    'deri' is jsonctmctree's per-edge pass (sdnderi), 'grad' its O(E) pre-order pass (sdngrad).
    Returns (lnL, edge_derivs or None, number of rows evaluated).
    """
    if compressed:
        scene, counts = compress_scene(scene)
        red = {"observation_indices": list(range(len(counts))), "weights": counts.tolist()}
        reqs = [{"property": "wnnlogl", "observation_reduction": red}]
        if edge:
            reqs.append({"property": "wdn" + edge, "observation_reduction": red})
        rows = len(counts)
    else:
        reqs = [{"property": "snnlogl"}]
        if edge:
            reqs.append({"property": "sdn" + edge})
        rows = len(scene["observed_data"]["iid_observations"])
    j = jsonctmctree.interface.process_json_in({"scene": scene, "requests": reqs})
    assert j["status"] == "feasible", j["status"]
    r = j["responses"]
    return float(r[0]), (np.asarray(r[1], float) if edge else None), rows
