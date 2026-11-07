"""M0 fair baselines for the BEAGLE decisions D1 and D3 (plan §9, §12).

    .venv/bin/python benchmarks/m0_baselines.py isigc [--procs 8]
    .venv/bin/python benchmarks/m0_baselines.py psigc [--procs 10] [--skip-status-quo]

isigc: MG94 on YGL062W_YBR218C (3721 states, 1129 codons), the IS-IGC benchmark pair.
    1 process: lnL and lnL + edge derivatives, as coded (sdnderi) and with the O(E)
    pre-order pass (sdngrad), uncompressed and compressed. Then the fair full-gradient
    baseline: compressed lnL + edge pass, plus 12 central differences for the 6 model
    parameters, over a process pool.
psigc: PSJSGeneconv HKY "One rate" with rate variation on YDR418W_YEL054C (1,461 scenes).
    Full objective: status quo (1 process, uncompressed), compressed 1 process, compressed
    over a process pool. Then the fair lnL + gradient baseline over the pool: compressed edge
    derivatives plus IGCexpansion's own per-parameter forward differences (delta 1e-8).

Each pool worker builds its own model once; pool start-up is excluded from timings.
BLAS is pinned to one thread per process. Results go to benchmarks/results/.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import platform  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), "tests")]

import numpy as np  # noqa: E402

import cases  # noqa: E402
import jct_compressed as jc  # noqa: E402

RESULTS = os.path.join(HERE, "results")
_W = {}  # per-process model


def _now():
    return time.perf_counter()


def _best(fn, repeats):
    """(min wall seconds, last result)."""
    best, out = float("inf"), None
    for _ in range(repeats):
        t0 = _now()
        out = fn()
        best = min(best, _now() - t0)
    return best, out


def _env(procs):
    return {"python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine(),
            "platform": platform.platform(), "cpu_count": os.cpu_count(), "procs": procs,
            "date": time.strftime("%Y-%m-%d %H:%M")}


# ---------------------------------------------------------------- IS-IGC (YGL062W)

def _isigc_init():
    d = tempfile.mkdtemp()
    with cases.quiet():
        _W["obj"] = cases.mg94_ygl062w(d)
    _W["x0"] = np.array(_W["obj"].x, float)


def _isigc_task(task):
    obj = _W["obj"]
    kind, x = task
    with cases.quiet():
        obj.update_by_x(np.asarray(x))
        scene = obj.get_scene()
    if kind == "edge":
        ll, ed, _ = jc.run(scene, "grad", True)
        return ll, ed.tolist()
    ll, _, _ = jc.run(scene, None, True)
    return ll, None


def isigc(procs, repeats):
    R = {"case": "MG94 YGL062W_YBR218C, 3721 states", "env": _env(procs)}
    _isigc_init()
    obj, x0 = _W["obj"], _W["x0"]
    with cases.quiet():
        t_scene, scene = _best(obj.get_scene, repeats)
    R["codons"] = len(scene["observed_data"]["iid_observations"])
    R["patterns"] = len(jc.compress_scene(scene)[1])
    R["sec_get_scene"] = t_scene
    one = {}
    for label, edge, comp in [("lnL", None, False), ("lnL_compressed", None, True),
                              ("lnL+sdnderi", "deri", False), ("lnL+sdngrad", "grad", False),
                              ("lnL+wdngrad_compressed", "grad", True)]:
        t, (ll, ed, rows) = _best(lambda: jc.run(scene, edge, comp), repeats)
        one[label] = {"sec": t, "lnL": ll, "rows": rows}
        print(label, one[label], file=sys.stderr, flush=True)
    R["one_process"] = one
    R["lnL_rel_diff_compressed"] = abs(one["lnL_compressed"]["lnL"] - one["lnL"]["lnL"]) / abs(one["lnL"]["lnL"])

    # Fair full gradient: 1 edge task + 12 central-difference lnL tasks over the pool.
    nproc_params = len(x0) - len(obj.edge_to_blen)
    h = 1e-5  # IGCexpansion's central-difference step for model parameters
    tasks = [("edge", x0.tolist())]
    for i in range(nproc_params):
        for s in (+1, -1):
            x = x0.copy()
            x[i] += s * h
            tasks.append(("lnl", x.tolist()))
    ctx = mp.get_context("spawn")
    with ctx.Pool(procs, initializer=_isigc_init) as pool:
        pool.map(_isigc_task, [("lnl", x0.tolist())] * procs)  # warm-up: every worker built
        t_pool, out = _best(lambda: pool.map(_isigc_task, tasks, chunksize=1), repeats)
    ll0, ed = out[0]
    g_model = [(out[1 + 2 * i][0] - out[2 + 2 * i][0]) / (2 * h) for i in range(nproc_params)]
    R["full_gradient_pool"] = {
        "sec": t_pool, "tasks": len(tasks), "procs": procs, "lnL": ll0,
        "dlnL_dx_model_central_fd": g_model, "dlnL_dlogt_edges": ed,
        "note": "compressed wdngrad edge pass + 2 x %d compressed lnL, chunksize 1" % nproc_params,
    }
    print("pool", t_pool, file=sys.stderr, flush=True)
    return R


# ---------------------------------------------------------------- PS-IGC (YDR418W)

def _psigc_init():
    d = tempfile.mkdtemp()
    with cases.quiet():
        obj = cases.psjs_ydr418w(d)
        # The constructor leaves every branch length at 1.0; apply self.x so that lnL and
        # gradient are evaluated at the same, well-defined point (see tests/cases.py).
        obj.unpack_x(np.array(obj.x, float))
    _W["obj"] = obj


def _scenes(obj):
    return [(n, cs) for n in sorted(obj.data.space_list) for cs in obj.data.space_codon_site_pair[n]]


def _psigc_lnl(group):
    obj = _W["obj"]
    with cases.quiet():
        return sum(jc.run(obj.get_scene(n, cs), None, True)[0] for n, cs in group)


def _psigc_grad(args):
    """Compressed version of PSJSGeneconv.loglikelihood_and_gradient_for_one_n over a group."""
    group, x = args
    obj = _W["obj"]
    x = np.asarray(x, float)
    delta = 1e-8
    with cases.quiet():
        obj.unpack_x(x)
        ll, ed = 0.0, 0.0
        for n, cs in group:
            a, e, _ = jc.run(obj.get_scene(n, cs), "deri", True)
            ll += a
            ed = ed + e
        m = len(x) - _n_rate_params(obj)
        other = np.zeros(m)
        for i in range(m):
            if obj.force is not None and i in obj.force:
                continue
            xd = x.copy()
            xd[i] += delta
            obj.unpack_x(xd)
            other[i] = (sum(jc.run(obj.get_scene(n, cs), None, True)[0] for n, cs in group) - ll) / delta
        obj.unpack_x(x)
    return ll, np.asarray(ed).tolist(), other.tolist()


def _n_rate_params(obj):
    k = len(obj.tree.edge_list)
    return k if obj.root_by_dup else k - 1


def _merge_edges(obj, ed):
    """IGCexpansion's mapping from edge derivatives to rate parameters (merged N0 edges)."""
    if obj.root_by_dup:
        return list(ed)
    el = obj.tree.edge_list
    merged = sum(ed[i] * (10.0 if "D1" in el[i] else 10.0 / 9.0) for i in range(len(ed)) if "N0" in el[i])
    return [merged] + [ed[i] for i in range(len(ed)) if "N0" not in el[i]]


def _groups(obj, scenes, ngroups):
    """Cost-balanced groups (greedy, cost = compressed rows), largest first for dynamic dispatch."""
    with cases.quiet():
        cost = [len(jc.compress_scene(obj.get_scene(n, cs))[1]) + 20 for n, cs in scenes]
    order = np.argsort(cost)[::-1]
    bins = [[] for _ in range(ngroups)]
    load = np.zeros(ngroups)
    for j in order:
        b = int(np.argmin(load))
        bins[b].append(scenes[j])
        load[b] += cost[j]
    return [b for b in bins if b], cost


def psigc(procs, repeats, skip_status_quo):
    R = {"case": "PSJSGeneconv HKY One rate, rate variation, YDR418W_YEL054C", "env": _env(procs)}
    _psigc_init()
    obj = _W["obj"]
    scenes = _scenes(obj)
    groups, cost = _groups(obj, scenes, 4 * procs)
    with cases.quiet():
        rows = [len(obj.get_scene(n, cs)["observed_data"]["iid_observations"]) for n, cs in scenes]
    R["scenes"] = len(scenes)
    R["rows_total"] = int(sum(rows))
    R["rows_compressed_total"] = int(sum(c - 20 for c in cost))
    R["x"] = list(map(float, obj.x))

    t0 = _now()
    with cases.quiet():
        for n, cs in scenes:
            obj.get_scene(n, cs)
    R["sec_get_scene_all_1proc"] = _now() - t0
    print("get_scene all", R["sec_get_scene_all_1proc"], file=sys.stderr, flush=True)

    if not skip_status_quo:
        t0 = _now()
        with cases.quiet():
            ll_sq = sum(obj._loglikelihood(n, cs)[0] for n, cs in scenes)
        R["status_quo_lnL_1proc"] = {"sec": _now() - t0, "lnL": float(ll_sq)}
        print("status quo", R["status_quo_lnL_1proc"], file=sys.stderr, flush=True)

    t0 = _now()
    ll_c = _psigc_lnl(scenes)
    R["compressed_lnL_1proc"] = {"sec": _now() - t0, "lnL": float(ll_c)}
    print("compressed 1p", R["compressed_lnL_1proc"], file=sys.stderr, flush=True)

    ctx = mp.get_context("spawn")
    with ctx.Pool(procs, initializer=_psigc_init) as pool:
        pool.map(_psigc_lnl, [[scenes[0]]] * procs)
        t, parts = _best(lambda: pool.map(_psigc_lnl, groups, chunksize=1), repeats)
        R["compressed_lnL_pool"] = {"sec": t, "lnL": float(sum(parts)), "groups": len(groups)}
        print("compressed pool lnL", R["compressed_lnL_pool"], file=sys.stderr, flush=True)

        x = list(map(float, obj.x))
        t, parts = _best(lambda: pool.map(_psigc_grad, [(g, x) for g in groups], chunksize=1), 1)
    ll = sum(p[0] for p in parts)
    ed = np.sum([p[1] for p in parts], axis=0)
    other = np.sum([p[2] for p in parts], axis=0)
    g = np.concatenate([other, _merge_edges(obj, ed)])
    R["compressed_lnL_and_gradient_pool"] = {
        "sec": t, "lnL": float(ll), "dlnL_dx": g.tolist(),
        "calls": len(scenes) * (1 + sum(1 for i in range(len(other)) if not (obj.force and i in obj.force))),
    }
    print("compressed pool grad", t, file=sys.stderr, flush=True)

    # Check the compressed gradient against IGCexpansion's own per-distance gradient at two distances.
    chk = {}
    with ctx.Pool(1, initializer=_psigc_init) as pool:
        for n in sorted(obj.data.space_list)[:1] + sorted(obj.data.space_list)[-1:]:
            grp = [(m, cs) for m, cs in scenes if m == n]
            with cases.quiet():
                f_sq, g_sq = obj.loglikelihood_and_gradient_for_one_n(n)
            llc, edc, oc = pool.apply(_psigc_grad, ((grp, x),))
            gc = np.concatenate([oc, _merge_edges(obj, edc)])
            chk[str(n)] = {"lnL_abs_diff": abs(-f_sq - llc),
                           "grad_max_abs_diff": float(np.max(np.abs(-np.asarray(g_sq) - gc))),
                           "grad_max_abs": float(np.max(np.abs(gc)))}
    R["gradient_check_vs_status_quo"] = chk
    print("check", chk, file=sys.stderr, flush=True)
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["isigc", "psigc"])
    ap.add_argument("--procs", type=int, default=None)
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--skip-status-quo", action="store_true")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.which == "isigc":
        R = isigc(a.procs or 8, a.repeats)
    else:
        R = psigc(a.procs or 10, a.repeats, a.skip_status_quo)
    tag = platform.node().split(".")[0]
    path = os.path.join(RESULTS, "m0_%s_%s.json" % (a.which, tag))
    with open(path, "w") as f:
        json.dump(R, f, indent=1)
        f.write("\n")
    print(path)


if __name__ == "__main__":
    main()
