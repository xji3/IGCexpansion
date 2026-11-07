"""Drop-in ``process_json_in`` that evaluates jsonctmctree requests with BEAGLE (plan §3).

Handled requests: ``?nnlogl`` and ``?dnderi`` / ``?dngrad`` with observation code s, d or w
(``w`` uses ``observation_reduction``). Any other request, or a scene the compiler does not
accept, sends the *whole call* to jsonctmctree, so results never mix engines within a call.

Modes (``IGC_ENGINE`` or ``set_mode``):
  jsonctmctree  the original engine (default; install() is then a no-op)
  beagle        BEAGLE for handled calls, jsonctmctree for the rest
  shadow        both; jsonctmctree's answer is returned and every disagreement beyond
                tolerance is recorded in ``SHADOW_LOG`` and warned about

Other settings: ``IGC_BEAGLE_RESOURCE`` (cpu | gpu), ``IGC_BEAGLE_THREADS`` (CPU threads,
default 1), ``IGC_BEAGLE_MAX_GB`` (instance cache budget per process, default 4).
"""

import collections
import os
import time
import warnings

import numpy as np

import IGCexpansion  # noqa: F401  (jsonctmctree shims)
import jsonctmctree.interface as _jif

from ._lib import get_lib
from .scene import BeagleScene, UnsupportedScene, compile_scene, scene_values

ORIGINAL = _jif.process_json_in
MODES = ("jsonctmctree", "beagle", "shadow")
SHADOW_RTOL = 1e-9          # lnL (relative)
SHADOW_ATOL_DERI = 1e-7     # edge derivatives (absolute, scaled by max(1, |d|))
SHADOW_LOG = collections.deque(maxlen=1000)
SHADOW_MAX = {}             # request properties -> (calls, largest scaled difference), all calls
STATS = collections.Counter()

_state = {"mode": os.environ.get("IGC_ENGINE", "jsonctmctree"), "installed": False}
_cache = collections.OrderedDict()   # (signature, resource, threads, pid) -> BeagleScene


class Unhandled(Exception):
    pass


def set_mode(mode):
    if mode not in MODES:
        raise ValueError("IGC_ENGINE must be one of %s" % (MODES,))
    _state["mode"] = mode


def get_mode():
    return _state["mode"]


def _settings():
    res = os.environ.get("IGC_BEAGLE_RESOURCE", "cpu")
    threads = int(os.environ.get("IGC_BEAGLE_THREADS", "1"))
    max_bytes = float(os.environ.get("IGC_BEAGLE_MAX_GB", "4")) * 2 ** 30
    return res, threads, max_bytes


def clear_cache():
    while _cache:
        _cache.popitem(last=False)[1].inst.close()


def _parse(requests):
    out = []
    for r in requests:
        p = r.get("property", "")
        if len(p) != 7 or p[2] != "n":
            raise Unhandled(p)
        obs, edge, base = p[0], p[1], p[3:]
        if obs not in "sdw":
            raise Unhandled(p)
        if base == "logl" and edge == "n":
            kind = "logl"
        elif base in ("deri", "grad") and edge == "d":
            kind = "deri"
        else:
            raise Unhandled(p)
        if set(r) - {"property", "observation_reduction"}:
            raise Unhandled(p)
        red = None
        if obs == "w":
            red = r.get("observation_reduction")
            if red is None:
                raise Unhandled(p)
            red = (np.asarray(red["observation_indices"], dtype=np.int64), np.asarray(red["weights"], dtype=float))
        out.append((obs, kind, red))
    return out


def _get_scene(scene):
    res, threads, max_bytes = _settings()
    cs = compile_scene(scene)
    key = (cs.signature, res, threads, os.getpid())
    bs = _cache.pop(key, None)
    if bs is None:
        STATS["instances_created"] += 1
        bs = BeagleScene(get_lib(), cs, resource=res, threads=threads)
    _cache[key] = bs
    while len(_cache) > 1 and sum(b.cs.nbytes for b in _cache.values()) > max_bytes:
        _cache.popitem(last=False)[1].inst.close()
    return bs


def _beagle(j_in):
    reqs = _parse(j_in["requests"])
    scene = j_in["scene"]
    try:
        bs = _get_scene(scene)
        trips, distn, rates = scene_values(scene, bs.cs)
    except UnsupportedScene as e:
        raise Unhandled(str(e))
    bs.set_values(trips, distn, rates)
    total, site = bs.loglik()
    if not np.all(np.isfinite(site)):
        return {"status": "infeasible", "responses": None}
    need_deri = any(k == "deri" for _, k, _ in reqs)
    per_site = any(k == "deri" and o != "s" for o, k, _ in reqs)
    if need_deri:
        dtot, dper = bs.edge_derivatives(per_site=per_site)
    responses = []
    for obs, kind, red in reqs:
        if kind == "logl":
            if obs == "s":
                responses.append(float(np.sum(site)))
            elif obs == "d":
                responses.append(site.tolist())
            else:
                responses.append(float(np.dot(site[red[0]], red[1])))
        else:
            if obs == "s":
                responses.append(dtot.tolist())
            elif obs == "d":
                responses.append(dper.tolist())
            else:
                responses.append(np.tensordot(red[1], dper[red[0]], axes=1).tolist())
    return {"status": "feasible", "responses": responses}


def _compare(a, b):
    """Largest scaled difference between two jsonctmctree outputs (0 if identical status)."""
    if a["status"] != b["status"]:
        return float("inf")
    if a["responses"] is None:
        return 0.0
    worst = 0.0
    for x, y in zip(a["responses"], b["responses"]):
        x, y = np.asarray(x, float), np.asarray(y, float)
        if x.shape != y.shape:
            return float("inf")
        worst = max(worst, float(np.max(np.abs(x - y) / np.maximum(1.0, np.abs(y)))) if x.size else 0.0)
    return worst


def process_json_in(j_in, debug=False, seed=None):
    mode = _state["mode"]
    if mode == "jsonctmctree":
        return ORIGINAL(j_in, debug=debug, seed=seed)
    t0 = time.perf_counter()
    try:
        out = _beagle(j_in)
        STATS["beagle_calls"] += 1
    except Unhandled as e:
        STATS["fallback_calls"] += 1
        STATS["fallback:" + str(e)[:60]] += 1
        return ORIGINAL(j_in, debug=debug, seed=seed)
    t_beagle = time.perf_counter() - t0
    if mode == "beagle":
        return out
    t0 = time.perf_counter()
    ref = ORIGINAL(j_in, debug=debug, seed=seed)
    t_ref = time.perf_counter() - t0
    err = _compare(out, ref)
    tol = min(SHADOW_RTOL, SHADOW_ATOL_DERI)
    ok = err <= max(SHADOW_RTOL, SHADOW_ATOL_DERI) if np.isfinite(err) else False
    props = [r["property"] for r in j_in["requests"]]
    SHADOW_LOG.append({"properties": props, "max_scaled_diff": err, "ok": ok,
                       "sec_beagle": t_beagle, "sec_jsonctmctree": t_ref})
    k = ",".join(props)
    n, m = SHADOW_MAX.get(k, (0, 0.0))
    SHADOW_MAX[k] = (n + 1, max(m, err))
    STATS["sec_beagle"] += t_beagle
    STATS["sec_jsonctmctree"] += t_ref
    if not ok:
        STATS["shadow_mismatch"] += 1
        warnings.warn("BEAGLE/jsonctmctree mismatch %.3g (tol %.1g) on %s"
                      % (err, tol, [r["property"] for r in j_in["requests"]]), RuntimeWarning)
    return ref


def install(mode=None):
    """Route jsonctmctree.interface.process_json_in through this engine (idempotent)."""
    if mode is not None:
        set_mode(mode)
    if _state["mode"] == "jsonctmctree":
        uninstall()
        return False
    _jif.process_json_in = process_json_in
    _state["installed"] = True
    return True


def uninstall():
    _jif.process_json_in = ORIGINAL
    _state["installed"] = False
