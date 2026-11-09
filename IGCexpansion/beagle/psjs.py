"""PS-IGC (PSJSGeneconv) objective and gradient on resident worker processes (plan M2).

The engine is selectable:

  legacy        PSJSGeneconv's own methods, unchanged (default)
  jsonctmctree  resident workers, jsonctmctree with pattern compression
  beagle        resident workers, one BEAGLE action instance per (n, codon-site pair) scene
                ("BEAGLE-PS")

Select with ``IGC_PSJS_ENGINE`` (plus ``IGC_PSJS_PROCS``, ``IGC_PSJS_FD``) before importing
IGCexpansion, or per object with ``attach(obj, engine=..., processes=...)``. Both worker
engines share everything except the likelihood call, so they can be compared directly.

What the workers do differently from the legacy code:
* Q(n, cs) is built by ``QBuilder`` from precomputed structure: every rate is
  Q_mut[i, j] * r(site) + c_a * a(n) + c_b * b(n), with a = init/p (1 - (1-p)^n) and
  b = init/p (1-p)^n (PSJSModel.cal_IGC_transition_rate), instead of 3,360 Python calls per scene.
* Observations are compressed to unique rows per scene.
* Each worker owns a fixed, cost-balanced subset of the scenes; only x and the results cross
  process boundaries.
* The merged-N0 gradient uses unit weights (d lnL/dx0 = d_N0D1 + d_N0out, since both edges are
  x0 + constant in log space); the legacy code weights them 10 and 10/9 (plan §2.1).
* Model-parameter derivatives (``fd``): with BEAGLE the default is ``'gl'``, exact to quadrature
  tolerance (Gauss-Legendre Frechet integrals on BEAGLE's edge-derivative kernel, with the
  analytic complex-step differential matrices ``PSJSDifferentials``); ``'central'`` (h = 1e-5,
  the jsonctmctree default) and ``'forward'`` (the legacy 1e-8 step) are finite differences.
* BEAGLE instances are shared per worker and pattern-capacity bucket (ratio 2^(1/4)); each scene's
  tips are re-uploaded, padded with weight-0 all-ones patterns, before it is evaluated.
"""

import atexit
import itertools
import multiprocessing as mp
import os
import pickle
import time
import weakref

import numpy as np

ENGINES = ("legacy", "jsonctmctree", "beagle")
_POOLS = weakref.WeakKeyDictionary()     # PSJSGeneconv object -> PSJSPool


# ----------------------------------------------------------------------------- Q(n, cs)
def _other_pos(pos):
    return {0: (1, 2, 3), 1: (0, 3, 2), 2: (3, 0, 1), 3: (2, 1, 0)}[pos]


class QBuilder:
    """Vectorised values of PSJSModel's PM and IGC process definitions (HKY, 'One rate').

    The structure (row/column states, in the legacy enumeration order) depends only on the
    codon-site pair; the values are recomputed from the model's current parameters.
    """

    def __init__(self, psjsmodel, cs_list):
        if psjsmodel.pm_model != "HKY" or psjsmodel.IGC_pm != "One rate" or psjsmodel.n_js != 2:
            raise NotImplementedError("QBuilder covers HKY + 'One rate' with two paralogs")
        self.model = psjsmodel
        self.struct = {cs: (self._pm_struct(cs), self._igc_struct(cs)) for cs in cs_list}

    @staticmethod
    def _pm_struct(cs):
        rows, cols, q, site = [], [], [], []
        for sf in itertools.product(range(4), repeat=2):
            sf = sf + sf
            for i in range(2):
                for nt in range(4):
                    if sf[i] == nt:
                        continue
                    st = list(sf)
                    st[i] = st[i + 2] = nt
                    rows.append(sf)
                    cols.append(tuple(st))
                    q.append(4 * sf[i] + nt)
                    site.append(cs[i])
        return dict(rows=np.array(rows), cols=np.array(cols), q=np.array(q), site=np.array(site),
                    ca=np.zeros(len(q)), cb=np.zeros(len(q)))

    @staticmethod
    def _igc_struct(cs):
        ext = tuple(cs) * 2
        rows, cols, q, site, ca, cb = [], [], [], [], [], []
        for sf in itertools.product(range(4), repeat=4):
            for i in range(4):
                for nt in range(4):
                    if nt == sf[i]:
                        continue
                    st = list(sf)
                    st[i] = nt
                    so, po, poo = _other_pos(i)
                    a = 1.0 if nt == sf[po] else 0.0
                    b = 1.0 if (a and sf[so] == sf[poo]) else 0.0
                    rows.append(sf)
                    cols.append(tuple(st))
                    q.append(4 * sf[i] + nt)
                    site.append(ext[i])
                    ca.append(a)
                    cb.append(b)
            if sf[0] != sf[2] and sf[1] != sf[3]:
                for st in ((sf[2], sf[3], sf[2], sf[3]), (sf[0], sf[1], sf[0], sf[1])):
                    rows.append(sf)
                    cols.append(st)
                    q.append(0)
                    site.append(0)            # no point-mutation part
                    ca.append(0.0)
                    cb.append(1.0)
        return dict(rows=np.array(rows), cols=np.array(cols), q=np.array(q), site=np.array(site),
                    ca=np.array(ca), cb=np.array(cb))

    def coefficients(self):
        """(Q_mut flat, site multipliers [0, 1, r2, r3], init, p) from the model's current state."""
        m = self.model
        pm = m.PMModel
        mult = np.array([0.0, 1.0, pm.parameters.get("r2", np.nan), pm.parameters.get("r3", np.nan)])
        e = np.exp(np.asarray(m.x_IGC, dtype=float))
        if m.IGC_force is not None:
            for k, v in m.IGC_force.items():
                e[k] = v
        return np.asarray(pm.Q_mut, dtype=float).ravel(), mult, float(e[0]), float(e[1])

    @staticmethod
    def igc_ab(init, p, n):
        return init / p * (1 - (1 - p) ** n), init / p * (1 - p) ** n

    def values(self, cs, n, coef):
        qm, mult, init, p = coef
        a, b = self.igc_ab(init, p, n)
        pm, igc = self.struct[cs]
        out = []
        for s, (A, B) in ((pm, (0.0, 0.0)), (igc, (a, b))):
            out.append(qm[s["q"]] * mult[s["site"]] * (s["site"] > 0) + s["ca"] * A + s["cb"] * B)
        return out


def psjs_coefficients(x_js, psjsmodel):
    """(Q_mut flat[16], site multipliers [0, 1, r2, r3], init, p, root distn over (a,b,a,b)[16])
    as a pure function of x_js, mirroring PMModel (HKY) and PSJSModel ('One rate') including
    forced values. Works with complex x for complex-step differentiation."""
    m = psjsmodel
    pm_force = m.PMModel.force or {}
    igc_force = m.IGC_force or {}
    x = np.asarray(x_js)
    npm = len(m.x_pm)
    xp = [pm_force[i] if i in pm_force else np.exp(x[i]) for i in range(3)]
    pa, pc, pg, pt = xp[0] * xp[1], (1 - xp[0]) * xp[2], xp[0] * (1 - xp[1]), (1 - xp[0]) * (1 - xp[2])
    pi = [pa, pc, pg, pt]
    kappa = pm_force[3] if 3 in pm_force else np.exp(x[3])
    one = 1.0 + 0.0 * kappa
    Qb = [[0 * one, one, kappa, one], [one, 0 * one, one, kappa], [kappa, one, 0 * one, one], [one, kappa, one, 0 * one]]
    Qb = [[Qb[i][j] * pi[j] for j in range(4)] for i in range(4)]
    rate = sum(pi[i] * sum(Qb[i]) for i in range(4))
    qm = np.array([Qb[i][j] / rate for i in range(4) for j in range(4)])
    if m.rate_variation:
        r2 = pm_force[4] if 4 in pm_force else np.exp(x[4])
        r3 = pm_force[5] if 5 in pm_force else np.exp(x[5])
    else:
        r2 = r3 = np.nan
    mult = np.array([0 * one, one, r2 * one, r3 * one])
    xi = [igc_force[i] if i in igc_force else np.exp(x[npm + i]) for i in range(2)]
    pp = np.array([pi[a] * pi[b] for a in range(4) for b in range(4)])
    return qm, mult, xi[0], xi[1], pp / pp.sum()


def psjs_root_states():
    return [(a, b, a, b) for a in range(4) for b in range(4)]


class PSJSDifferentials:
    """Differential matrices dQ/dx and root-prior derivatives dpi/dx for PS-IGC scenes, by
    complex-step differentiation (exact to rounding) of QBuilder.values through psjs_coefficients."""

    def __init__(self, qb, shape, free, h=1e-30):
        self.qb, self.shape, self.free, self.h = qb, shape, list(free), h
        S = int(np.prod(shape))
        self.root_idx = np.ravel_multi_index(np.asarray(psjs_root_states()).T, shape)
        self.S = S

    def at(self, x_js):
        """Per-parameter coefficient derivatives at x_js (computed once per x)."""
        x_js = np.asarray(x_js, dtype=float)
        self.dcoef = []
        for i in self.free:
            xc = x_js.astype(complex)
            xc[i] += 1j * self.h
            self.dcoef.append(psjs_coefficients(xc, self.qb.model))

    def scene(self, cs, n):
        """D[k][g] triplets and dpi[k] for scene (n, cs) at the x given to ``at``."""
        pm, igc = self.qb.struct[cs]
        D, dpi = [], []
        for c in self.dcoef:
            vals = self.qb.values(cs, n, c[:4])
            Dk = []
            for st, v in zip((pm, igc), vals):
                dv = np.imag(v) / self.h
                Dk.append(_generator(st["rows"], st["cols"], dv, self.shape) if np.any(dv != 0) else None)
            D.append(Dk)
            d = np.zeros(self.S)
            d[self.root_idx] = np.imag(c[4]) / self.h
            dpi.append(d if np.any(d != 0) else None)
        return D, dpi


def _generator(rows, cols, v, shape):
    """Generator triplets from fixed structure and current values (zeros dropped, -rowsum diagonal)."""
    S = int(np.prod(shape))
    keep = v != 0.0
    r = np.ravel_multi_index(rows[keep].T, shape)
    c = np.ravel_multi_index(cols[keep].T, shape)
    vv = v[keep]
    rowsum = np.bincount(r, weights=vv, minlength=S)
    d = np.flatnonzero(rowsum)
    return (np.ascontiguousarray(np.concatenate([r, d]), dtype=np.int32),
            np.ascontiguousarray(np.concatenate([c, d]), dtype=np.int32),
            np.ascontiguousarray(np.concatenate([vv, -rowsum[d]]), dtype=np.float64))


# ----------------------------------------------------------------------------- scenes
def scene_list(obj):
    """All (n, cs) scenes of the objective; cs is None without rate variation."""
    rv = obj.psjsmodel.rate_variation
    out = []
    for n in sorted(obj.data.space_list):
        if rv:
            out += [(n, tuple(cs)) for cs in obj.data.space_codon_site_pair[n]]
        else:
            out.append((n, None))
    return out


def _observations(obj, n, cs):
    if obj.iid_observations is None:
        obj.cal_iid_observations()
    return obj.iid_observations[cs][n] if cs is not None else obj.iid_observations[n]


def partition(costs, k):
    """Greedy longest-processing-time partition of scene indices into k groups."""
    groups = [[] for _ in range(k)]
    load = np.zeros(k)
    for j in np.argsort(costs)[::-1]:
        g = int(np.argmin(load))
        groups[g].append(int(j))
        load[g] += costs[j]
    return [sorted(g) for g in groups], load


# ----------------------------------------------------------------------------- worker
class _Worker:
    """Holds one group of scenes; evaluates lnL, edge derivatives and FD model gradients."""

    def __init__(self, obj, scenes, engine, threads, resource, gl=None):
        import IGCexpansion  # noqa: F401  (shims)

        self.obj, self.engine = obj, engine
        self.scenes = scenes
        cs_all = sorted({cs if cs is not None else (1, 1) for _, cs in scenes})
        self.qb = QBuilder(obj.psjsmodel, cs_all)
        with _quiet():
            first = obj.get_scene(*scenes[0]) if obj.psjsmodel.rate_variation else obj.get_scene(scenes[0][0])
        self.shape = list(first["state_space_shape"])
        self.template = first
        self.items = []
        self.instances = {}
        if engine == "beagle":
            from ._lib import get_lib
            from .scene import BeagleScene, compile_scene, pattern_data
            lib = get_lib()
            self.topology = compile_scene(first)
        for n, cs in scenes:
            iid = np.asarray(_observations(obj, n, cs))
            uniq, counts = np.unique(iid, axis=0, return_counts=True)
            sc = dict(first)
            sc["observed_data"] = dict(first["observed_data"], iid_observations=uniq.tolist())
            item = {"n": n, "cs": cs if cs is not None else (1, 1), "scene": sc, "counts": counts}
            if engine == "beagle":
                # Shared instances per pattern-capacity bucket; each scene's tips are re-uploaded
                # (padded with weight-0 all-ones patterns) before it is evaluated.
                cap = _bucket(len(counts))
                item["data"] = pattern_data(self.topology, iid, capacity=cap)
                if cap not in self.instances:
                    self.instances[cap] = BeagleScene(lib, self.topology, capacity=cap, resource=resource,
                                                      threads=threads, gl=gl)
                item["bs"] = self.instances[cap]
            self.items.append(item)
        self.differentials = None
        if gl:
            self.differentials = PSJSDifferentials(self.qb, self.shape, gl["free"])

    def _set_x(self, x):
        with _quiet():
            self.obj.unpack_x(np.asarray(x, dtype=float))
        tj = self.obj.tree.tree_json
        self.rates = np.asarray(tj["edge_rate_scaling_factors"], dtype=float)
        self.coef = self.qb.coefficients()
        if self.differentials is not None:
            self.differentials.at(self.obj.psjsmodel.x_js)
        states, distn = self.obj.get_prior()
        self.distn = np.zeros(int(np.prod(self.shape)))
        np.add.at(self.distn, np.ravel_multi_index(np.asarray(states).T, self.shape), distn)
        self.root_states, self.root_probs = states, distn

    def _eval(self, item, edges):
        vals = self.qb.values(item["cs"], item["n"], self.coef)
        if self.engine == "beagle":
            bs = item["bs"]
            bs.load_patterns(item["data"])
            pm, igc = self.qb.struct[item["cs"]]
            trips = [_generator(s["rows"], s["cols"], v, self.shape) for s, v in zip((pm, igc), vals)]
            bs.set_values(trips, self.distn, self.rates)
            ll, site = bs.loglik()
            if not np.all(np.isfinite(site)):
                return -np.inf, None
            if not edges:
                return ll, None
            ed = bs.edge_derivatives()[0]
            if edges == "gl":
                D, dpi = self.differentials.scene(item["cs"], item["n"])
                return ll, (ed, bs.gl_gradient(D, dpi)[0])
            return ll, ed
        from .engine import ORIGINAL
        pm, igc = self.qb.struct[item["cs"]]
        sc = dict(item["scene"])
        sc["tree"] = dict(sc["tree"], edge_rate_scaling_factors=self.rates.tolist())
        sc["root_prior"] = {"states": self.root_states, "probabilities": self.root_probs}
        sc["process_definitions"] = [dict(row_states=s["rows"].tolist(), column_states=s["cols"].tolist(),
                                          transition_rates=v.tolist()) for s, v in zip((pm, igc), vals)]
        red = {"observation_indices": list(range(len(item["counts"]))), "weights": item["counts"].tolist()}
        reqs = [{"property": "wnnlogl", "observation_reduction": red}]
        if edges:
            reqs.append({"property": "wdnderi", "observation_reduction": red})
        out = ORIGINAL({"scene": sc, "requests": reqs})
        if out["status"] != "feasible":
            return -np.inf, None
        return float(out["responses"][0]), (np.asarray(out["responses"][1], float) if edges else None)

    def lnl(self, x):
        self._set_x(x)
        return sum(self._eval(it, False)[0] for it in self.items)

    def evaluate(self, x, grad, fd, h, free):
        t0 = time.perf_counter()
        x = np.asarray(x, dtype=float)
        self._set_x(x)
        ll, ed, gm = 0.0, 0.0, 0.0
        mode = "gl" if (grad and fd == "gl") else grad
        for it in self.items:
            a, e = self._eval(it, mode)
            ll += a
            if grad:
                if e is None:
                    ed = ed + np.nan
                elif fd == "gl":
                    ed, gm = ed + e[0], gm + e[1]
                else:
                    ed = ed + e
        out = {"ll": ll, "sec": 0.0}
        if grad:
            if fd == "gl":
                g = np.asarray(gm, dtype=float)
            else:
                g = np.zeros(len(free))
                for k, i in enumerate(free):
                    xp = x.copy()
                    xp[i] += h
                    if fd == "central":
                        xm = x.copy()
                        xm[i] -= h
                        g[k] = (self.lnl(xp) - self.lnl(xm)) / (2 * h)
                    else:
                        g[k] = (self.lnl(xp) - ll) / h
                self._set_x(x)
            out.update(edges=np.asarray(ed, dtype=float), model=g)
        out["sec"] = time.perf_counter() - t0
        out["rss_mb"] = _rss_mb()
        if self.instances:
            out["gl_clipped"] = sum(getattr(b, "gl_clipped", 0) for b in self.instances.values())
        return out


def _bucket(P, lo=8, ratio=2 ** 0.25):
    """Pattern-capacity bucket: the smallest lo * ratio^k >= P (rounded up). Ratio 2^(1/4) wastes
    ~9% of pattern work on padding on average (powers of two waste ~30%)."""
    if P <= lo:
        return lo
    k = int(np.ceil(np.log(P / lo) / np.log(ratio) - 1e-12))
    return int(np.ceil(lo * ratio ** k))


class _quiet:
    def __enter__(self):
        import contextlib
        import io
        self._c = contextlib.redirect_stdout(io.StringIO())
        return self._c.__enter__()

    def __exit__(self, *a):
        return self._c.__exit__(*a)


def _worker_main(conn, obj_path, scenes, engine, threads, resource, gl=None):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(v, "1")
    try:
        with open(obj_path, "rb") as f:
            obj = pickle.load(f)
        w = _Worker(obj, scenes, engine, threads, resource, gl)
        conn.send(("ready", _rss_mb()))
    except Exception as e:  # report construction failures to the parent
        import traceback
        conn.send(("error", traceback.format_exc()))
        return
    while True:
        msg = conn.recv()
        if msg[0] == "stop":
            break
        try:
            conn.send(("ok", w.evaluate(*msg[1:])))
        except Exception:
            import traceback
            conn.send(("error", traceback.format_exc()))


def _rss_mb():
    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return r / 2 ** 20 if os.uname().sysname == "Darwin" else r / 1024.0
    except Exception:
        return float("nan")


# ----------------------------------------------------------------------------- pool
class PSJSPool:
    def __init__(self, obj, engine="beagle", processes=None, threads=1, resource="cpu", fd=None,
                 h=None, merged_n0="fixed", distances=None, K=6, J_max=4):
        """``distances``: restrict the objective to these site distances n (tests, benchmarks)."""
        if engine not in ("jsonctmctree", "beagle"):
            raise ValueError("pool engine must be 'jsonctmctree' or 'beagle'")
        if fd is None:
            fd = "gl" if engine == "beagle" else "central"
        if fd not in ("gl", "central", "forward"):
            raise ValueError("fd (model-parameter gradient) must be 'gl', 'central' or 'forward'")
        if fd == "gl" and engine != "beagle":
            raise ValueError("exact (GL) parameter gradients need engine='beagle'")
        if merged_n0 not in ("fixed", "legacy"):
            raise ValueError("merged_n0 must be 'fixed' or 'legacy'")
        self.obj, self.engine, self.fd, self.merged_n0 = obj, engine, fd, merged_n0
        self.h = h if h is not None else (1e-5 if fd == "central" else 1e-8)
        k = processes or int(os.environ.get("IGC_PSJS_PROCS", "0")) or max(1, (os.cpu_count() or 2) // 2)
        self.scenes = scene_list(obj)
        if distances is not None:
            keep = set(distances)
            self.scenes = [sc for sc in self.scenes if sc[0] in keep]
        k = min(k, len(self.scenes))
        costs = []
        for n, cs in self.scenes:
            iid = np.asarray(_observations(obj, n, cs))
            costs.append(len(np.unique(iid, axis=0)) + 8.0)     # + fixed per-scene overhead
        self.groups, self.load = partition(np.asarray(costs), k)
        x_rate_len = len(obj.tree.edge_list) if obj.root_by_dup else len(obj.tree.edge_list) - 1
        m = len(obj.x) - x_rate_len
        self.m = m
        self.free = [i for i in range(m) if not (obj.force is not None and i in obj.force)]
        self.gl = dict(m=len(self.free), K=K, J_max=J_max, free=self.free) if fd == "gl" else None
        # The model goes to the workers through a file: a large pickle in the spawn arguments can
        # block Process.start() forever if a child dies while bootstrapping.
        import tempfile
        fd_, self._blob = tempfile.mkstemp(prefix="igc_psjs_", suffix=".pkl")
        with os.fdopen(fd_, "wb") as f:
            pickle.dump(obj, f)
        ctx = mp.get_context("spawn")
        self.conns, self.procs = [], []
        t0 = time.perf_counter()
        for g in self.groups:
            a, b = ctx.Pipe()
            p = ctx.Process(target=_worker_main, args=(b, self._blob, [self.scenes[j] for j in g], engine, threads,
                                                       resource, self.gl), daemon=True)
            p.start()
            self.conns.append(a)
            self.procs.append(p)
        self.worker_rss_mb = []
        for c, pr in zip(self.conns, self.procs):
            try:
                tag, val = self._recv(c, pr)
            except RuntimeError:
                self.close()
                raise
            if tag != "ready":
                self.close()
                raise RuntimeError("PS-IGC worker failed to start:\n" + val)
            self.worker_rss_mb.append(val)
        os.unlink(self._blob)
        self.sec_start = time.perf_counter() - t0
        self.last = {}
        atexit.register(self.close)

    def _run(self, x, grad):
        for c in self.conns:
            c.send(("eval", list(map(float, x)), grad, self.fd, self.h, self.free))
        res = []
        for c, pr in zip(self.conns, self.procs):
            tag, val = self._recv(c, pr)
            if tag != "ok":
                raise RuntimeError("PS-IGC worker failed:\n" + val)
            res.append(val)
        self.last = {"worker_sec": [r["sec"] for r in res], "worker_rss_mb": [r["rss_mb"] for r in res],
                     "gl_clipped": sum(r.get("gl_clipped", 0) for r in res)}
        return res

    @staticmethod
    def _recv(conn, proc):
        """Receive from a worker; never block forever on a worker that has died."""
        while not conn.poll(1.0):
            if not proc.is_alive():
                break
        try:
            return conn.recv()
        except EOFError:
            proc.join(timeout=5)
            raise RuntimeError("PS-IGC worker %d exited (code %s). Under 'spawn' the main script must be "
                               "importable (not read from stdin) and guarded by if __name__ == '__main__'."
                               % (proc.pid, proc.exitcode))

    def loglikelihood(self, x):
        return float(sum(r["ll"] for r in self._run(x, False)))

    def loglikelihood_and_gradient(self, x):
        """(lnL, d lnL / dx) at x, in PSJSGeneconv's x layout."""
        res = self._run(x, True)
        ll = float(sum(r["ll"] for r in res))
        ed = np.sum([r["edges"] for r in res], axis=0)
        gm = np.zeros(self.m)
        gm[self.free] = np.sum([r["model"] for r in res], axis=0)
        return ll, np.concatenate([gm, self._rate_grad(ed)])

    def _rate_grad(self, ed):
        obj = self.obj
        if obj.root_by_dup:
            return np.asarray(ed)
        el = obj.tree.edge_list
        w = (lambda e: 1.0) if self.merged_n0 == "fixed" else (lambda e: 10.0 if "D1" in e else 10.0 / 9.0)
        merged = sum(ed[i] * w(el[i]) for i in range(len(ed)) if "N0" in el[i])
        return np.concatenate([[merged], [ed[i] for i in range(len(ed)) if "N0" not in el[i]]])

    def close(self):
        blob = getattr(self, "_blob", None)
        if blob and os.path.exists(blob):
            os.unlink(blob)
        for c in getattr(self, "conns", []):
            try:
                c.send(("stop",))
            except Exception:
                pass
        for p in getattr(self, "procs", []):
            p.join(timeout=5)
            if p.is_alive():
                p.terminate()
        self.conns, self.procs = [], []


# ----------------------------------------------------------------------------- attach / install
def _settings():
    return dict(engine=os.environ.get("IGC_PSJS_ENGINE", "legacy"),
                processes=int(os.environ.get("IGC_PSJS_PROCS", "0")) or None,
                fd=os.environ.get("IGC_PSJS_FD") or None,
                threads=int(os.environ.get("IGC_BEAGLE_THREADS", "1")),
                resource=os.environ.get("IGC_BEAGLE_RESOURCE", "cpu"))


def get_pool(obj, **kw):
    pool = _POOLS.get(obj)
    if pool is None:
        s = _settings()
        s.update(kw)
        engine = s.pop("engine")
        if engine == "legacy":
            return None
        pool = PSJSPool(obj, engine=engine, **s)
        _POOLS[obj] = pool
    return pool


def attach(obj, engine="beagle", **kw):
    """Use a resident worker pool with ``engine`` for this object's objective (or detach with 'legacy')."""
    if engine not in ENGINES:
        raise ValueError("engine must be one of %s" % (ENGINES,))
    detach(obj)
    if engine != "legacy":
        _POOLS[obj] = PSJSPool(obj, engine=engine, **kw)
    return _POOLS.get(obj)


def detach(obj):
    pool = _POOLS.pop(obj, None)
    if pool is not None:
        pool.close()


_LEGACY = {}


def install():
    """Patch PSJSGeneconv so its objective methods use the pool selected for each object."""
    from IGCexpansion.PSJSGeneconv import PSJSGeneconv

    if _LEGACY:
        return
    _LEGACY["loglikelihood_and_gradient"] = PSJSGeneconv.loglikelihood_and_gradient
    _LEGACY["objective_wo_gradient"] = PSJSGeneconv.objective_wo_gradient

    def loglikelihood_and_gradient(self, display=False):
        pool = _POOLS.get(self) or get_pool(self)
        if pool is None:
            return _LEGACY["loglikelihood_and_gradient"](self, display)
        ll, g = pool.loglikelihood_and_gradient(self.x)
        self.ll = ll
        if display:
            print("Sum log likelihood = ", ll)
            print("All derivatives = ", g)
            print("Current exp x array = ", np.exp(self.x))
        return -ll, -g

    def objective_wo_gradient(self, display, x):
        pool = _POOLS.get(self) or get_pool(self)
        if pool is None:
            return _LEGACY["objective_wo_gradient"](self, display, x)
        self.unpack_x(x)
        ll = pool.loglikelihood(self.x)
        self.ll = ll
        if display:
            print("log likelihood = ", ll)
            print("Current x array = ", self.x)
        self.save_iteration(ll, x, [])
        return -ll

    PSJSGeneconv.loglikelihood_and_gradient = loglikelihood_and_gradient
    PSJSGeneconv.objective_wo_gradient = objective_wo_gradient


def uninstall():
    from IGCexpansion.PSJSGeneconv import PSJSGeneconv

    for k, f in _LEGACY.items():
        setattr(PSJSGeneconv, k, f)
    _LEGACY.clear()


# ----------------------------------------------------------------------------- Godambe (2-D)
def _x_from_2d(obj, x2d):
    """PSJSGeneconv's 2-D reparameterisation x = [ln tau, logit p] (``_sitewise_lnL_2d_x``)."""
    a, b = float(x2d[0]), float(x2d[1])
    igc = [a + b - np.log1p(np.exp(b)), b - np.log1p(np.exp(b))]
    return np.concatenate([obj.psjsmodel.x_pm, igc, obj.x[len(obj.psjsmodel.x_js):]])


def pair_scores_2d(obj, x2d, worker=None):
    """{site pair: d lnL_pair / d (ln tau, logit p)} from exact GL per-site scores."""
    npm = len(obj.psjsmodel.x_pm)
    free = [npm, npm + 1]                      # log init, log p
    if worker is None:
        worker = _Worker(obj, scene_list(obj), "beagle", 1, "cpu", dict(m=2, K=6, free=free))
    x = _x_from_2d(obj, x2d)
    worker._set_x(x)
    sig = 1.0 / (1.0 + np.exp(-float(x2d[1])))
    out = {}
    for it in worker.items:
        bs = it["bs"]
        bs.load_patterns(it["data"])
        vals = worker.qb.values(it["cs"], it["n"], worker.coef)
        pm, igc = worker.qb.struct[it["cs"]]
        bs.set_values([_generator(s["rows"], s["cols"], v, worker.shape) for s, v in zip((pm, igc), vals)],
                      worker.distn, worker.rates)
        bs.loglik()
        D, dpi = worker.differentials.scene(it["cs"], it["n"])
        _, per = bs.gl_gradient(D, dpi, per_site=True)            # [nobs, 2] in (log init, log p)
        g2 = np.stack([per[:, 0], (per[:, 0] + per[:, 1]) * (1.0 - sig)], axis=1)
        n, cs = it["n"], it["cs"]
        pairs = obj.data.space_idx_pairs[cs][n] if obj.psjsmodel.rate_variation else obj.data.space_idx_pairs[n]
        for p, gv in zip(pairs, g2):
            out[p] = gv
    return out, worker


def godambe_2d(obj, x2d, step=1e-5):
    """Per-pair gradients and Hessians in (ln tau, logit p) and the Godambe matrix H J^-1 H,
    as PSJSGeneconv.get_Godambe_matrix but with exact (GL) scores; Hessians are central
    differences of the exact scores (relative step ``step``, as the legacy code)."""
    x2d = np.asarray(x2d, dtype=float)
    g, w = pair_scores_2d(obj, x2d)
    hx = np.abs(x2d * step)
    cols = []
    for i in range(2):
        xp, xm = x2d.copy(), x2d.copy()
        xp[i] += hx[i]
        xm[i] -= hx[i]
        gp, _ = pair_scores_2d(obj, xp, w)
        gm, _ = pair_scores_2d(obj, xm, w)
        cols.append({p: (gp[p] - gm[p]) / (2 * hx[i]) for p in g})
    hess = {p: np.stack([cols[0][p], cols[1][p]], axis=1) for p in g}     # [i, j] = d g_i / d x_j
    hess = {p: 0.5 * (h + h.T) for p, h in hess.items()}
    G = np.array([g[p] for p in g])
    Hs = np.array([hess[p] for p in g])
    H = -Hs.sum(axis=0) / len(Hs)
    J = sum(np.outer(u, u) for u in G) / len(G)
    return g, hess, H @ np.linalg.solve(J, H)
