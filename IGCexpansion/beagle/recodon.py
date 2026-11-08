"""IS-IGC (ReCodonGeneconv, MG94 or HKY) objective with exact gradients on BEAGLE (plan M4).

``loglikelihood_and_gradient(obj)`` replaces ReCodonGeneconv's method of the same name:

* lnL and the edge derivatives come from one BEAGLE action instance resident on the object
  (pattern-compressed tips uploaded once).
* The process-parameter derivatives are exact to quadrature tolerance: Gauss-Legendre Frechet
  integrals on BEAGLE's edge-derivative kernel (``BeagleScene.gl_gradient``) with the analytic
  Jacobian ``IsIgcJacobian``. The legacy code uses central differences with relative step 1e-6
  (2m extra likelihoods); jsonctmctree's ``sdngrad`` is not used.
* f = -lnL and g = -[d/dx_process, d/dlog t] keep the legacy layout. Forced process
  parameters get 0, as in the legacy code.

BEAGLE is the default engine (D3). ``IGC_RECODON_ENGINE=legacy`` before importing
IGCexpansion restores the legacy code (jsonctmctree + central differences); the variable is read
at every call, so it can also be changed at run time. If the BEAGLE library does not load, the
legacy code is used with one warning. ``attach(obj, processes=...)`` selects BEAGLE for one
object regardless. ``IGC_RECODON_PROCS`` > 1 evaluates a single pair on a ``RecodonPool`` that
splits its GL edges over worker processes (use threads first: ``IGC_BEAGLE_THREADS``);
``IGC_BEAGLE_RESOURCE=gpu`` uses the GPU. Homo_Omega models are not covered by the analytic
Jacobian; they use ``scene_fd_jacobian``.

JointAnalysis (both objective methods) runs on a ``RecodonPool`` (by default one process per
pair) instead of the legacy per-call ``multiprocessing.Pool``: workers build and keep their
pair objects and BEAGLE instances, and only x crosses process boundaries. Like any spawn-based
pool, the driver must be a file with an ``if __name__ == "__main__":`` guard.
"""

import itertools
import os
import weakref

import numpy as np

from ._lib import get_lib
from .jacobian import scene_fd_jacobian
from .scene import BeagleScene, compile_scene, scene_values

ENGINES = ("legacy", "beagle")
_STATE = weakref.WeakKeyDictionary()       # ReCodonGeneconv -> _Resident


# ----------------------------------------------------------------------------- generator values
def _flat_generator(r, c, v, S):
    """Generator triplets from flat (r, c) structure and values (zeros dropped, -rowsum diagonal)."""
    keep = v != 0.0
    r, c, vv = r[keep], c[keep], v[keep]
    rowsum = np.bincount(r, weights=vv, minlength=S)
    d = np.flatnonzero(rowsum)
    return (np.ascontiguousarray(np.concatenate([r, d]), dtype=np.int32),
            np.ascontiguousarray(np.concatenate([c, d]), dtype=np.int32),
            np.ascontiguousarray(np.concatenate([vv, -rowsum[d]]), dtype=np.float64))


class IsIgcBuilder:
    """Vectorised ReCodonGeneconv process values and root prior as functions of x_process
    (log scale), complex-safe so that complex-step differentiation is exact to rounding.

    Rates (n = 61 codons or 4 nucleotides, Qb = the normalised single-copy matrix):
      geneconv  (a,b)->(a,c), c != b:  Qb[b,c] + [a != b][c == a] T(a,b)
                (a,b)->(c,b), c != a:  Qb[a,c] + [a != b][c == b] T(a,b)
      basic     (a,a)->(c,c), c != a:  Qb[a,c]
    T(a,b) = tau if a, b are synonymous (always for HKY), else the IGC nonsynonymous
    contribution (omega*tau, IGC_Omega*tau or Tau_Omega). This is get_MG94Geneconv_and_MG94 and
    get_HKYGeneconv; ``check`` compares against the object's own scene.
    """

    def __init__(self, obj):
        if obj.Model not in ("MG94", "HKY"):
            raise NotImplementedError(obj.Model)
        if obj.Model == "MG94" and obj.use_Homo_Omega():
            raise NotImplementedError("Homo_Omega")
        self.model = obj.Model
        self.force = {i: v for i, v in (obj.Force or {}).items() if i < len(obj.x) - len(obj.edge_to_blen)}
        self.igc_omega = obj.Model == "MG94" and obj.use_IGC_Omega()
        self.tau_omega = self.igc_omega and obj.IGC_Omega is None
        if self.model == "MG94":
            units = obj.codon_nonstop
            nts = np.array([["ACGT".index(ch) for ch in u] for u in units])          # [61, 3]
            aa = np.array([obj.codon_table[u] for u in units])
        else:
            units = list("ACGT")
            nts = np.arange(4)[:, None]
            aa = np.array(units)
        n = self.n = len(units)
        self.S = n * n
        self.nts = nts
        # single-copy structure: one-position differences
        diff = nts[:, None, :] != nts[None, :, :]
        one = diff.sum(axis=2) == 1
        pos = np.argmax(diff, axis=2)
        ia, ib = np.nonzero(one)
        self.qi, self.qj = ia, ib
        p = pos[ia, ib]
        self.to_nt = nts[ib, p]
        fr = nts[ia, p]
        self.ts = ((fr + self.to_nt) % 2 == 0) & (fr != self.to_nt)      # A<->G (0,2), C<->T (1,3)
        self.ns = aa[ia] != aa[ib]
        # paired-state structure (geneconv process)
        rows, cols, qidx, tflag = [], [], [], []
        nonsyn = aa[:, None] != aa[None, :]
        A, B, C = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
        A, B, C = A.ravel(), B.ravel(), C.ravel()
        for first in (False, True):          # (a,b)->(a,c) with c != b ; (a,b)->(c,b) with c != a
            if not first:
                m = C != B
                r, c_ = A * n + B, A * n + C
                q = B * n + C
                t = (A != B) & (C == A)
            else:
                m = C != A
                r, c_ = A * n + B, C * n + B
                q = A * n + C
                t = (A != B) & (C == B)
            rows.append(r[m]); cols.append(c_[m]); qidx.append(q[m]); tflag.append(t[m])
        self.g_r = np.concatenate(rows)
        self.g_c = np.concatenate(cols)
        self.g_q = np.concatenate(qidx)
        tf = np.concatenate(tflag)
        ab_nonsyn = nonsyn[self.g_r // n, self.g_r % n]
        self.g_tau = tf & ~ab_nonsyn
        self.g_tns = tf & ab_nonsyn
        a = np.repeat(np.arange(n), n)
        c = np.tile(np.arange(n), n)
        m = a != c
        self.b_r, self.b_c, self.b_q = (a * n + a)[m], (c * n + c)[m], (a * n + c)[m]
        self.root = np.arange(n) * n + np.arange(n)

    def params(self, xp):
        v = np.exp(np.asarray(xp))
        for i, f in self.force.items():
            v[i] = f
        pi = np.stack([v[0] * v[1], (1 - v[0]) * v[2], v[0] * (1 - v[1]), (1 - v[0]) * (1 - v[2])])
        kappa = v[3]
        if self.model == "HKY":
            return pi, kappa, 1.0, v[4], v[4]
        omega = v[4]
        if self.igc_omega:
            tau = v[6]
            tns = v[5] if self.tau_omega else v[5] * tau
        else:
            tau = v[5]
            tns = omega * tau
        return pi, kappa, omega, tau, tns

    def values(self, xp):
        """(basic rates, geneconv rates, root distribution over the n*n paired states)."""
        pi, kappa, omega, tau, tns = self.params(xp)
        prior = np.prod(pi[self.nts], axis=1)
        prior = prior / prior.sum()
        q = pi[self.to_nt] * np.where(self.ts, kappa, 1.0) * np.where(self.ns, omega, 1.0)
        Q = np.zeros(self.n * self.n, dtype=q.dtype)
        Q[self.qi * self.n + self.qj] = q
        rowsum = Q.reshape(self.n, self.n).sum(axis=1)
        Q = Q / (prior @ rowsum)
        geneconv = Q[self.g_q] + np.where(self.g_tau, tau, 0.0) + np.where(self.g_tns, tns, 0.0)
        basic = Q[self.b_q]
        distn = np.zeros(self.S, dtype=prior.dtype)
        distn[self.root] = prior
        return basic, geneconv, distn

    def generators(self, xp):
        basic, geneconv, distn = self.values(xp)
        return [_flat_generator(self.b_r, self.b_c, basic, self.S),
                _flat_generator(self.g_r, self.g_c, geneconv, self.S)], distn

    def check(self, obj, cs, tol=1e-13):
        """Largest relative difference to the object's own scene (raises above ``tol``)."""
        import scipy.sparse as sp

        trips, distn, _ = scene_values(obj.get_scene(), cs)
        mine, d2 = self.generators(np.asarray(obj.x_process, float))
        err = float(np.max(np.abs(distn - d2)))
        for (r, c, v), (r2, c2, v2) in zip(trips, mine):
            M = sp.csr_matrix((v, (r, c)), shape=(self.S, self.S))
            M2 = sp.csr_matrix((v2, (r2, c2)), shape=(self.S, self.S))
            scale = max(1e-300, abs(M).max())
            err = max(err, abs(M - M2).max() / scale)
        if err > tol:
            raise AssertionError("IsIgcBuilder differs from the scene by %.3g" % err)
        return err


class IsIgcJacobian:
    """dQ_g/dx_k (generator triplets) and dpi/dx_k by complex-step differentiation of
    ``IsIgcBuilder.values`` (exact to rounding)."""

    def __init__(self, builder, free, h=1e-30):
        self.b, self.free, self.h = builder, list(free), h

    def at(self, xp):
        b = self.b
        D, dpi = [], []
        for i in self.free:
            xc = np.asarray(xp, dtype=float).astype(complex)
            xc[i] += 1j * self.h
            basic, geneconv, distn = b.values(xc)
            Dk = []
            for r, c, v in ((b.b_r, b.b_c, basic), (b.g_r, b.g_c, geneconv)):
                dv = np.imag(v) / self.h
                Dk.append(_flat_generator(r, c, dv, b.S) if np.any(dv != 0) else None)
            D.append(Dk)
            d = np.imag(distn) / self.h
            dpi.append(d if np.any(d != 0) else None)
        return D, dpi


# ----------------------------------------------------------------------------- resident instance
class _Resident:
    """BEAGLE state kept on one ReCodonGeneconv object (per process)."""

    def __init__(self, obj, resource="cpu", threads=1, gl=None):
        self.pid = os.getpid()
        scene = obj.get_scene()
        self.cs = compile_scene(scene)
        self.m = len(obj.x) - len(obj.edge_to_blen)
        self.free = [i for i in range(self.m) if not (obj.Force and i in obj.Force)]
        g = dict(m=max(1, len(self.free)))
        g.update(gl or {})
        self.bs = BeagleScene(get_lib(), self.cs, resource=resource, threads=threads, gl=g)
        self.obs_key = _obs_key(scene)
        try:
            self.builder = IsIgcBuilder(obj)
            self.builder_err = self.builder.check(obj, self.cs)
            self.jac = IsIgcJacobian(self.builder, self.free)
        except NotImplementedError:
            self.builder = self.jac = None


def _obs_key(scene):
    od = scene["observed_data"]
    t = scene["tree"]
    return (tuple(od["nodes"]), tuple(od["variables"]), len(od["iid_observations"]),
            tuple(t["row_nodes"]), tuple(t["column_nodes"]), tuple(t["edge_processes"]))


def _resident(obj, **kw):
    st = _STATE.get(obj)
    if st is None or st.pid != os.getpid():
        s = _settings()
        s.pop("engine")
        s.update(kw)
        st = _Resident(obj, **s)
        _STATE[obj] = st
    return st


def evaluate(obj, grad=True, per_site=False, edges=None, primary=True, **kw):
    """(lnL, d lnL/dx in the legacy layout or None, per-site scores [nsites, len(x)] or None)
    at obj.x. The model state of ``obj`` is left at obj.x.

    ``edges``/``primary`` split one gradient over several instances (``RecodonPool``): GL runs
    on ``edges`` only (all if None); the edge derivatives and the root-prior term are added only
    when ``primary``. Without ``primary`` the edge part of g is 0.
    """
    obj.update_by_x()
    st = _resident(obj, **kw)
    scene = obj.get_scene()
    if _obs_key(scene) != st.obs_key:
        raise RuntimeError("the scene's tree or observations changed; detach(obj) first")
    bs = st.bs
    trips, distn, rates = scene_values(scene, st.cs)
    bs.set_values(trips, distn, rates)
    ll, site = bs.loglik()
    if not np.all(np.isfinite(site)):
        return -np.inf, (np.full(len(obj.x), np.nan) if grad else None), None
    if not grad:
        return float(ll), None, None
    E = len(obj.edge_to_blen)
    if primary:
        ed, ed_site = bs.edge_derivatives(per_site=per_site)
    else:
        ed, ed_site = np.zeros(E), (np.zeros((len(site), E)) if per_site else None)
    other = np.zeros(st.m)
    other_site = np.zeros((len(site), st.m)) if per_site else None
    if st.free and (edges is None or len(edges) or primary):
        if st.jac is not None:
            D, dpi = st.jac.at(np.asarray(obj.x_process, float))
        else:
            x0 = np.array(obj.x, float)

            def scene_at(xx):
                obj.update_by_x(xx)
                return obj.get_scene()
            D, dpi = scene_fd_jacobian(scene_at, x0, st.free, st.cs)
            obj.update_by_x(x0)
        g, per = bs.gl_gradient(D, dpi, per_site=per_site, edges=edges)
        if edges is not None and primary:
            g2, per2 = bs.gl_gradient([[None] * st.cs.G for _ in D], dpi, per_site=per_site)
            g, per = g + g2, (per + per2 if per_site else None)
        other[st.free] = g
        if per_site:
            other_site[:, st.free] = per
    g = np.concatenate([other, ed])
    per = np.concatenate([other_site, ed_site], axis=1) if per_site else None
    return float(ll), g, per


def edge_costs(obj, gl=None):
    """Relative GL cost of each edge of the object's scene (for splitting it over processes):
    quadrature nodes plus the Taylor work t ||Q - mu I||, at the object's current x."""
    from .scene import BeagleScene, gl_nodes, generator_triplets

    scene = obj.get_scene()
    shape = [int(v) for v in scene["state_space_shape"]]
    S = int(np.prod(shape))
    trips = [generator_triplets(shape, p) for p in scene["process_definitions"]]
    norms = [BeagleScene._norm_inf_shifted(t, S) for t in trips]
    # one action costs ~ (Taylor steps) x nnz(Q) x patterns
    P = len(np.unique(np.asarray(scene["observed_data"]["iid_observations"]), axis=0))
    tr = scene["tree"]
    out = []
    for t, g in zip(tr["edge_rate_scaling_factors"], tr["edge_processes"]):
        g = int(g)
        J, K = gl_nodes(float(t), norms[g], gl)
        out.append(2.0 * (J * K + float(t) * norms[g]) * (len(trips[g][2]) + S) * P * 1e-6)
    return np.asarray(out)



# ----------------------------------------------------------------------------- resident pool
def plan_units(costs, processes):
    """Split pairs into work units (pair, edges, primary), one list of units per process.

    ``costs[i]``: per-edge GL costs of pair i. With at least as many pairs as processes, whole
    pairs are balanced over processes (LPT). Otherwise each pair gets processes in proportion to
    its GL cost and its edges are LPT-partitioned over them; every unit repeats the pair's
    post- and pre-order passes, and the first one also adds lnL, edge derivatives and dpi.
    """
    from .psjs import partition

    n = len(costs)
    tot = np.array([float(np.sum(c)) for c in costs])
    if n >= processes:
        groups, _ = partition(tot + 1.0, processes)
        return [[(i, None, True) for i in g] for g in groups if len(g)]
    share = np.maximum(1, np.floor(processes * tot / tot.sum()).astype(int))
    while share.sum() < processes:
        share[np.argmax(tot / share)] += 1
    while share.sum() > processes:
        j = np.argmax(np.where(share > 1, share / np.maximum(tot, 1e-300), -np.inf))
        share[j] -= 1
    out = []
    for i in range(n):
        k = int(min(share[i], max(1, int(np.sum(np.asarray(costs[i]) > 0)))))
        if k == 1:
            out.append([(i, None, True)])
            continue
        groups, load = partition(np.asarray(costs[i], float) + 1e-9, k)
        order = np.argsort(load)            # the lightest group also carries lnL, edges, dpi
        groups = [groups[j] for j in order if len(groups[j])]
        out += [[(i, g, j == 0)] for j, g in enumerate(groups)]
    return out


def _pool_worker(conn, blob, units, threads, resource):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(v, "1")
    import pickle
    import time
    from .psjs import _quiet, _rss_mb
    try:
        with open(blob, "rb") as f:
            objs = pickle.load(f)
        need = sorted({u[0] for u in units})
        with _quiet():
            for i in need:
                objs[i].update_by_x()
                _resident(objs[i], threads=threads, resource=resource)
        conn.send(("ready", _rss_mb()))
    except Exception:
        import traceback
        conn.send(("error", traceback.format_exc()))
        return
    while True:
        msg = conn.recv()
        if msg[0] == "stop":
            break
        try:
            _, xs, grad = msg
            t0 = time.perf_counter()
            res = []
            with _quiet():
                for i, edges, primary in units:
                    objs[i].update_by_x(np.asarray(xs[i], float))
                    if not grad and not primary:
                        res.append((i, 0.0, None))
                        continue
                    ll, g, _ = evaluate(objs[i], grad=grad, edges=edges, primary=primary)
                    res.append((i, ll if primary else 0.0, g))
            conn.send(("ok", dict(res=res, sec=time.perf_counter() - t0, rss_mb=_rss_mb())))
        except Exception:
            import traceback
            conn.send(("error", traceback.format_exc()))


class RecodonPool:
    """Resident BEAGLE workers for one or several ReCodonGeneconv objects (IS-IGC pairs).

    Each worker builds its pair objects once (from a pickle file) and keeps their BEAGLE
    instances; a call sends only the x of each pair. ``processes`` workers with ``threads``
    BEAGLE threads each. With fewer pairs than processes, a pair's GL edges are split over
    several workers (``plan_units``), which lowers the latency of a single-pair gradient.
    """

    def __init__(self, objs, processes=None, threads=1, resource="cpu"):
        import multiprocessing as mp
        import pickle
        import tempfile
        import time

        self.objs = list(objs)
        k = processes or max(1, (os.cpu_count() or 2) // 2)
        with _quiet_ctx():
            costs = []
            for o in self.objs:
                o.update_by_x()
                costs.append(edge_costs(o))
        self.units = plan_units(costs, k)
        fd_, self._blob = tempfile.mkstemp(prefix="igc_recodon_", suffix=".pkl")
        with os.fdopen(fd_, "wb") as f:
            pickle.dump(self.objs, f)
        ctx = mp.get_context("spawn")
        self.conns, self.procs = [], []
        t0 = time.perf_counter()
        try:
            for u in self.units:
                a, b = ctx.Pipe()
                p = ctx.Process(target=_pool_worker, args=(b, self._blob, u, threads, resource), daemon=True)
                p.start()
                self.conns.append(a)
                self.procs.append(p)
            self.worker_rss_mb = []
            for c, pr in zip(self.conns, self.procs):
                tag, val = self._recv(c, pr)
                if tag != "ready":
                    raise RuntimeError("IS-IGC worker failed to start:\n" + val)
                self.worker_rss_mb.append(val)
        except BaseException:
            self.close()
            raise
        os.unlink(self._blob)
        self.sec_start = time.perf_counter() - t0
        self.last = {}
        import atexit
        atexit.register(self.close)

    @staticmethod
    def _recv(conn, proc):
        from .psjs import PSJSPool
        return PSJSPool._recv(conn, proc)

    def evaluate(self, xs, grad=True):
        """([lnL_i], [d lnL_i / dx_i] or None) for the pairs at their parameter vectors ``xs``."""
        xs = [list(map(float, x)) for x in xs]
        for c in self.conns:
            c.send(("eval", xs, grad))
        n = len(self.objs)
        ll = np.zeros(n)
        gs = [np.zeros(len(x)) for x in xs] if grad else None
        secs = []
        for c, pr in zip(self.conns, self.procs):
            tag, val = self._recv(c, pr)
            if tag != "ok":
                raise RuntimeError("IS-IGC worker failed:\n" + val)
            secs.append(val["sec"])
            for i, l, g in val["res"]:
                ll[i] += l
                if grad:
                    gs[i] += g
        self.last = {"worker_sec": secs}
        return ll, gs

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


def _quiet_ctx():
    from .psjs import _quiet
    return _quiet()


# ----------------------------------------------------------------------------- selection
def _settings():
    return dict(engine=os.environ.get("IGC_RECODON_ENGINE", "beagle"),
                threads=int(os.environ.get("IGC_BEAGLE_THREADS", "1")),
                resource=os.environ.get("IGC_BEAGLE_RESOURCE", "cpu"))


def _procs():
    return int(os.environ.get("IGC_RECODON_PROCS", "0")) or None


_ATTACHED = weakref.WeakKeyDictionary()     # ReCodonGeneconv or JointAnalysis -> settings
_POOLS = weakref.WeakKeyDictionary()        # ReCodonGeneconv or JointAnalysis -> RecodonPool
_LEGACY = {}


def attach(obj, processes=None, threads=None, resource=None):
    """Use BEAGLE with exact gradients for this ReCodonGeneconv or JointAnalysis object,
    whatever IGC_RECODON_ENGINE says. ``processes`` > 1 (or any JointAnalysis) uses a resident
    ``RecodonPool``; otherwise a ReCodonGeneconv evaluates in this process."""
    detach(obj)
    s = _settings()
    s.pop("engine")
    if threads is not None:
        s["threads"] = threads
    if resource is not None:
        s["resource"] = resource
    s["processes"] = processes
    _ATTACHED[obj] = s
    _patch()
    return _pool_for(obj)


def detach(obj):
    _ATTACHED.pop(obj, None)
    pool = _POOLS.pop(obj, None)
    if pool is not None:
        pool.close()
    st = _STATE.pop(obj, None)
    if st is not None:
        st.bs.inst.close()


_AVAILABLE = []


def available():
    """Whether the BEAGLE action library loads (checked once per process). When it does not,
    the default engine falls back to the legacy code with one warning."""
    if not _AVAILABLE:
        try:
            get_lib()
            _AVAILABLE.append(True)
        except Exception as e:  # missing or broken BEAGLE install
            import warnings
            warnings.warn("BEAGLE action library not available (%s); IS-IGC uses the legacy jsonctmctree "
                          "code. Set IGC_BEAGLE_PREFIX, or IGC_RECODON_ENGINE=legacy to silence this." % e)
            _AVAILABLE.append(False)
    return _AVAILABLE[0]


def _active(obj):
    if obj in _ATTACHED:
        return True
    return _settings()["engine"] == "beagle" and available()


def _pool_for(obj):
    """The object's pool (created on first use), or None for in-process evaluation."""
    pool = _POOLS.get(obj)
    if pool is not None:
        return pool
    s = dict(_ATTACHED.get(obj) or {})
    if not s:
        s = _settings()
        s.pop("engine")
        s["processes"] = _procs()
    joint = hasattr(obj, "geneconv_list")
    k = s.pop("processes")
    if joint:
        objs = obj.geneconv_list
        k = k or len(objs)
    else:
        objs = [obj]
        if not k or k <= 1:
            return None
    pool = RecodonPool(objs, processes=k, **s)
    _POOLS[obj] = pool
    return pool


def _single(obj, grad):
    """(lnL, g) for a ReCodonGeneconv at obj.x through its pool or in process."""
    obj.update_by_x()
    pool = _pool_for(obj)
    if pool is None:
        s = _ATTACHED.get(obj) or {}
        kw = {k: v for k, v in s.items() if k in ("threads", "resource") and v is not None}
        ll, g, _ = evaluate(obj, grad=grad, **kw)
        return ll, g
    ll, gs = pool.evaluate([obj.x], grad=grad)
    return float(ll[0]), (gs[0] if grad else None)


def _joint(self, x, grad=True):
    """JointAnalysis (f, g) at x with the legacy sign and layout, from its resident pool."""
    self.update_by_x(x)
    pool = _pool_for(self)
    lls, gs = pool.evaluate([g.x for g in self.geneconv_list], grad=grad)
    self.ll = list(-lls)
    for o, l in zip(self.geneconv_list, lls):
        o.ll = float(l)
    f = -float(np.sum(lls))
    if not grad:
        return f, None
    sh = self.shared_parameters
    gs = [-np.asarray(v) for v in gs]
    uniq = np.concatenate([[v[i] for i in range(len(v)) if i not in sh] for v in gs])
    shared = [[v[i] for i in range(len(v)) if i in sh] for v in gs]
    return f, np.concatenate((uniq, np.sum(shared, axis=0)))


def _patch():
    from IGCexpansion.CodonGeneconv import ReCodonGeneconv
    from IGCexpansion.JointAnalysis import JointAnalysis

    if _LEGACY:
        return
    _LEGACY["loglikelihood_and_gradient"] = ReCodonGeneconv.loglikelihood_and_gradient
    _LEGACY["_loglikelihood"] = ReCodonGeneconv._loglikelihood
    _LEGACY["joint_objective_and_gradient"] = JointAnalysis.objective_and_gradient
    _LEGACY["joint_objective_and_gradient_multi_threaded"] = JointAnalysis.objective_and_gradient_multi_threaded

    def loglikelihood_and_gradient(self, display=False):
        if not _active(self):
            return _LEGACY["loglikelihood_and_gradient"](self, display)
        ll, g = _single(self, True)
        self.ll = ll
        if display:
            k = len(self.edge_to_blen)
            print("log likelihood = ", ll)
            print("Edge derivatives = ", g[len(g) - k:])
            print("other derivatives:", g[: len(g) - k])
            print("Current x array = ", self.x)
            if self.Force is not None:
                print("Forced parameter:", self.Force)
        return -ll, -g

    def _loglikelihood(self, store=True, edge_derivative=False):
        if not _active(self):
            return _LEGACY["_loglikelihood"](self, store, edge_derivative)
        if store:
            self.scene_ll = self.get_scene()
        ll, g = _single(self, edge_derivative)
        self.ll = ll
        return ll, (list(g[len(g) - len(self.edge_to_blen):]) if edge_derivative else [])

    def _joint_wrapper(name, parallel):
        def objective(self, x):
            if not _active(self):
                return _LEGACY[name](self, x)
            f, g = _joint(self, x)
            self.auto_save += 1
            if self.auto_save == JointAnalysis.auto_save_step:
                self.save_x()
                self.auto_save = 0
            print("log likelihood = ", f)
            print("Current x array = ", self.x)
            print(("Gradient = " if parallel else "Derivatives = "), g)
            return f, g
        return objective

    ReCodonGeneconv.loglikelihood_and_gradient = loglikelihood_and_gradient
    ReCodonGeneconv._loglikelihood = _loglikelihood
    JointAnalysis.objective_and_gradient = _joint_wrapper("joint_objective_and_gradient", False)
    JointAnalysis.objective_and_gradient_multi_threaded = _joint_wrapper(
        "joint_objective_and_gradient_multi_threaded", True)


def install():
    """Patch ReCodonGeneconv and JointAnalysis; objects use BEAGLE when IGC_RECODON_ENGINE=beagle
    or after attach(obj)."""
    _patch()


def uninstall():
    from IGCexpansion.CodonGeneconv import ReCodonGeneconv
    from IGCexpansion.JointAnalysis import JointAnalysis

    for k, f in _LEGACY.items():
        if k.startswith("joint_"):
            setattr(JointAnalysis, k[len("joint_"):], f)
        else:
            setattr(ReCodonGeneconv, k, f)
    _LEGACY.clear()
