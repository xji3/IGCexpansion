"""ActionInstance: one BEAGLE action instance with the checks BEAGLE itself does not make.

* Implementation: the instance must be ``CPU-Action-Double`` (CPU) or ``Action-CUDA-Double``
  (GPU); on CPU it must also report COMPUTATION_ACTION and SCALING_MANUAL. (The CUDA action
  instance does not report the ACTION flag, so it is checked by name.)
* Threads: THREADING_CPP is requested only for ``threads > 1``, and ``setCPUThreadCount`` is
  called right after creation, before any data call.
* Indices: BEAGLE does not range-check buffer, key, slot or scale indices; every call here does.
* Slots: eigen slots ``[0, generators)`` hold generators (``set_generator``); the rest are spare
  slots for derivative matrices (``set_differential``) and never overwrite a generator.
* Epochs: post-order results are valid only for the current values, pre-order only after a post
  pass at those values, and edge derivatives only after both.
* Ownership: an instance is used only by the process and thread that created it (no pickling;
  a forked child must build its own), and is finalized when garbage-collected.
"""

import os
import threading
import weakref

import numpy as np

from ._lib import LOCK, NONE, BeagleInstanceDetails, BeagleOperation, dptr, f64, i32, iptr
import ctypes as C

IMPL = {"cpu": "CPU-Action-Double", "gpu": "Action-CUDA-Double"}


class EpochError(RuntimeError):
    pass


def _finalize(lib, handle, pid):
    if os.getpid() == pid:
        with LOCK:
            lib.lib.beagleFinalizeInstance(handle)


class ActionInstance:
    def __init__(self, lib, *, tips, partials, states, patterns, generators, spare_slots=0, keys,
                 categories=1, scale_buffers=0, resource="cpu", device=None, threads=1):
        if resource not in IMPL:
            raise ValueError("resource must be 'cpu' or 'gpu'")
        F = lib.F
        self.lib, self.resource = lib, resource
        self.tips, self.partials, self.states, self.patterns = tips, partials, states, patterns
        self.generators, self.eigen = generators, generators + spare_slots
        self.categories, self.scale_buffers, self.threads = categories, scale_buffers, int(threads)
        self.keys = keys
        pref = F["COMPUTATION_ACTION"] | F["PRECISION_DOUBLE"] | F["SCALING_MANUAL"] | F["SCALERS_LOG"]
        rl, nrl = None, 0
        if resource == "gpu":
            req = F["PROCESSOR_GPU"] | F["FRAMEWORK_CUDA"]
            gpus = lib.gpu_resources()
            if not gpus and device is None:
                raise RuntimeError("no GPU resource in this BEAGLE build/host")
            rl, nrl = (C.c_int * 1)(gpus[0][0] if device is None else device), 1
            matrix_count = keys          # GPU transition keys index matrix buffers
        else:
            req = F["PROCESSOR_CPU"]
            if self.threads > 1:
                pref |= F["THREADING_CPP"]
            if keys > 2 * partials:
                raise ValueError("CPU action keys index the 2*partials key map: %d > %d" % (keys, 2 * partials))
            matrix_count = 0             # CPU action stores no dense matrices
        info = BeagleInstanceDetails()
        with LOCK:
            h = lib.lib.beagleCreateInstance(tips, partials, 0, states, patterns, self.eigen, matrix_count,
                                             categories, scale_buffers, rl, nrl, pref, req, C.byref(info))
            if h < 0:
                raise RuntimeError("beagleCreateInstance failed: %d (%s)" % (h, lib.codes.get(h)))
            self.impl = info.implName.decode()
            self.resource_name = info.resourceName.decode()
            self.flags = int(info.flags)
            ok = self.impl == IMPL[resource]
            if resource == "cpu":
                ok = ok and bool(self.flags & F["COMPUTATION_ACTION"]) and bool(self.flags & F["SCALING_MANUAL"])
            ok = ok and not (self.flags & F["PRECISION_SINGLE"])
            if not ok:
                lib.lib.beagleFinalizeInstance(h)
                raise RuntimeError("wrong BEAGLE implementation: %s (flags %#x), expected %s"
                                   % (self.impl, self.flags, IMPL[resource]))
            if resource == "cpu" and self.threads > 1:
                lib.call("beagleSetCPUThreadCount", h, self.threads)
        self.handle = h
        self._pid, self._tid = os.getpid(), threading.get_ident()
        self._finalizer = weakref.finalize(self, _finalize, lib, h, self._pid)
        self.values_epoch, self.post_epoch, self.pre_epoch = 0, -1, -1

    # ------------------------------------------------------------------ guards
    def __getstate__(self):
        raise TypeError("ActionInstance holds a native BEAGLE handle and cannot be pickled")

    def close(self):
        self._finalizer()

    def _own(self):
        if os.getpid() != self._pid:
            raise RuntimeError("BEAGLE instance used in a different process than the one that created it")
        if threading.get_ident() != self._tid:
            raise RuntimeError("BEAGLE instance used from a different thread than the one that created it")
        if not self._finalizer.alive:
            raise RuntimeError("BEAGLE instance already finalized")

    def _call(self, name, *args):
        self._own()
        with LOCK:
            return self.lib.call(name, self.handle, *args)

    @staticmethod
    def _range(what, idx, hi, lo=0):
        idx = np.atleast_1d(np.asarray(idx))
        if idx.size and (idx.min() < lo or idx.max() >= hi):
            raise IndexError("%s index out of range [%d, %d): %s" % (what, lo, hi, idx[(idx < lo) | (idx >= hi)][:5]))

    def _len(self, what, a, n):
        if a.size != n:
            raise ValueError("%s: expected %d values, got %d" % (what, n, a.size))

    def _values_changed(self):
        self.values_epoch += 1

    # ------------------------------------------------------------------ data
    def set_tip_partials(self, tip, a):
        self._range("tip", tip, self.tips)
        a = f64(a)
        self._len("tip partials", a, self.patterns * self.states)
        self._call("beagleSetTipPartials", int(tip), dptr(a))
        self._values_changed()

    def set_partials(self, buf, a):
        """Full buffer: categories x patterns x states (used for the ONES buffer and the root pre)."""
        self._range("partials buffer", buf, self.partials, lo=self.tips)
        a = f64(a)
        self._len("partials", a, self.categories * self.patterns * self.states)
        self._call("beagleSetPartials", int(buf), dptr(a))

    def get_partials(self, buf):
        self._range("partials buffer", buf, 2 * self.partials)
        out = np.zeros(self.categories * self.patterns * self.states)
        self._call("beagleGetPartials", int(buf), NONE, dptr(out))
        return out.reshape(self.categories, self.patterns, self.states)

    def set_pattern_weights(self, w):
        w = f64(w)
        self._len("pattern weights", w, self.patterns)
        self._call("beagleSetPatternWeights", dptr(w))
        self._values_changed()

    def set_categories(self, rates, weights):
        r, w = f64(rates), f64(weights)
        self._len("category rates", r, self.categories)
        self._len("category weights", w, self.categories)
        self._call("beagleSetCategoryRates", dptr(r))
        self._call("beagleSetCategoryWeights", 0, dptr(w))
        self._values_changed()

    def set_state_frequencies(self, pi):
        pi = f64(pi)
        self._len("state frequencies", pi, self.states)
        self._call("beagleSetStateFrequencies", 0, dptr(pi))
        self._values_changed()

    def _sparse(self, name, slot, r, c, v):
        r, c, v = i32(r), i32(c), f64(v)
        if not (len(r) == len(c) == len(v)):
            raise ValueError("triplet lengths differ")
        self._range("sparse row", r, self.states)
        self._range("sparse column", c, self.states)
        self._call(name, int(slot), iptr(r), iptr(c), dptr(v), len(v))

    def set_generator(self, slot, r, c, v):
        """Generator Q for process ``slot`` (diagonal included by the caller; BEAGLE adds none)."""
        self._range("generator slot", slot, self.generators)
        self._sparse("beagleSetSparseMatrix", slot, r, c, v)
        self._values_changed()

    def set_differential(self, slot, r, c, v):
        """Derivative matrix in a spare slot; leaves every generator and all likelihood state untouched."""
        self._range("spare slot", slot, self.eigen, lo=self.generators)
        self._sparse("beagleSetSparseDifferentialMatrix", slot, r, c, v)

    def update_transitions(self, generator, keys, lengths):
        """e^{Q_g t} for each key (keys are an index space of their own, not buffers)."""
        self._range("generator slot", generator, self.generators)
        keys, lengths = i32(keys), f64(lengths)
        if len(keys) == 0:
            return
        self._range("transition key", keys, self.keys)
        if len(lengths) != len(keys) or np.any(~np.isfinite(lengths)) or np.any(lengths < 0):
            raise ValueError("edge lengths must be finite and >= 0, one per key")
        self._call("beagleUpdateTransitionMatrices", int(generator), iptr(keys), None, None, dptr(lengths), len(keys))
        self._values_changed()

    # ------------------------------------------------------------------ passes
    def _ops(self, ops):
        ops = np.asarray(ops, dtype=np.int64).reshape(-1, 7)
        self._range("destination buffer", ops[:, 0], self.partials, lo=self.tips)
        self._range("child buffer", ops[:, [3, 5]], self.partials)
        self._range("transition key", ops[:, [4, 6]], self.keys)
        sw = ops[:, 1][ops[:, 1] != NONE]
        self._range("scale buffer", sw, self.scale_buffers)
        if np.any(ops[:, 2] != NONE):
            raise ValueError("read-scale indices are not used on the action path")
        return (BeagleOperation * len(ops))(*[BeagleOperation(*map(int, o)) for o in ops])

    def update_partials(self, ops):
        """Post-order ops (dest, writeScale, NONE, child1, key1, child2, key2)."""
        arr = self._ops(ops)
        self._call("beagleUpdatePartials", arr, len(arr), NONE)
        self.post_epoch = self.values_epoch

    def update_pre_partials(self, ops):
        """Pre-order ops (pre(n), NONE, NONE, pre(parent), key(n), post(sibling), key(sibling))."""
        if self.post_epoch != self.values_epoch:
            raise EpochError("pre-order pass needs a post-order pass at the current values")
        arr = self._ops(ops)
        self._call("beagleUpdatePrePartials", arr, len(arr), NONE)
        self.pre_epoch = self.post_epoch

    def accumulate_scaling(self, scale_indices, cumulative):
        self._range("scale buffer", list(scale_indices) + [cumulative], self.scale_buffers)
        s = i32(scale_indices)
        self._call("beagleResetScaleFactors", int(cumulative))
        if len(s):
            self._call("beagleAccumulateScaleFactors", iptr(s), len(s), int(cumulative))

    def root_log_likelihood(self, root, cumulative=NONE):
        """Total weighted lnL and per-pattern site lnL at ``root``."""
        if self.post_epoch != self.values_epoch:
            raise EpochError("root likelihood needs a post-order pass at the current values")
        self._range("root buffer", root, self.partials)
        if cumulative != NONE:
            self._range("scale buffer", cumulative, self.scale_buffers)
        out = np.zeros(1)
        r, z, cc = i32([root]), i32([0]), i32([cumulative])
        self._call("beagleCalculateRootLogLikelihoods", iptr(r), iptr(z), iptr(z), iptr(cc), 1, dptr(out))
        site = np.zeros(self.patterns)
        self._call("beagleGetSiteLogLikelihoods", dptr(site))
        return float(out[0]), site

    def edge_derivatives(self, post, pre, deriv_slots, per_pattern=False):
        """sum_s w_s (q' e^{Qt} D p)/(q' e^{Qt} p) for each (post, pre, D) triple; optionally per pattern."""
        if not (self.pre_epoch == self.post_epoch == self.values_epoch):
            raise EpochError("edge derivatives need post- and pre-order passes at the current values")
        post, pre, d = i32(post), i32(pre), i32(deriv_slots)
        n = len(post)
        if not (len(pre) == len(d) == n):
            raise ValueError("post, pre and slot lists differ in length")
        self._range("post buffer", post, self.partials)
        self._range("pre buffer", pre, self.partials)
        self._range("derivative slot", d, self.eigen)
        tot = np.zeros(n)
        per = np.zeros(n * self.patterns) if per_pattern else None
        cw = i32([0])
        self._call("beagleCalculateEdgeDerivatives", iptr(post), iptr(pre), iptr(d), iptr(cw), n,
                   dptr(per) if per_pattern else None, dptr(tot), None)
        return tot, (per.reshape(n, self.patterns) if per_pattern else None)
