"""ctypes binding to libhmsbeagle with the action plugin (only the calls the action path needs).

The install prefix is ``$IGC_BEAGLE_PREFIX`` or ``~/opt/beagle-action-07fc8b6`` (the pinned
beagle-dev ``action`` commit). Flag values are parsed from the installed ``beagle.h``, because
they differ between branches (COMPUTATION_ACTION is bit 40 on old forks and bit 42 here).

BEAGLE does no index checking and is not safe to drive from several Python threads (normest1
uses a shared, unseeded engine), and ctypes releases the GIL. Every call sequence therefore
holds ``LOCK``; parallelism is by processes only.
"""

import ctypes as C
import glob
import os
import re
import threading

import numpy as np

NONE = -1
DP = C.POINTER(C.c_double)
IP = C.POINTER(C.c_int)

DEFAULT_PREFIX = os.path.join(os.path.expanduser("~"), "opt", "beagle-action-07fc8b6")
LOCK = threading.RLock()


class BeagleError(RuntimeError):
    def __init__(self, func, code, name=None):
        self.func, self.code, self.name = func, code, name
        super().__init__("%s returned %d%s" % (func, code, " (%s)" % name if name else ""))


class BeagleInstanceDetails(C.Structure):
    _fields_ = [("resourceNumber", C.c_int), ("resourceName", C.c_char_p), ("implName", C.c_char_p),
                ("implDescription", C.c_char_p), ("flags", C.c_longlong)]


class BeagleResource(C.Structure):
    _fields_ = [("name", C.c_char_p), ("description", C.c_char_p), ("supportFlags", C.c_longlong),
                ("requiredFlags", C.c_longlong)]


class BeagleResourceList(C.Structure):
    _fields_ = [("list", C.POINTER(BeagleResource)), ("length", C.c_int)]


class BeagleOperation(C.Structure):
    _fields_ = [(n, C.c_int) for n in ("destinationPartials", "destinationScaleWrite", "destinationScaleRead",
                                       "child1Partials", "child1TransitionMatrix",
                                       "child2Partials", "child2TransitionMatrix")]


def dptr(a):
    return a.ctypes.data_as(DP)


def iptr(a):
    return a.ctypes.data_as(IP)


def i32(x):
    return np.ascontiguousarray(x, dtype=np.int32)


def f64(x):
    return np.ascontiguousarray(x, dtype=np.float64)


def parse_header(path):
    """Flags (enum entries and the #define'd ones) and return codes from beagle.h."""
    txt = open(path).read()
    flags = {}
    for m in re.finditer(r"^\s*BEAGLE_FLAG_(\w+)\s*=\s*\(?\s*1\s*(?:U?LL)?\s*<<\s*(\d+)", txt, re.M):
        flags[m.group(1)] = 1 << int(m.group(2))
    for m in re.finditer(r"#define\s+BEAGLE_FLAG_(\w+)\s+\(\s*1\s*(?:U?LL)?\s*<<\s*(\d+)\s*\)", txt):
        flags[m.group(1)] = 1 << int(m.group(2))
    codes = {int(m.group(2)): m.group(1)
             for m in re.finditer(r"^\s*BEAGLE_(SUCCESS|ERROR_\w+)\s*=\s*(-?\d+)", txt, re.M)}
    return flags, codes


_SIGNATURES = {
    "beagleSetTipPartials": [C.c_int, C.c_int, DP],
    "beagleSetPartials": [C.c_int, C.c_int, DP],
    "beagleGetPartials": [C.c_int, C.c_int, C.c_int, DP],
    "beagleSetPatternWeights": [C.c_int, DP],
    "beagleSetCategoryRates": [C.c_int, DP],
    "beagleSetCategoryWeights": [C.c_int, C.c_int, DP],
    "beagleSetStateFrequencies": [C.c_int, C.c_int, DP],
    "beagleSetSparseMatrix": [C.c_int, C.c_int, IP, IP, DP, C.c_int],
    "beagleSetSparseDifferentialMatrix": [C.c_int, C.c_int, IP, IP, DP, C.c_int],
    "beagleUpdateTransitionMatrices": [C.c_int, C.c_int, IP, IP, IP, DP, C.c_int],
    "beagleUpdatePartials": [C.c_int, C.POINTER(BeagleOperation), C.c_int, C.c_int],
    "beagleUpdatePrePartials": [C.c_int, C.POINTER(BeagleOperation), C.c_int, C.c_int],
    "beagleResetScaleFactors": [C.c_int, C.c_int],
    "beagleAccumulateScaleFactors": [C.c_int, IP, C.c_int, C.c_int],
    "beagleCalculateRootLogLikelihoods": [C.c_int, IP, IP, IP, IP, C.c_int, DP],
    "beagleGetSiteLogLikelihoods": [C.c_int, DP],
    "beagleCalculateEdgeDerivatives": [C.c_int, IP, IP, IP, IP, C.c_int, DP, DP, DP],
    "beagleSetCPUThreadCount": [C.c_int, C.c_int],
    "beagleFinalizeInstance": [C.c_int],
}


class BeagleLib:
    """One loaded libhmsbeagle (per process)."""

    def __init__(self, prefix=None, preload_plugin=None):
        self.prefix = prefix = os.path.abspath(prefix or os.environ.get("IGC_BEAGLE_PREFIX") or DEFAULT_PREFIX)
        libdir = os.path.join(prefix, "lib")
        darwin = os.uname().sysname == "Darwin"
        cands = [os.path.join(libdir, "libhmsbeagle.dylib")] if darwin else \
            [os.path.join(libdir, "libhmsbeagle.so")] + sorted(glob.glob(os.path.join(libdir, "libhmsbeagle.so.*")))
        path = next((p for p in cands if os.path.exists(p)), None)
        if path is None:
            raise OSError("libhmsbeagle not found under %s (set IGC_BEAGLE_PREFIX)" % libdir)
        self.path = path
        self.lib = L = C.CDLL(path, mode=C.RTLD_GLOBAL)
        # BEAGLE dlopens plugins by bare leaf name and a failed load is silent; on a relocated
        # prefix the absolute-path preload makes the action plugin available anyway.
        if preload_plugin is None:
            preload_plugin = os.environ.get("IGC_BEAGLE_PRELOAD", "0") == "1"
        self.preloaded = []
        if preload_plugin:
            for p in sorted(glob.glob(os.path.join(libdir, "libhmsbeagle-cpu-action*"))):
                if re.search(r"\.\d+\.so$|\.so(\.\d+)+$", p):
                    C.CDLL(p, mode=C.RTLD_GLOBAL)
                    self.preloaded.append(p)
        header = os.path.join(prefix, "include", "libhmsbeagle-1", "libhmsbeagle", "beagle.h")
        self.F, self.codes = parse_header(header)
        if "COMPUTATION_ACTION" not in self.F:
            raise OSError("%s has no BEAGLE_FLAG_COMPUTATION_ACTION: not an action-enabled BEAGLE" % header)
        L.beagleCreateInstance.argtypes = [C.c_int] * 9 + [IP, C.c_int, C.c_longlong, C.c_longlong,
                                                           C.POINTER(BeagleInstanceDetails)]
        L.beagleCreateInstance.restype = C.c_int
        L.beagleGetResourceList.argtypes = []
        L.beagleGetResourceList.restype = C.POINTER(BeagleResourceList)
        for name, args in _SIGNATURES.items():
            fn = getattr(L, name)
            fn.argtypes = args
            fn.restype = C.c_int
        with LOCK:
            rl = L.beagleGetResourceList().contents
            self.resources = [(i, rl.list[i].name.decode(), int(rl.list[i].supportFlags)) for i in range(rl.length)]

    def call(self, name, *args):
        rc = getattr(self.lib, name)(*args)
        if rc != 0:
            raise BeagleError(name, rc, self.codes.get(rc))
        return rc

    def gpu_resources(self):
        return [(i, n) for i, n, f in self.resources if f & self.F["PROCESSOR_GPU"]]


_LIBS = {}


def get_lib(prefix=None):
    """The process-wide BeagleLib for ``prefix`` (default: $IGC_BEAGLE_PREFIX or DEFAULT_PREFIX)."""
    key = os.path.abspath(prefix or os.environ.get("IGC_BEAGLE_PREFIX") or DEFAULT_PREFIX)
    with LOCK:
        if key not in _LIBS:
            _LIBS[key] = BeagleLib(key)
        return _LIBS[key]
