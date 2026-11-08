# Shims for jsonctmctree on current NumPy/SciPy; must run before any submodule imports it.
from . import _jct_compat  # noqa: F401

import os as _os

if _os.environ.get("IGC_ENGINE", "jsonctmctree") != "jsonctmctree":
    # Route jsonctmctree calls through the BEAGLE engine (process workers inherit this on import).
    from . import beagle as _beagle

    _beagle.install()

if _os.environ.get("IGC_PSJS_ENGINE", "legacy") != "legacy":
    # PS-IGC objective on resident workers with the selected engine (jsonctmctree | beagle).
    from .beagle import psjs as _psjs

    _psjs.install()
