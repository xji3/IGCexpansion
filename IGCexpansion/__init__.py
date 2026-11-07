# Shims for jsonctmctree on current NumPy/SciPy; must run before any submodule imports it.
from . import _jct_compat  # noqa: F401
