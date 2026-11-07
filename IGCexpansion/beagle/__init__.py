"""BEAGLE action engine for IGCexpansion (likelihoods and edge derivatives).

Enable with ``IGC_ENGINE=beagle`` (or ``shadow`` to cross-check against jsonctmctree) before
importing IGCexpansion, or call ``IGCexpansion.beagle.install('beagle')``. jsonctmctree remains
the default engine, the oracle, and the fallback for anything this engine does not handle.
"""

from .engine import STATS, SHADOW_LOG, SHADOW_MAX, clear_cache, get_mode, install, process_json_in, set_mode, uninstall  # noqa: F401
