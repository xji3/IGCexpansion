"""Shared helpers for the BEAGLE tests: skip cleanly when no action-enabled BEAGLE is installed."""

import pytest


def lib_or_skip():
    try:
        from IGCexpansion.beagle._lib import get_lib

        return get_lib()
    except OSError as e:
        pytest.skip("BEAGLE not available: %s" % e)
