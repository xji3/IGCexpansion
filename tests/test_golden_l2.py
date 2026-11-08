"""L2: jsonctmctree log-likelihoods at fixed points match the checked-in golden values.

Fast cases run by default; ``-m slow`` runs the slow ones (3721-state YGL062W, which needs the
IGCTract alignments or IGC_MAFFT_DIR, and the full PS-IGC objective); ``-m ""`` runs all.
"""

import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases  # noqa: E402

GOLDEN = json.load(open(os.path.join(HERE, "golden", "l2.json")))["values"]
RTOL = 1e-9


def _params():
    for name in cases.CASES:
        marks = []
        if name in cases.SLOW:
            marks.append(pytest.mark.slow)
        if name in cases.EXTERNAL:
            marks.append(pytest.mark.external)
        yield pytest.param(name, marks=marks, id=name)


@pytest.mark.parametrize("name", list(_params()))
def test_golden(name, tmp_path):
    try:
        got = cases.evaluate(name, str(tmp_path))
    except FileNotFoundError as e:
        pytest.skip("missing data: %s" % e)
    want = GOLDEN[name]
    assert set(got) == set(want)
    for k in want:
        assert got[k] == pytest.approx(want[k], rel=RTOL), k


@pytest.mark.parametrize("name", sorted(cases.SPIKE))
def test_golden_matches_spike(name):
    """The checked-in values agree with the BEAGLE spike's recorded outputs."""
    for k, v in cases.SPIKE[name].items():
        assert GOLDEN[name][k] == pytest.approx(v, rel=1e-12), (name, k)


def test_axis_order_is_hash_independent(tmp_path):
    """Paralog axes follow the given paralog order, not set() iteration order."""
    with cases.quiet():
        obj = cases.hky_ylr406c(str(tmp_path))
    want = [0 if n.endswith("YLR406C") else 1 for n in obj.observable_names]
    assert obj.observable_axes == want


def test_compat_shims_applied():
    import numpy.testing
    import scipy.sparse

    import IGCexpansion  # noqa: F401

    assert hasattr(numpy.testing, "Tester")
    assert hasattr(scipy.sparse.spmatrix, "A")
