"""L3b: PS-IGC Godambe ingredients with exact (GL) per-pair scores (plan §6.2, §11).

Legacy fixtures test/PSJS_HKY_rv_YDR418W_YEL054C_nonclock_{gradient,hessian}.txt hold per-pair
2-D (ln tau, logit p) gradients and 4th-order FD Hessians on the 6-site test alignment at
x = log(4, 0.3), IGC on. Their rows are in an older dict order, so rows are matched by pair.

Findings pinned here:
* GL per-pair scores reproduce the fixture gradients to ~1e-10.
* GL-based Hessians (central differences of exact scores) are stable to ~1e-9 across steps.
* The legacy Hessians (second differences of lnL) carry errors of the size of the entries
  themselves, so the legacy Godambe matrix is unreliable; the GL one is not.
"""

import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from beagle_util import lib_or_skip  # noqa: E402
import cases  # noqa: E402

X2D = np.log([4.0, 0.3])


@pytest.fixture(scope="module")
def obj():
    lib_or_skip()
    from IGCexpansion.PSJSGeneconv import PSJSGeneconv
    T, d = cases.T, tempfile.mkdtemp()
    x_js = np.log([0.5, 0.5, 0.5, 4.35588244, 1.5, 5.0, 0.3, 1.0 / 30.0])
    with cases.quiet():
        return PSJSGeneconv(T + "/YDR418W_YEL054C_test.fasta", T + "/YDR418W_YEL054C_GeneToOrlg.txt", None, True,
                            True, T + "/YeastTree.newick", T + "/YeastTestDupLost.txt", x_js, "HKY", "One rate", True,
                            {"D1": 0}, ["kluyveri", "castellii", "bayanus", "kudriavzevii", "mikatae", "paradoxus",
                                        "cerevisiae"], os.path.join(d, "s.txt"), os.path.join(d, "l.txt"),
                            force=None, space_list=None)


@pytest.fixture(scope="module")
def gl(obj):
    from IGCexpansion.beagle.psjs import godambe_2d
    return {st: godambe_2d(obj, X2D, step=st) for st in (1e-4, 1e-5)}


def _match(fixture_g, g):
    """Fixture row -> pair: optimal one-to-one assignment on the gradient rows (two pairs have
    nearly identical gradients, so nearest-row matching is not enough)."""
    from scipy.optimize import linear_sum_assignment
    pairs = list(g)
    G = np.array([g[p] for p in pairs])
    cost = np.max(np.abs(fixture_g[:, None, :] - G[None, :, :]), axis=2)
    rows, cols = linear_sum_assignment(cost)
    return [(pairs[c], float(cost[r, c])) for r, c in sorted(zip(rows, cols))]


def test_scores_match_fixture(gl):
    g, _, _ = gl[1e-5]
    fg = np.loadtxt(os.path.join(cases.T, "PSJS_HKY_rv_YDR418W_YEL054C_nonclock_gradient.txt"))
    m = _match(fg, g)
    assert len(g) == len(fg) == 15
    assert len({p for p, _ in m}) == 15                       # a bijection
    assert max(e for _, e in m) < 1e-8


def test_gl_hessian_step_stable(gl):
    (_, h4, G4), (_, h5, G5) = gl[1e-4], gl[1e-5]
    assert max(np.max(np.abs(h4[p] - h5[p])) for p in h5) < 1e-8
    np.testing.assert_allclose(G4, G5, rtol=1e-6, atol=1e-8)
    assert np.all(np.linalg.eigvalsh(G5) > 0)


def test_legacy_hessians_are_noise_limited(gl):
    """Report-level: legacy FD Hessians differ from the GL ones by about their own size."""
    g, h, _ = gl[1e-5]
    fg = np.loadtxt(os.path.join(cases.T, "PSJS_HKY_rv_YDR418W_YEL054C_nonclock_gradient.txt"))
    fh = np.loadtxt(os.path.join(cases.T, "PSJS_HKY_rv_YDR418W_YEL054C_nonclock_hessian.txt"))
    diffs = [np.max(np.abs(fh[i].reshape(2, 2) - h[p])) for i, (p, _) in enumerate(_match(fg, g))]
    typical = np.median(np.abs(fh))
    assert max(diffs) < 1e-4                                  # same quantity ...
    assert max(diffs) > 0.1 * typical                         # ... but legacy noise ~ entry size
