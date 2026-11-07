"""L0: the BEAGLE binding and instance wrapper on their own (plan §11).

* The actiontest example from beagle-lib (3 taxa, 5 sites, 2 rate categories): logL matches
  its R reference, and the edge derivatives (2 categories) match a dense scipy expm.
* Slot isolation: a derivative matrix in a spare slot leaves lnL bitwise unchanged.
* The checks BEAGLE does not make: implementation and flags, index ranges, epochs, ownership.
"""

import os
import pickle
import sys
import threading

import numpy as np
import pytest
import scipy.linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from beagle_util import lib_or_skip  # noqa: E402

from IGCexpansion.beagle._lib import NONE  # noqa: E402
from IGCexpansion.beagle.instance import ActionInstance, EpochError  # noqa: E402

SEQ = {0: "GAGTC", 1: "GAGGT", 2: "AAATC"}  # human, chimp, gorilla
RATES = np.array([0.14251623900062188, 1.857483760999378])
FREQS = np.array([0.1, 0.3, 0.2, 0.4])
QROWS = np.array([[-1.285714, 0.4285712, 0.2857142, 0.5714284],
                  [0.142857, -0.9999997, 0.2857143, 0.5714284],
                  [0.142857, 0.4285714, -1.1428568, 0.5714284],
                  [0.142857, 0.4285713, 0.2857142, -0.8571426]])
T = {0: 0.6, 1: 0.6, 2: 1.3, 3: 0.7}       # node -> length of the edge above it
R_REFERENCE = -21.827146738233282


def _tip(s):
    a = np.zeros((len(s), 4))
    for i, ch in enumerate(s):
        a[i, "ACGT".index(ch)] = 1.0
    return a


def _actiontest(lib, spare=0, threads=1):
    inst = ActionInstance(lib, tips=3, partials=20, states=4, patterns=5, generators=1, spare_slots=spare,
                          keys=12, categories=2, scale_buffers=7, threads=threads)
    for n, s in SEQ.items():
        inst.set_tip_partials(n, _tip(s))
    inst.set_categories(RATES, [0.5, 0.5])
    inst.set_pattern_weights(np.ones(5))
    inst.set_state_frequencies(FREQS)
    r, c = np.nonzero(np.ones((4, 4)))
    inst.set_generator(0, r, c, QROWS[r, c])
    inst.update_transitions(0, [0, 1, 2, 3], [T[0], T[1], T[2], T[3]])
    return inst


def _post(inst):
    inst.update_partials([(3, 0, NONE, 0, 0, 1, 1), (4, 1, NONE, 2, 2, 3, 3)])
    inst.accumulate_scaling([0, 1], 2)
    return inst.root_log_likelihood(4, 2)


def _pre(inst):
    inst.set_partials(5, np.tile(FREQS, 2 * 5))
    inst.update_pre_partials([(6, NONE, NONE, 5, 3, 2, 2), (7, NONE, NONE, 5, 2, 3, 3),
                              (8, NONE, NONE, 6, 1, 0, 0), (9, NONE, NONE, 6, 0, 1, 1)])


def _dense(dnode=None):
    """Site likelihoods (and d/dt of the edge above dnode) for the actiontest tree, by expm."""
    tips = {n: _tip(s).T for n, s in SEQ.items()}       # 4 x 5
    L = np.zeros(5)
    dL = np.zeros(5)
    for r, w in zip(RATES, [0.5, 0.5]):
        P = {n: scipy.linalg.expm(QROWS * r * T[n]) for n in T}
        dP = {n: r * QROWS @ P[n] for n in T}

        def root_vec(d):
            def branch(n, below):
                return (dP[n] if n == d else P[n]) @ below
            n3 = branch(0, tips[0]) * branch(1, tips[1])
            return branch(2, tips[2]) * branch(3, n3)
        L += w * (FREQS @ root_vec(None))
        dL += w * (FREQS @ root_vec(dnode))
    return L, dL


def test_actiontest_logl_matches_r():
    lib = lib_or_skip()
    inst = _actiontest(lib)
    total, site = _post(inst)
    assert total == pytest.approx(R_REFERENCE, rel=1e-10)
    L, _ = _dense()
    np.testing.assert_allclose(site, np.log(L), rtol=1e-12)


@pytest.mark.parametrize("node", [0, 1, 2, 3])
def test_two_category_edge_derivatives_match_dense(node):
    lib = lib_or_skip()
    inst = _actiontest(lib)
    _post(inst)
    _pre(inst)
    post = {0: 0, 1: 1, 2: 2, 3: 3}[node]
    pre = {0: 9, 1: 8, 2: 7, 3: 6}[node]       # actiontest: pre 8 = node 1 (key 1), pre 9 = node 0
    tot, per = inst.edge_derivatives([post], [pre], [0], per_pattern=True)
    L, dL = _dense(node)
    np.testing.assert_allclose(per[0], dL / L, rtol=1e-10, atol=1e-12)
    assert tot[0] == pytest.approx(np.sum(dL / L), rel=1e-10, abs=1e-12)


def test_spare_slot_isolation_and_same_derivatives():
    lib = lib_or_skip()
    inst = _actiontest(lib, spare=4)
    a, sa = _post(inst)
    r, c = np.nonzero(np.ones((4, 4)))
    inst.set_differential(4, r, c, QROWS[r, c] * 3.0)      # something different from Q
    b, sb = _post(inst)
    assert a == b and np.array_equal(sa, sb)               # bitwise
    inst.set_differential(4, r, c, QROWS[r, c])
    _pre(inst)
    t0, _ = inst.edge_derivatives([0, 1, 2, 3], [9, 8, 7, 6], [0, 0, 0, 0])
    t4, _ = inst.edge_derivatives([0, 1, 2, 3], [9, 8, 7, 6], [4, 4, 4, 4])
    np.testing.assert_allclose(t4, t0, rtol=1e-13)


def test_flags_and_threads():
    lib = lib_or_skip()
    F = lib.F
    inst = _actiontest(lib)
    assert inst.impl == "CPU-Action-Double"
    assert inst.flags & F["COMPUTATION_ACTION"] and inst.flags & F["SCALING_MANUAL"]
    assert not inst.flags & F["THREADING_CPP"]
    a, _ = _post(inst)
    inst2 = _actiontest(lib, threads=2)
    assert inst2.flags & F["THREADING_CPP"]
    b, _ = _post(inst2)
    assert b == pytest.approx(a, rel=1e-14)


def test_range_checks():
    lib = lib_or_skip()
    inst = _actiontest(lib, spare=1)
    r = c = np.arange(4)
    with pytest.raises(IndexError):
        inst.set_generator(1, r, c, np.ones(4))              # spare slot is not a generator slot
    with pytest.raises(IndexError):
        inst.set_differential(0, r, c, np.ones(4))           # never overwrite a generator
    with pytest.raises(IndexError):
        inst.set_differential(2, r, c, np.ones(4))           # beyond eigenCount
    with pytest.raises(IndexError):
        inst.set_generator(0, [0, 4], [0, 0], [1.0, 1.0])    # state out of range
    with pytest.raises(IndexError):
        inst.update_transitions(0, [12], [0.1])              # key out of range
    with pytest.raises(ValueError):
        inst.update_transitions(0, [1], [-0.1])
    with pytest.raises(IndexError):
        inst.set_tip_partials(3, np.ones(20))
    with pytest.raises(IndexError):
        inst.update_partials([(20, NONE, NONE, 0, 0, 1, 1)])
    with pytest.raises(IndexError):
        inst.update_partials([(3, 7, NONE, 0, 0, 1, 1)])     # scale buffer out of range
    with pytest.raises(ValueError):
        inst.set_tip_partials(0, np.ones(19))


def test_epochs():
    lib = lib_or_skip()
    inst = _actiontest(lib)
    with pytest.raises(EpochError):
        inst.root_log_likelihood(4, 2)                       # no post pass yet
    _post(inst)
    with pytest.raises(EpochError):
        inst.edge_derivatives([0], [8], [0])                 # no pre pass yet
    _pre(inst)
    inst.edge_derivatives([0], [9], [0])
    inst.update_transitions(0, [0], [0.61])                  # values changed
    with pytest.raises(EpochError):
        inst.edge_derivatives([0], [9], [0])
    with pytest.raises(EpochError):
        _pre(inst)


def test_ownership():
    lib = lib_or_skip()
    inst = _actiontest(lib)
    with pytest.raises(TypeError):
        pickle.dumps(inst)
    err = []

    def other():
        try:
            inst.set_pattern_weights(np.ones(5))
        except RuntimeError as e:
            err.append(e)
    th = threading.Thread(target=other)
    th.start()
    th.join()
    assert err
    inst.close()
    with pytest.raises(RuntimeError):
        inst.set_pattern_weights(np.ones(5))
