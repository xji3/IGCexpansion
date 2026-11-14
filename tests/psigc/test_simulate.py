import os
import numpy as np
import pytest
from IGCexpansion import psigc
from IGCexpansion.psigc import simulate as S

# bundle of the 2017 PS-IGC simulation study (alignments and truth parameters): part of the paper's
# supplementary folder psigc_reproduce, outside this repository (default: next to the repository)
BUNDLE = os.environ.get('PSIGC_SIM2017_BUNDLE', os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    'psigc_reproduce', 'data', 'sim2017_bundle'))
pytestmark = [pytest.mark.external,
              pytest.mark.skipif(not os.path.isdir(BUNDLE), reason='2017 simulation bundle not found (set PSIGC_SIM2017_BUNDLE)')]


@pytest.fixture(scope='module')
def setup():
    import json
    tree = psigc.Tree(open(os.path.join(BUNDLE, 'YeastTree.newick')).read())
    si = np.loadtxt(os.path.join(BUNDLE, 'seq_index.txt'), dtype=int)
    man = [json.loads(l) for l in open(os.path.join(BUNDLE, 'manifest.jsonl'))]
    return tree, si, man


def par_for(setup, cell):
    tree, si, man = setup
    return psigc.unpack_x(np.array([m for m in man if m['cell'] == cell][0]['truth_x']), tree)


def test_reproducible_streams(setup):
    tree, si, _ = setup
    par = par_for(setup, 50.0)
    a = S.simulate(par, tree, si, S.rng_for(1, 3, 50, 7)[0])
    b = S.simulate(par, tree, si, S.rng_for(1, 3, 50, 7)[0])
    c = S.simulate(par, tree, si, S.rng_for(1, 3, 50, 8)[0])
    assert all(np.array_equal(a.tips[k], b.tips[k]) for k in a.tips)
    assert any(not np.array_equal(a.tips[k], c.tips[k]) for k in a.tips)


def test_shapes_and_tracts(setup):
    tree, si, _ = setup
    par = par_for(setup, 3.0)
    res = S.simulate(par, tree, si, S.rng_for(2, 0, 0, 0)[0])
    assert res.L == 492 and len(res.columns) == 489
    for sp in tree.post_tips:
        assert res.tips[sp].shape == (2, 492)
    for sp in tree.pre_tips:
        assert res.tips[sp].shape == (1, 492)
    for e in res.igc_log:
        _, _, fr, to, length, start, stop, nd = e
        assert 0 <= start <= stop <= 491 and stop - start + 1 <= length and fr != to
    ds = S.to_dataset(res, tree, si)
    assert ds.L == 489 and ds.n_pairs == 119316


def test_event_counts_match_poisson_mean(setup):
    tree, si, _ = setup
    par = par_for(setup, 50.0)
    N = [S.simulate(par, tree, si, S.rng_for(3, 0, 50, r)[0]).summary()['N'] for r in range(400)]
    lam = 2 * par.eta * (491 + 1 / par.p) * sum(par.blen[e] for e in tree.post_edges)
    assert abs(np.mean(N) - lam) < 4 * np.sqrt(lam / 400)
