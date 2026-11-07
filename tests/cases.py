"""Fixed evaluation points for the jsonctmctree golden values (plan L2).

Each case builds an IGCexpansion object at a fixed parameter vector and returns its
log-likelihood from the current engine. The points are the ones used by the BEAGLE
feasibility spike, so the values can be compared with its recorded outputs.

IGCexpansion prints progress to stdout; ``quiet()`` keeps that out of fixture files.
"""

import contextlib
import io
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(REPO, "test")
# Yeast alignments of the IGCTract study, outside this repository.
MAFFT = os.environ.get(
    "IGC_MAFFT_DIR",
    os.path.join(os.path.dirname(REPO), "IGCTract", "MafftAlignment"),
)


@contextlib.contextmanager
def quiet():
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def _set_x(obj, kappa, tau, omega=None, blen_seed=1):
    """Spike convention: log kappa, log omega, log tau, then branch lengths ~ U(0.02, 0.25)."""
    rng = np.random.default_rng(blen_seed)
    x = np.array(obj.x, dtype=float)
    k = len(obj.edge_to_blen)
    nproc = len(x) - k
    x[3] = np.log(kappa)
    if omega is not None:
        x[4] = np.log(omega)
    x[nproc - 1] = np.log(tau)
    x[nproc:] = np.log(rng.uniform(0.02, 0.25, size=k))
    obj.update_by_x(x)
    return x


def _recodon(model, newick, aln, paralog, save_dir):
    from IGCexpansion.CodonGeneconv import ReCodonGeneconv

    return ReCodonGeneconv(
        newick, aln, paralog, Model=model, Force=None, clock=False,
        save_path=save_dir + "/", save_name=os.path.join(save_dir, "case_save.txt"),
    )


def hky_ylr406c(save_dir):
    obj = _recodon("HKY", T + "/YeastTree.newick", T + "/YLR406C_YDL075W_test_input.fasta",
                   ["YLR406C", "YDL075W"], save_dir)
    _set_x(obj, kappa=3.0, tau=2.0)
    return obj


def mg94_ylr406c(save_dir):
    obj = _recodon("MG94", T + "/YeastTree.newick", T + "/YLR406C_YDL075W_test_input.fasta",
                   ["YLR406C", "YDL075W"], save_dir)
    _set_x(obj, kappa=3.0, tau=2.0, omega=0.3)
    return obj


def mg94_ygl062w(save_dir):
    """Longest yeast pair (1129 codons, 1013 unique patterns); the IS-IGC benchmark case."""
    pair = "YGL062W_YBR218C"
    aln = os.path.join(MAFFT, pair, pair + "_input.fasta")
    if not os.path.exists(aln):
        raise FileNotFoundError(aln)
    obj = _recodon("MG94", T + "/YeastTree.newick", aln, ["__YGL062W", "__YBR218C"], save_dir)
    _set_x(obj, kappa=3.0, tau=2.0, omega=0.3, blen_seed=3)
    return obj


def psjs_ydr418w(save_dir):
    """PS-IGC (PSJSGeneconv, HKY, one rate) on YDR418W_YEL054C, as in PSJSGeneconv's own test."""
    from IGCexpansion.PSJSGeneconv import PSJSGeneconv

    x_js = np.log([0.5, 0.5, 0.5, 4.35588244, 1.5, 5.0, 0.3, 1.0 / 30.0])
    return PSJSGeneconv(
        T + "/YDR418W_YEL054C_input.fasta", T + "/YDR418W_YEL054C_GeneToOrlg.txt", None, True, True,
        T + "/YeastTree.newick", T + "/YeastTestDupLost.txt", x_js, "HKY", "One rate", True,
        {"D1": 0}, ["kluyveri", "castellii", "bayanus", "kudriavzevii", "mikatae", "paradoxus", "cerevisiae"],
        os.path.join(save_dir, "psjs_save.txt"), os.path.join(save_dir, "psjs_log.txt"),
        force={6: 0.0}, space_list=None,
    )


# PS-IGC pairwise scenes checked by PSJSGeneconv's own test: (n, codon_site_pair).
PSJS_SCENES = [(1, (1, 2)), (1, (2, 3)), (2, (1, 3))]

# Legacy quirk: PSJSGeneconv.__init__ leaves every branch length at 1.0 and does not apply
# the rate part of self.x; the first unpack_x(self.x) changes the likelihood. Package-test
# values (psjs_ydr418w) are at the constructor state; psjs_ydr418w_at_x is at unpack_x(self.x),
# the point that fitting code (objective_and_gradient) actually evaluates.


def evaluate(name, save_dir):
    """Return {key: lnL} for one case."""
    with quiet():
        if name in ("psjs_ydr418w", "psjs_ydr418w_at_x"):
            obj = psjs_ydr418w(save_dir)
            if name == "psjs_ydr418w_at_x":
                obj.unpack_x(np.array(obj.x, float))
            out = {}
            for n, pair in PSJS_SCENES:
                v = obj._loglikelihood(n, pair)
                out["%d_%d_%d" % (n, *pair)] = float(v[0] if isinstance(v, (tuple, list)) else v)
            return out
        obj = {"hky_ylr406c": hky_ylr406c, "mg94_ylr406c": mg94_ylr406c,
               "mg94_ygl062w": mg94_ygl062w}[name](save_dir)
        ll, _ = obj._loglikelihood(store=True, edge_derivative=False)
        return {"lnL": float(ll)}


CASES = ["hky_ylr406c", "mg94_ylr406c", "psjs_ydr418w", "psjs_ydr418w_at_x", "mg94_ygl062w"]
SLOW = {"mg94_ygl062w"}
EXTERNAL = {"mg94_ygl062w"}

# Values recorded by the spike (2026-11-06; PYTHONHASHSEED=0, numpy 2.5, scipy 1.18).
SPIKE = {
    "hky_ylr406c": {"lnL": -1821.1487039278636},
    "mg94_ylr406c": {"lnL": -1328.8276129742574},
    "mg94_ygl062w": {"lnL": -20782.709884346277},
    "psjs_ydr418w": {"1_1_2": -5251.584676106582, "1_2_3": -5693.500662621693, "2_1_3": -5431.747983027612},
}
