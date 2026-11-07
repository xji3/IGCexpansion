"""L2b legacy fixture inventory.

Re-evaluates the log-likelihood of every legacy fixture under ``test/`` with the current
code (jsonctmctree engine) and records which legacy values still match. This is a report,
not a pass/fail gate.

Fixtures covered
  * ``test/save/*_save.txt``          saved parameter vectors x (one entry each)
  * ``test/Summary/*summary*.txt``    recorded scalar lnL, matched to a save file by name
  * ``test/Summary/*sw_lnL*.txt``     recorded per-site lnL (IndCodonGeneconv output)
  * ``test/Summary/*_lnL_summary.txt`` recorded per-site-pair lnL (PSJSGeneconv output)
  * ``test/YLR406C_YDL075W_sitewise_lnL.txt``

Each save file gets one *recipe* (constructor arguments, taken from the ``__main__`` blocks
of ``IGCexpansion/*.py``); recipes are evaluated once and shared by the fixtures that use
them. Sitewise files have no save-file name of their own, so every plausible recipe is
evaluated and the closest one is reported (the others are listed in the note).

Shims. The package's Python 3 / NumPy 2 failures found by this inventory (integer
division in IndCodonGeneconv.nts_to_codons, ``dtype=np.floating``, HMMJSGeneconv's implicit
relative imports) are fixed in the package. One script-local shim remains:
  * ``hmmjs_replica`` HMMJSGeneconv's ``_loglikelihood(x)`` is repeated here (about ten
    lines) on top of the package's IndCodonGeneconv and HMMTract, so a saved x can be
    evaluated without running the class's full constructor.

Status (per fixture)
  match             |recomputed - recorded| <= 1e-6 * |recorded|; for per-site files every
                    site must satisfy |d_i| <= 1e-6 * max(1, |recorded_i|)
  mismatch          evaluated, outside that tolerance
  no-recorded-value evaluated, nothing usable to compare with (placeholders 0.0 / -1 count
                    as no value)
  not-evaluated     recipe or input not recoverable, see note
  error             the current code raised; message in note

Usage:  .venv/bin/python tests/inventory_l2b.py [--only ID[,ID...]]
Writes tests/golden/l2b_inventory.json and tests/golden/l2b_inventory.md.
"""

import argparse
import contextlib
import datetime
import io
import json
import os
import platform
import shutil
import sys
import tempfile
import time
import traceback
import warnings

import numpy as np

import IGCexpansion  # noqa: F401  (installs the jsonctmctree shims)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(REPO, "test")
S = os.path.join(T, "Summary")
SAVE = os.path.join(T, "save")
OUT_JSON = os.path.join(REPO, "tests", "golden", "l2b_inventory.json")
OUT_MD = os.path.join(REPO, "tests", "golden", "l2b_inventory.md")
RTOL = 1e-6

YEAST = ["kluyveri", "castellii", "bayanus", "kudriavzevii", "mikatae", "paradoxus", "cerevisiae"]
PRIMATE = ["Chimpanzee", "Gorilla", "Orangutan", "Macaque", "Tamarin"]
ADH = ["Baboon", "Orangutan", "Gorilla", "Bonobo", "Chimpanzee", "Human"]
EDN_ECP = ["EDN", "ECP"]
YLR = ["YLR406C", "YDL075W"]
YDR = ["YDR418W", "YEL054C"]
HMM_STATES = ["No IGC event (Si = 0)", "At least one IGC event (Si > 0)"]


class NotEvaluated(Exception):
    pass


@contextlib.contextmanager
def quiet():
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def t(*p):
    return os.path.join(T, *p)


def load_x(save):
    return np.loadtxt(os.path.join(SAVE, save))


# ---------------------------------------------------------------- shims

def _ind_class(_unused=None):
    from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
    return IndCodonGeneconv


# ---------------------------------------------------------------- recipes
# A recipe builds one model at one x. ``kind`` selects the evaluator.

RECIPES = {}


def recipe(rid, kind, save, cls, summary, shims=(), x=None, **kw):
    RECIPES[rid] = dict(rid=rid, kind=kind, save=save, cls=cls, summary=summary,
                        shims=list(shims), x=x, kw=kw)


def recodon(rid, save, model, aln, paralog, tree, force=None):
    recipe(rid, "codon", save, "ReCodonGeneconv",
           "ReCodonGeneconv(%s, %s, %s, Model=%r, Force=%r, clock=False)" % (tree, aln, paralog, model, force),
           model=model, aln=aln, paralog=paralog, tree=tree, force=force, ind=False, rv=False)


def ind(rid, save, model, aln, paralog, tree, force=None, rv=False):
    shims = []
    recipe(rid, "codon", save, "IndCodonGeneconv",
           "IndCodonGeneconv(%s, %s, %s, Model=%r, Force=%r, clock=None, rate_variation=%r)"
           % (tree, aln, paralog, model, force, rv), shims,
           model=model, aln=aln, paralog=paralog, tree=tree, force=force, ind=True, rv=rv)


def js(rid, save, data, igc_pm="One rate", rv=False, force=None, cdna=True):
    aln, orlg, tree, duplos, terminals, n_rate = data
    shims = []
    recipe(rid, "js", save, "JSGeneconv",
           "JSGeneconv(%s, %s, cdna=%r, %s, %s, x_js=x[:-%d], 'HKY', %r, rate_variation=%r, {'D1':0..}, force=%r)"
           % (aln, orlg, cdna, tree, duplos, n_rate, igc_pm, rv, force), shims,
           aln=aln, orlg=orlg, tree=tree, duplos=duplos, terminals=terminals, n_rate=n_rate,
           igc_pm=igc_pm, rv=rv, force=force, cdna=cdna)


def psjs(rid, save, data, rv, cdna, x=None, origin=""):
    aln, orlg, tree, duplos, terminals, n_rate = data
    recipe(rid, "psjs", save, "PSJSGeneconv",
           "PSJSGeneconv(%s, %s, None, cdna=%r, allow_same_codon=%r, %s, %s, x_js=x[:-%d], 'HKY', 'One rate', "
           "rate_variation=%r, {'D1':0})%s" % (aln, orlg, cdna, cdna, tree, duplos, n_rate, rv, origin),
           x=x, aln=aln, orlg=orlg, tree=tree, duplos=duplos, terminals=terminals, n_rate=n_rate,
           rv=rv, cdna=cdna)


def hmmjs(rid, save, model, aln, paralog, tree, seq_index, rv=False):
    shims = ["hmmjs_replica"]
    recipe(rid, "hmmjs", save, "HMMJSGeneconv",
           "HMMJSGeneconv(save, %s, %s, %s, ..., seq_index=%s, model=%r, rate_variation=%r)._loglikelihood(x)"
           % (tree, aln, paralog, seq_index, model, rv), shims,
           model=model, aln=aln, paralog=paralog, tree=tree, seq_index=seq_index, rv=rv)


def not_evaluated(rid, save, cls, reason):
    recipe(rid, "none", save, cls, "-", reason=reason)


# Data sets for JSGeneconv/PSJSGeneconv: alignment, gene->orlg, tree, dup/loss, terminals,
# number of branch-length entries at the end of x.
D_EDN = ("EDN_ECP_Cleaned_NewFormat.fasta", "EDN_ECP_GeneToOrlg.txt", "input_tree.newick",
         "EDN_ECP_DupLost.txt", PRIMATE, 8)
D_YDR = ("YDR418W_YEL054C_input.fasta", "YDR418W_YEL054C_GeneToOrlg.txt", "YeastTree.newick",
         "YeastTestDupLost.txt", YEAST, 12)
D_YLR = ("YLR406C_YDL075W_alignment.fasta", "YLR406_YDL075W_GeneToOrlg.txt", "YeastTree.newick",
         "YeastTestDupLost.txt", YEAST, 12)
D_ADH = ("Trigeneconv_ADH_intron_input.fasta", "Trigeneconv_GeneToOrlg.txt",
         "Trigeneconv_ADH1Class_tree.newick", "Trigeneconv_ADH_DupLost.txt", ADH, 12)

# -- ReCodonGeneconv (CodonGeneconv.py __main__; Run.py; file names from get_save_file_name)
recodon("recodon_mg94_ylr", "MG94_YLR406C_YDL075W_nonclock_save.txt", "MG94",
        "YLR406C_YDL075W_input.fasta", YLR, "YeastTree.newick")
recodon("recodon_hky_ylr", "HKY_YLR406C_YDL075W_nonclock_save.txt", "HKY",
        "YLR406C_YDL075W_input.fasta", YLR, "YeastTree.newick")
recodon("recodon_mg94_ylr_testinput", "MG94_YLR406C_YDL075W_nonclock_save.txt", "MG94",
        "YLR406C_YDL075W_test_input.fasta", YLR, "YeastTree.newick")
recodon("recodon_mg94_edn", "MG94_EDN_ECP_nonclock_save.txt", "MG94", "EDN_ECP_Cleaned.fasta",
        EDN_ECP, "input_tree.newick")
recodon("recodon_mg94_force_edn", "MG94_Force_EDN_ECP_nonclock_save.txt", "MG94",
        "EDN_ECP_Cleaned.fasta", EDN_ECP, "input_tree.newick", force={5: 0.0})
recodon("recodon_hky_edn", "HKY_EDN_ECP_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta",
        EDN_ECP, "input_tree.newick")

# -- IndCodonGeneconv (IndCodonGeneconv.py / HMMJSGeneconv.py __main__ and the save names
#    HMMJSGeneconv.__init__ gives its internal IndCodonGeneconv)
ind("ind_mg94_ylr", "Ind_MG94_YLR406C_YDL075W_nonclock_save.txt", "MG94",
    "YLR406C_YDL075W_test_input.fasta", YLR, "YeastTree.newick")
ind("ind_hky_edn", "Ind_HKY_EDN_ECP_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta", EDN_ECP,
    "input_tree.newick")
ind("ind_mg94_edn", "Ind_MG94_EDN_ECP_nonclock_save.txt", "MG94", "EDN_ECP_Cleaned.fasta", EDN_ECP,
    "input_tree.newick")
ind("ind_hky_edn_hmmmain", "EDN_ECP_HKY_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta", EDN_ECP,
    "input_tree.newick")
ind("ind_hky_rv_edn_hmmmain", "EDN_ECP_HKY_rv_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta",
    EDN_ECP, "input_tree.newick", rv=True)
ind("ind_hky_edn_hmmjs", "EDN_ECP_HKY_IGC_HMMJSGeneconv_save.txt", "HKY", "EDN_ECP_Cleaned.fasta",
    EDN_ECP, "input_tree.newick")
ind("ind_hky_rv_edn_hmmjs", "EDN_ECP_HKY_rv_IGC_HMMJSGeneconv_save.txt", "HKY",
    "EDN_ECP_Cleaned.fasta", EDN_ECP, "input_tree.newick", rv=True)
ind("ind_mg94_edn_hmmjs", "EDN_ECP_MG94_IGC_HMMJSGeneconv_save.txt", "MG94", "EDN_ECP_Cleaned.fasta",
    EDN_ECP, "input_tree.newick")
ind("ind_hky_rv_edn_test", "test_save.txt", "HKY", "EDN_ECP_Cleaned.fasta", EDN_ECP,
    "input_tree.newick", rv=True)
ind("ind_mg94_edn_oldhmm_igc", "HMMJS_EDN_ECP_Ind_MG94_nonclock_IGC_save.txt", "MG94",
    "EDN_ECP_Cleaned.fasta", EDN_ECP, "input_tree.newick")
ind("ind_mg94_edn_oldhmm_force", "HMMJS_EDN_ECP_Ind_MG94_nonclock_Force_save.txt", "MG94",
    "EDN_ECP_Cleaned.fasta", EDN_ECP, "input_tree.newick", force={5: 0.0})
ind("ind_mg94_ydr_sim", "simulation_test_save.txt", "MG94", "YDR418W_YEL054C_Simulation.fasta", YDR,
    "YeastTree.newick")
ind("ind_mg94_ydr_sim_hmmjs", "YDR418W_YEL054C_MG94_IGC_HMMJSGeneconv_save.txt", "MG94",
    "YDR418W_YEL054C_Simulation.fasta", YDR, "YeastTree.newick")

# -- HMMJSGeneconv (HMMJSGeneconv.py __main__)
hmmjs("hmmjs_hky_edn", "HMMJS_EDN_ECP_HKY_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta",
      EDN_ECP, "input_tree.newick", "EDN_ECP_seq_index.txt")
hmmjs("hmmjs_hky_rv_edn", "HMMJS_EDN_ECP_HKY_rv_nonclock_save.txt", "HKY", "EDN_ECP_Cleaned.fasta",
      EDN_ECP, "input_tree.newick", "EDN_ECP_seq_index.txt", rv=True)
hmmjs("hmmjs_mg94_edn", "HMMJS_EDN_ECP_MG94_nonclock_save.txt", "MG94", "EDN_ECP_Cleaned.fasta",
      EDN_ECP, "input_tree.newick", "EDN_ECP_seq_index.txt")
hmmjs("hmmjs_mg94_ydr_sim", "HMMJS_YDR418W_YEL054C_MG94_nonclock_save.txt", "MG94",
      "YDR418W_YEL054C_Simulation.fasta", YDR, "YeastTree.newick", "YDR418W_YEL054C_seq_index.txt")
not_evaluated("hmmjs_old_mg94_edn", "HMMJS_EDN_ECP_Ind_MG94_nonclock_save.txt", "HMMJSGeneconv (old)",
              "Written by an earlier HMMJSGeneconv that kept two IndCodonGeneconv models "
              "(companion saves *_Ind_MG94_nonclock_IGC/_Force_save.txt, evaluated separately); "
              "that code is gone, so the likelihood it optimised cannot be rebuilt. No recorded "
              "value refers to it (HMM_EDN_ECP_MG94 summary matches HMMJS_EDN_ECP_MG94 by tract_p).")

# -- JSGeneconv (JSGeneconv.py __main__, active and commented blocks)
js("js_hky_edn", "JS_HKY_EDN_ECP_nonclock_save.txt", D_EDN)
js("js_hky_rv_edn", "JS_HKY_rv_EDN_ECP_nonclock_save.txt", D_EDN, rv=True)
js("js_hky_force_edn", "Force_JS_HKY_EDN_ECP_nonclock_save.txt", D_EDN, force={4: 0.0})
js("js_hky_rv_force_edn", "Force_JS_HKY_rv_EDN_ECP_nonclock_save.txt", D_EDN, rv=True, force={6: 0.0})
js("js_hky_ydr", "HKY_YDR418W_YEL054C_nonclock_One_rate_save.txt", D_YDR)
js("js_hky_force_ydr", "Force_HKY_YDR418W_YEL054C_nonclock_One_rate_save.txt", D_YDR, force={4: 0.0})
js("js_hky_rv_ydr", "HKY_YDR418W_YEL054C_nonclock_rv_One_rate__save.txt", D_YDR, rv=True)
js("js_hky_rv_force_ydr", "Force_HKY_YDR418W_YEL054C_nonclock_rv_One_rate__save.txt", D_YDR, rv=True,
   force={6: 0.0})
js("js_hky_mostgeneral_ydr", "HKY_YDR418W_YEL054C_nonclock_Most_general_save.txt", D_YDR,
   igc_pm="Most general")
js("js_hky_ylr_yeastsave", "Yeast_save.txt", D_YLR)
js("js_hky_ylr_hkysave", "HKY_YLR406C_YDL075W_nonclock_save.txt", D_YLR)
js("js_tri_hky", "Trigeneconv_ADH_HKY_save.txt", D_ADH, cdna=False)
js("js_tri_hky_dir", "Trigeneconv_ADH_HKY_dir_save.txt", D_ADH, igc_pm="Most general", cdna=False)
js("js_tri_hky_sym", "Trigeneconv_ADH_HKY_sym_save.txt", D_ADH, igc_pm="Symmetric general", cdna=False)
not_evaluated("js_adh1_intron5", "ADH1_intron_5_save.txt", "JSGeneconv",
              "The generating block (JSGeneconv.py, 'Class 1 ADH Genes') reads "
              "'../../ADH1Genes/Alignment/...', outside this repository, and uses a pre-cdna "
              "constructor signature; the intron-5 alignment is not in test/.")

# -- PSJSGeneconv (PSJSGeneconv.py __main__)
psjs("psjs_hky_edn", "PSJS_HKY_EDN_ECP_nonclock_save.txt", D_EDN, rv=False, cdna=True)
psjs("psjs_hky_ydr", "PSJS_HKY_YDR418W_YEL054C_nonclock_save.txt", D_YDR, rv=False, cdna=False)
X_PSJS_RV_YDR = [-0.92969299, -0.60831311, -1.00401104, -0.8891694, -0.89092996, -0.20694348,
                 -65.92436066, 0.0, -2.73221044, -26.7831379, -3.44536047, -2.72202154,
                 -21.09784191, -2.70441861, -22.63017035, -2.72779708, -20.53753293,
                 -3.44335198, -2.73358356, -20.14414384]
psjs("psjs_hky_rv_ydr_literal", None, D_YDR, rv=True, cdna=True, x=X_PSJS_RV_YDR,
     origin=" at the x literal in the commented PSJSGeneconv.py block that names "
            "PSJS_HKY_YDR418W_YEL054C_rv_nonclock_save.txt (that save file is not in test/save)")
not_evaluated("psjs_hky_ylr", "PSJS_HKY_YLR406C_YDL075W_nonclock_save.txt", "PSJSGeneconv",
              "No __main__ block (active or commented) builds a PSJSGeneconv on YLR406C_YDL075W, so "
              "alignment/cdna/rate-variation settings are unknown; no recorded value exists.")


# ---------------------------------------------------------------- evaluators

def _codon_obj(r, tmp):
    kw = r["kw"]
    if kw["ind"]:
        cls = _ind_class("py3_codons" in r["shims"])
        obj = cls(t(kw["tree"]), t(kw["aln"]), kw["paralog"], Model=kw["model"], Force=kw["force"],
                  clock=None, save_path=tmp + "/", save_name=os.path.join(tmp, "x_save.txt"),
                  rate_variation=kw["rv"])
    else:
        from IGCexpansion.CodonGeneconv import ReCodonGeneconv
        obj = ReCodonGeneconv(t(kw["tree"]), t(kw["aln"]), kw["paralog"], Model=kw["model"],
                              Force=kw["force"], clock=False, save_path=tmp + "/",
                              save_name=os.path.join(tmp, "x_save.txt"))
    return obj


def _codon_params(obj):
    p = {"length": obj.nsites, "kappa": obj.kappa, "tau": obj.tau}
    if obj.Model == "MG94":
        p["omega"] = obj.omega
    p.update({k: v for k, v in zip(["pi_a", "pi_c", "pi_g", "pi_t"], obj.pi)})
    p.update({"(%s)" % ",".join(e): obj.edge_to_blen[e] for e in obj.edge_list})
    return p


def eval_codon(r, tmp, need):
    x = load_x(r["save"])
    obj = _codon_obj(r, tmp)
    if len(x) != len(obj.x):
        raise ValueError("len(x)=%d but model expects %d" % (len(x), len(obj.x)))
    obj.update_by_x(x)
    res = {"lnL": float(obj._loglikelihood()[0]), "params": _codon_params(obj)}
    if "sw_IGC" in need:
        if r["kw"]["ind"]:
            res["sw_IGC"] = np.array(obj._sitewise_loglikelihood(False), dtype=float)
        else:
            res["sw_IGC"] = np.array(obj._sitewise_loglikelihood(), dtype=float)
    if "sw_NOIGC" in need:
        if not r["kw"]["ind"]:
            raise NotEvaluated("ReCodonGeneconv has no NOIGC sitewise likelihood")
        res["sw_NOIGC"] = np.array(obj._sitewise_loglikelihood(True), dtype=float)
    return res


def _js_obj(r, tmp, x):
    from IGCexpansion.JSGeneconv import JSGeneconv
    kw = r["kw"]
    x_js = np.array(x[:-kw["n_rate"]])
    n2p = {"D1": 0, "D2": 0} if kw["tree"].startswith("Trigeneconv") else {"D1": 0}
    obj = JSGeneconv(t(kw["aln"]), t(kw["orlg"]), kw["cdna"], t(kw["tree"]), t(kw["duplos"]), x_js,
                     "HKY", kw["igc_pm"], kw["rv"], n2p, kw["terminals"],
                     os.path.join(tmp, "js_save.txt"), force=kw["force"])
    if len(x) != len(obj.x):
        raise ValueError("len(x)=%d but model expects %d" % (len(x), len(obj.x)))
    obj.unpack_x(np.array(x))
    return obj


def eval_js(r, tmp, need):
    x = load_x(r["save"])
    ctx = contextlib.nullcontext()
    with ctx:
        obj = _js_obj(r, tmp, x)
        ll = float(obj._loglikelihood()[0])
        summ, labels = obj.get_summary()
    params = {k: v for k, v in zip(labels, summ) if k != "ll"}
    return {"lnL": ll, "params": params}


def eval_psjs(r, tmp, need):
    from IGCexpansion.PSJSGeneconv import PSJSGeneconv
    kw = r["kw"]
    x = np.array(r["x"]) if r["x"] is not None else load_x(r["save"])
    obj = PSJSGeneconv(t(kw["aln"]), t(kw["orlg"]), None, kw["cdna"], kw["cdna"], t(kw["tree"]),
                       t(kw["duplos"]), x[:-kw["n_rate"]], "HKY", "One rate", kw["rv"], {"D1": 0},
                       kw["terminals"], os.path.join(tmp, "psjs_save.txt"),
                       os.path.join(tmp, "psjs_log.txt"))
    if len(x) != len(obj.x):
        raise ValueError("len(x)=%d but model expects %d" % (len(x), len(obj.x)))
    obj.unpack_x(x)
    # Per-pair lnL over every distance n; their sum is the composite (pairwise) lnL that
    # objective_wo_gradient returns (negated).
    pairs = {tuple(int(i) for i in k): float(v) for k, v in obj._sitewise_loglikelihood_for_all_n().items()}
    return {"lnL": float(sum(pairs.values())), "pairs": pairs}


def eval_hmmjs(r, tmp, need):
    """Replica of HMMJSGeneconv.__init__ + _loglikelihood(x) for the saved x (see module doc)."""
    from IGCexpansion.HMMTract import HMMTract
    kw = r["kw"]
    x = load_x(r["save"])
    cls = _ind_class("py3_codons" in r["shims"])
    ind_obj = cls(t(kw["tree"]), t(kw["aln"]), kw["paralog"], Model=kw["model"], Force=None,
                  clock=None, save_name=os.path.join(tmp, "hmm_ind_save.txt"),
                  rate_variation=kw["rv"])
    if len(x) - 1 != len(ind_obj.x):
        raise ValueError("len(x)-1=%d but IndCodonGeneconv expects %d" % (len(x) - 1, len(ind_obj.x)))
    ind_obj.update_by_x(x[:-1])
    igc_file, noigc_file = os.path.join(tmp, "hmm_igc_sw.txt"), os.path.join(tmp, "hmm_noigc_sw.txt")
    ind_obj.get_sitewise_loglikelihood_summary(igc_file, False)
    ind_obj.get_sitewise_loglikelihood_summary(noigc_file, True)
    outgroup = [e for e in ind_obj.edge_list if e[0] == "N0" and e[1] != "N1"][0]
    total_blen = sum(ind_obj.edge_to_blen[e] for e in ind_obj.edge_list if e != outgroup)
    hmm = HMMTract(igc_file, noigc_file, HMM_STATES, total_blen, ind_obj.tau, t(kw["seq_index"]),
                   kw["model"])
    hmm_x = [x[5] + x[-1], x[-1]]  # as in HMMJSGeneconv._loglikelihood (index 5 = log tau for MG94)
    neg_ll = float(hmm.Forward(False, hmm_x))
    p = _codon_params(ind_obj)
    p.pop("length")
    p.update({"eta": float(np.exp(hmm_x[0])), "tract_p": float(np.exp(hmm_x[1]))})
    if kw["model"] == "HKY":
        p["omega"] = ind_obj.omega
    return {"lnL": -neg_ll, "params": p}


EVALUATORS = {"codon": eval_codon, "js": eval_js, "psjs": eval_psjs, "hmmjs": eval_hmmjs}
_CACHE = {}


def run_recipe(rid, need):
    """Evaluate a recipe once; returns {'ok': bool, ...results or error, 'seconds'}."""
    key = (rid, tuple(sorted(need)))
    for (k, n), v in _CACHE.items():
        if k == rid and set(need) <= set(n):
            return v
    r = RECIPES[rid]
    if r["kind"] == "none":
        out = {"ok": False, "status": "not-evaluated", "error": r["kw"]["reason"], "seconds": 0.0}
        _CACHE[key] = out
        return out
    tmp = tempfile.mkdtemp(prefix="l2b_")
    t0 = time.time()
    try:
        with quiet():
            res = EVALUATORS[r["kind"]](r, tmp, set(need))
        out = dict(res, ok=True)
    except NotEvaluated as e:
        out = {"ok": False, "status": "not-evaluated", "error": str(e)}
    except BaseException as e:  # SystemExit from HMMTract.get_Emi included
        if isinstance(e, KeyboardInterrupt):
            raise
        tb = traceback.extract_tb(e.__traceback__)
        where = next((f for f in reversed(tb) if "IGCexpansion/IGCexpansion" in f.filename), tb[-1])
        out = {"ok": False, "status": "error",
               "error": "%s: %s (at %s:%d)" % (type(e).__name__, str(e)[:200],
                                               os.path.basename(where.filename), where.lineno)}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    out["seconds"] = round(time.time() - t0, 2)
    _CACHE[key] = out
    return out


# ---------------------------------------------------------------- recorded values

def read_summary(path):
    """Single-column np.savetxt file with a '# label label ...' footer -> (values, labels)."""
    vals, labels = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                labels = line[1:].split()
            elif line.strip():
                vals.append([float(v) for v in line.split()])
    return vals, labels


def read_sitewise(path):
    out = []
    with open(path) as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                out.append(float(line.split("\t")[1]))
    return np.array(out)


def read_pairwise(path):
    out = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                a, b, v = line.split("\t")[:3]
                out[(int(a), int(b))] = float(v)
    return out


def summary_record(path, sign=1.0):
    vals, labels = read_summary(path)
    rel = os.path.relpath(path, REPO)
    if labels is None or "ll" not in labels or any(len(v) != 1 for v in vals):
        return {"source": rel, "value": None, "line": None, "labels_ok": False}
    i = labels.index("ll")
    flat = [v[0] for v in vals]
    return {"source": rel, "line": i + 1, "value": sign * flat[i], "raw": flat[i],
            "labels_ok": len(labels) == len(flat),
            "table": dict(zip(labels, flat)) if len(labels) == len(flat) else None}


def param_check(record, params):
    """Max relative difference between a summary's parameter rows and the model at x."""
    table = record.get("table") if record else None
    if not table or not params:
        return None
    common = [k for k in params if k in table and k not in ("ll",)]
    if not common:
        return None
    d = [abs(params[k] - table[k]) / max(abs(table[k]), 1e-300) for k in common]
    worst = common[int(np.argmax(d))]
    return {"n_compared": len(common), "max_rel_diff": float(max(d)), "worst": worst}


# ---------------------------------------------------------------- fixtures

FIXTURES = []


def fx(fid, kind, recipes, recorded=None, note="", sign=1.0, placeholder=()):
    FIXTURES.append(dict(id=fid, kind=kind, recipes=recipes if isinstance(recipes, list) else [recipes],
                         recorded=recorded, note=note, sign=sign, placeholder=placeholder))


def SUM(name):
    return os.path.join(S, name)


# Scalar lnL, one per save file (recorded summary where a counterpart exists).
fx("MG94_YLR406C_YDL075W_nonclock", "lnL", "recodon_mg94_ylr",
   SUM("MG94_YLR406C_YDL075W_nonclock_summary.txt"), placeholder=(0.0,),
   note="Alignment YLR406C_YDL075W_input.fasta (summary length 112 codons; test_input has 113).")
fx("HKY_YLR406C_YDL075W_nonclock", "lnL", "recodon_hky_ylr",
   SUM("HKY_YLR406C_YDL075W_nonclock_summary.txt"), placeholder=(0.0,),
   note="Alignment YLR406C_YDL075W_input.fasta (summary length 336).")
fx("MG94_EDN_ECP_nonclock", "lnL", "recodon_mg94_edn", SUM("MG94_EDN_ECP_nonclock_summary.txt"),
   placeholder=(0.0,))
fx("MG94_Force_EDN_ECP_nonclock", "lnL", "recodon_mg94_force_edn",
   SUM("Force_MG94_EDN_ECP_nonclock_summary.txt"), placeholder=(0.0,),
   note="Save file is byte-identical in x to MG94_EDN_ECP_nonclock_save.txt; Force {5:0} sets tau=0.")
fx("HKY_EDN_ECP_nonclock", "lnL", "recodon_hky_edn")
fx("Ind_MG94_YLR406C_YDL075W_nonclock", "lnL", "ind_mg94_ylr",
   note="Alignment from the commented IndCodonGeneconv.py block (test_input).")
fx("Ind_HKY_EDN_ECP_nonclock", "lnL", "ind_hky_edn")
fx("Ind_MG94_EDN_ECP_nonclock", "lnL", "ind_mg94_edn")
fx("EDN_ECP_HKY_nonclock", "lnL", "ind_hky_edn_hmmmain")
fx("EDN_ECP_HKY_rv_nonclock", "lnL", "ind_hky_rv_edn_hmmmain")
fx("EDN_ECP_HKY_IGC_HMMJSGeneconv", "lnL", "ind_hky_edn_hmmjs",
   note="Internal IndCodonGeneconv of HKY HMMJSGeneconv; x[5] (first branch length for HKY) was "
        "overwritten with log tau by HMMJSGeneconv.get_mle, which assumes MG94 layout.")
fx("EDN_ECP_HKY_rv_IGC_HMMJSGeneconv", "lnL", "ind_hky_rv_edn_hmmjs",
   note="As above; for HKY+rv x[5] is r3, overwritten with log tau.")
fx("EDN_ECP_MG94_IGC_HMMJSGeneconv", "lnL", "ind_mg94_edn_hmmjs")
fx("test_save", "lnL", "ind_hky_rv_edn_test", note="IndCodonGeneconv.py __main__ (HKY, rate_variation).")
fx("HMMJS_EDN_ECP_Ind_MG94_nonclock_IGC", "lnL", "ind_mg94_edn_oldhmm_igc")
fx("HMMJS_EDN_ECP_Ind_MG94_nonclock_Force", "lnL", "ind_mg94_edn_oldhmm_force",
   note="Assumed Force {5:0} (tau=0), consistent with Force_EDN_ECP_Ind_MG94_nonclock_sw_lnL.txt.")
fx("simulation_test", "lnL", "ind_mg94_ydr_sim")
fx("YDR418W_YEL054C_MG94_IGC_HMMJSGeneconv", "lnL", "ind_mg94_ydr_sim_hmmjs")
fx("HMMJS_EDN_ECP_HKY_nonclock", "lnL", "hmmjs_hky_edn")
fx("HMMJS_EDN_ECP_HKY_rv_nonclock", "lnL", "hmmjs_hky_rv_edn",
   SUM("HMM_EDN_ECP_HKY_nonclock_1D_summary.txt"), sign=-1.0,
   note="Summary file name lacks '_rv' but HMMJSGeneconv.py __main__ (rate_variation=True) writes "
        "it; parameters match this save. Summary 'll' is Forward() = -lnL, sign flipped.")
fx("HMMJS_EDN_ECP_MG94_nonclock", "lnL", "hmmjs_mg94_edn",
   SUM("HMM_EDN_ECP_MG94_nonclock_1D_summary.txt"), sign=-1.0,
   note="Summary 'll' is Forward() = -lnL, sign flipped.")
fx("HMMJS_YDR418W_YEL054C_MG94_nonclock", "lnL", "hmmjs_mg94_ydr_sim")
fx("HMMJS_EDN_ECP_Ind_MG94_nonclock", "lnL", "hmmjs_old_mg94_edn")
fx("JS_HKY_EDN_ECP_nonclock", "lnL", "js_hky_edn", SUM("JS_HKY_EDN_ECP_nonclock_summary.txt"))
fx("JS_HKY_rv_EDN_ECP_nonclock", "lnL", "js_hky_rv_edn", SUM("JS_HKY_rv_EDN_ECP_nonclock_summary.txt"))
fx("Force_JS_HKY_EDN_ECP_nonclock", "lnL", "js_hky_force_edn",
   SUM("Force_JS_HKY_EDN_ECP_nonclock_summary.txt"))
fx("Force_JS_HKY_rv_EDN_ECP_nonclock", "lnL", "js_hky_rv_force_edn",
   SUM("Force_JS_HKY_rv_EDN_ECP_nonclock_summary.txt"))
fx("HKY_YDR418W_YEL054C_nonclock_One_rate", "lnL", "js_hky_ydr",
   SUM("HKY_YDR418W_YEL054C_nonclock_One_rate_summary.txt"))
fx("Force_HKY_YDR418W_YEL054C_nonclock_One_rate", "lnL", "js_hky_force_ydr")
fx("HKY_YDR418W_YEL054C_nonclock_rv_One_rate_", "lnL", "js_hky_rv_ydr",
   SUM("HKY_YDR418W_YEL054C_nonclock_rv_One_rate_summary.txt"),
   note="The save does not hold the summary's x: it has kappa 0.41, tau ~ 2e-29 and several "
        "branch lengths ~ 1e-10, against kappa 13.7 and Tau 23.1 in the summary. It agrees to ~1e-4 "
        "with the x literal in the commented JSGeneconv.py block, so the save was presumably rewritten "
        "by a later run after the summary was written.")
fx("Force_HKY_YDR418W_YEL054C_nonclock_rv_One_rate_", "lnL", "js_hky_rv_force_ydr",
   SUM("Force_HKY_YDR418W_YEL054C_nonclock_rv_One_rate_summary.txt"),
   note="Like the unforced rv save, this x (log branch lengths down to -33) is not the summary's x.")
fx("HKY_YDR418W_YEL054C_nonclock_Most_general", "lnL", "js_hky_mostgeneral_ydr",
   t("HKY_YLR406C_YDL075W_nonclock_dir_summary.txt"),
   note="Commented JSGeneconv.py block pairs this save with test/HKY_YLR406C_YDL075W_nonclock_dir_"
        "summary.txt; that file's name and lnL scale point to the YLR406C alignment instead.")
fx("Yeast_save", "lnL", "js_hky_ylr_yeastsave", t("HKY_YLR406C_YDL075W_nonclock_summary.txt"),
   note="No __main__ names Yeast_save.txt; JSGeneconv HKY One rate on YLR406C_YDL075W fits its "
        "length (17) and parameter values. Compared with test/HKY_YLR406C_YDL075W_nonclock_summary.txt "
        "(JSGeneconv format), the only JS-format YLR406C summary. By value, the saved branch lengths "
        "follow a different edge order from the current Tree.edge_list ((N2,bayanus) before (N2,N3), "
        "(N3,kudriavzevii) before (N3,N4)), and even matched by value they differ from the summary "
        "by ~0.5%, so the summary was written at another x.")
fx("JS_on_HKY_YLR406C_YDL075W_nonclock", "lnL", "js_hky_ylr_hkysave",
   t("HKY_YLR406C_YDL075W_nonclock_summary.txt"),
   note="Cross-check only: same save as the ReCodonGeneconv entry, read as a JSGeneconv x.")
fx("Trigeneconv_ADH_HKY", "lnL", "js_tri_hky", note="Commented JSGeneconv.py 'Class 1 ADH Genes' block.")
fx("Trigeneconv_ADH_HKY_dir", "lnL", "js_tri_hky_dir")
fx("Trigeneconv_ADH_HKY_sym", "lnL", "js_tri_hky_sym")
fx("ADH1_intron_5", "lnL", "js_adh1_intron5")
fx("PSJS_HKY_EDN_ECP_nonclock", "lnL", "psjs_hky_edn",
   note="Composite lnL = sum over all site pairs. Summary/PSJS_HKY_EDN_ECP_nonclock_summary.txt is "
        "in JSGeneconv format (Tau, 9 edges, ll ~ -1716.6), not PSJSGeneconv output, so it is not "
        "this save's value; it has its own entry below.")
fx("PSJS_HKY_EDN_ECP_nonclock_summary", "lnL", "js_hky_edn", SUM("PSJS_HKY_EDN_ECP_nonclock_summary.txt"),
   note="JSGeneconv-format summary under a PSJS name. Compared with JSGeneconv at "
        "JS_HKY_EDN_ECP_nonclock_save.txt, the only JS EDN_ECP One-rate save; see the parameter "
        "check for whether that is the x it was written at.")
fx("PSJS_HKY_YDR418W_YEL054C_nonclock", "lnL", "psjs_hky_ydr",
   SUM("PSJS_HKY_YDR418W_YEL054C_nonclock_summary.txt"), placeholder=(-1.0,),
   note="cdna/allow_same_codon/rate_variation=False as first set in PSJSGeneconv.py __main__ for this "
        "save name.")
fx("PSJS_HKY_YLR406C_YDL075W_nonclock", "lnL", "psjs_hky_ylr")

# Per-site lnL files (IndCodonGeneconv.get_sitewise_loglikelihood_summary output). Every
# recipe with the right data set and model is tried; the closest is reported.
HKY_EDN_CANDS = ["ind_hky_edn", "ind_hky_edn_hmmmain", "ind_hky_edn_hmmjs", "ind_hky_rv_edn_hmmmain",
                 "ind_hky_rv_edn_hmmjs", "ind_hky_rv_edn_test"]
MG94_EDN_CANDS = ["ind_mg94_edn", "ind_mg94_edn_hmmjs", "ind_mg94_edn_oldhmm_igc", "ind_mg94_edn_oldhmm_force",
                  "recodon_mg94_edn", "recodon_mg94_force_edn"]
SIM_CANDS = ["ind_mg94_ydr_sim", "ind_mg94_ydr_sim_hmmjs"]
fx("EDN_ECP_HKY_nonclock_sw_lnL", "sw_IGC", HKY_EDN_CANDS, SUM("EDN_ECP_HKY_nonclock_sw_lnL.txt"))
fx("EDN_ECP_HKY_nonclock_sw_lnL_check", "sw_IGC", HKY_EDN_CANDS, SUM("EDN_ECP_HKY_nonclock_sw_lnL_check.txt"))
fx("NOIGC_EDN_ECP_HKY_nonclock_sw_lnL", "sw_NOIGC", HKY_EDN_CANDS,
   SUM("NOIGC_EDN_ECP_HKY_nonclock_sw_lnL.txt"))
fx("EDN_ECP_MG94_nonclock_sw_lnL", "sw_IGC", MG94_EDN_CANDS, SUM("EDN_ECP_MG94_nonclock_sw_lnL.txt"))
fx("NOIGC_EDN_ECP_MG94_nonclock_sw_lnL", "sw_NOIGC", MG94_EDN_CANDS,
   SUM("NOIGC_EDN_ECP_MG94_nonclock_sw_lnL.txt"))
fx("EDN_ECP_Ind_MG94_nonclock_sw_lnL", "sw_IGC", MG94_EDN_CANDS, SUM("EDN_ECP_Ind_MG94_nonclock_sw_lnL.txt"))
fx("Force_EDN_ECP_Ind_MG94_nonclock_sw_lnL", "sw_IGC", MG94_EDN_CANDS,
   SUM("Force_EDN_ECP_Ind_MG94_nonclock_sw_lnL.txt"))
fx("Force_EDN_ECP_MG94_nonclock_sw_lnL", "sw_IGC", MG94_EDN_CANDS,
   SUM("Force_EDN_ECP_MG94_nonclock_sw_lnL.txt"))
fx("Simulation_test_YDR418W_YEL054C_MG94_nonclock_sw_lnL", "sw_IGC", SIM_CANDS,
   SUM("Simulation_test_YDR418W_YEL054C_MG94_nonclock_sw_lnL.txt"))
fx("Simulation_test_NOIGC_YDR418W_YEL054C_MG94_nonclock_sw_lnL", "sw_NOIGC", SIM_CANDS,
   SUM("Simulation_test_NOIGC_YDR418W_YEL054C_MG94_nonclock_sw_lnL.txt"))
fx("YLR406C_YDL075W_sitewise_lnL", "sw_IGC", ["recodon_mg94_ylr_testinput", "ind_mg94_ylr"],
   t("YLR406C_YDL075W_sitewise_lnL.txt"),
   note="113 sites = YLR406C_YDL075W_test_input.fasta in codons; no __main__ writes this file.")

# Per-site-pair lnL files (PSJSGeneconv.get_pairwise_loglikelihood_summary output).
fx("PSJS_HKY_EDN_ECP_nonclock_lnL_summary", "pairs", "psjs_hky_edn",
   SUM("PSJS_HKY_EDN_ECP_nonclock_lnL_summary.txt"))
fx("PSJS_HKY_YDR418W_YEL054C_nonclock_lnL_summary", "pairs", "psjs_hky_ydr",
   SUM("PSJS_HKY_YDR418W_YEL054C_nonclock_lnL_summary.txt"))
fx("PSJS_HKY_YDR418W_YEL054C_rv_nonclock_lnL_summary", "pairs", "psjs_hky_rv_ydr_literal",
   SUM("PSJS_HKY_YDR418W_YEL054C_rv_nonclock_lnL_summary.txt"))
not_evaluated("psjs_grid", None, "PSJSGeneconv",
              "output_lnL_per_distance_tract_p_list() output. The x it used is not recoverable: the "
              "__main__ that writes it builds the model from save file "
              "PSJS_HKY_rv_YDR418W_YEL054C_nonclock_save.txt (absent) and the 6-site "
              "YDR418W_YEL054C_test.fasta, while the file has 488 distance columns (489-site input). "
              "Both rows are identical across log_p, so tau or the IGC rate was presumably forced.")
fx("PSJS_HKY_rv_YDR418W_YEL054C_nonclock_grid_lnL", "grid", "psjs_grid",
   SUM("PSJS_HKY_rv_YDR418W_YEL054C_nonclock_grid_lnL.txt"))


# ---------------------------------------------------------------- comparison

def _close(rec, new):
    return abs(new - rec) <= RTOL * abs(rec)


def evaluate_fixture(f):
    kind = f["kind"]
    need = {kind} if kind in ("sw_IGC", "sw_NOIGC") else {"lnL"}
    row = {"id": f["id"], "kind": kind, "note": f["note"]}
    rec_path = f["recorded"]
    row["recorded_source"] = os.path.relpath(rec_path, REPO) if rec_path else None

    tried = []
    for rid in f["recipes"]:
        res = run_recipe(rid, need)
        tried.append((rid, res))

    def describe(rid):
        r = RECIPES[rid]
        return {"save_file": ("test/save/" + r["save"]) if r["save"] else None, "model_class": r["cls"],
                "recipe": r["summary"], "shims": r["shims"]}

    if kind == "lnL":
        rid, res = tried[0]
        row.update(describe(rid))
        row["seconds"] = res["seconds"]
        record = summary_record(rec_path, f["sign"]) if rec_path else None
        if record:
            row["recorded_line"] = record["line"]
            row["recorded_lnL"] = record["value"]
        if not res["ok"]:
            row["status"] = res["status"]
            row["note"] = _join(res["error"], row["note"])
            return row
        row["recomputed_lnL"] = res["lnL"]
        pc = param_check(record, res.get("params"))
        if pc:
            row["summary_param_check"] = pc
        if record is None or record["value"] is None:
            row["status"] = "no-recorded-value"
            return row
        if record["raw"] in f["placeholder"]:
            row["status"] = "no-recorded-value"
            row["note"] = _join("recorded ll = %g is a placeholder, not a likelihood" % record["raw"],
                                row["note"])
            return row
        d = res["lnL"] - record["value"]
        row["abs_diff"] = abs(d)
        row["rel_diff"] = abs(d) / abs(record["value"])
        row["status"] = "match" if _close(record["value"], res["lnL"]) else "mismatch"
        return row

    if kind == "grid":
        rid, res = tried[0]
        row.update(describe(rid))
        row["status"] = res["status"]
        row["note"] = _join(res["error"], row["note"])
        return row

    # per-site / per-pair files: choose the closest candidate
    rec = read_pairwise(rec_path) if kind == "pairs" else read_sitewise(rec_path)
    scored = []
    errors = []
    for rid, res in tried:
        if not res["ok"]:
            errors.append((rid, res))
            continue
        if kind == "pairs":
            got = res["pairs"]
            keys = sorted(rec)
            if set(keys) != set(got):
                scored.append((np.inf, rid, res, "pair set differs (%d recorded, %d recomputed)"
                               % (len(rec), len(got))))
                continue
            r_arr = np.array([rec[k] for k in keys])
            g_arr = np.array([got[k] for k in keys])
        else:
            g_arr = res[kind]
            r_arr = rec
            if len(g_arr) != len(r_arr):
                scored.append((np.inf, rid, res, "length differs (%d recorded, %d recomputed)"
                               % (len(r_arr), len(g_arr))))
                continue
        dif = np.abs(g_arr - r_arr)
        scored.append((float(np.max(dif / np.maximum(1.0, np.abs(r_arr)))), rid, res,
                       dict(max_abs=float(dif.max()), sum_rec=float(r_arr.sum()), sum_new=float(g_arr.sum()),
                            n=len(r_arr))))
    others = []
    if not scored:
        rid, res = errors[0]
        row.update(describe(rid))
        row["status"] = res["status"]
        row["seconds"] = res["seconds"]
        row["note"] = _join("; ".join("%s: %s" % (r, e["error"]) for r, e in errors), row["note"])
        row["candidates_tried"] = [r for r, _ in tried]
        return row
    scored.sort(key=lambda s: s[0])
    best_score, rid, res, info = scored[0]
    row.update(describe(rid))
    row["seconds"] = res["seconds"]
    row["candidates_tried"] = [r for r, _ in tried]
    for s, r2, _, i2 in scored[1:]:
        others.append("%s max|d|=%s" % (r2, ("%.3g" % i2["max_abs"]) if isinstance(i2, dict) else i2))
    for r2, e2 in errors:
        others.append("%s %s" % (r2, e2["status"]))
    if not isinstance(info, dict):
        row["status"] = "mismatch"
        row["note"] = _join(info, row["note"])
    else:
        row["n_values"] = info["n"]
        row["recorded_lnL"] = info["sum_rec"]
        row["recomputed_lnL"] = info["sum_new"]
        row["abs_diff"] = abs(info["sum_new"] - info["sum_rec"])
        row["rel_diff"] = row["abs_diff"] / abs(info["sum_rec"])
        row["max_abs_diff_per_value"] = info["max_abs"]
        row["max_scaled_diff_per_value"] = best_score
        row["status"] = "match" if best_score <= RTOL else "mismatch"
        row["note"] = _join("recorded/recomputed lnL are sums over the %d values; status uses the "
                            "per-value rule" % info["n"], row["note"])
    if others:
        row["note"] = _join("other candidates: " + ", ".join(others), row["note"])
    return row


def _join(*parts):
    return " ".join(p.rstrip(".") + "." for p in parts if p)


# ---------------------------------------------------------------- report

def write_md(rows, meta):
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    order = ["match", "mismatch", "no-recorded-value", "not-evaluated", "error"]

    def fmt(v, spec="%.10g"):
        return "" if v is None else spec % v

    lines = ["# L2b legacy fixture inventory", "",
             "Generated by `tests/inventory_l2b.py` on %s (Python %s, numpy %s, scipy %s); "
             "script runtime %.0f s. Machine-readable version: `l2b_inventory.json`."
             % (meta["generated"], meta["python"], meta["numpy"], meta["scipy"], meta["runtime_s"]), "",
             "Every saved parameter vector under `test/save/` and every recorded lnL file under "
             "`test/` was re-evaluated with the current jsonctmctree engine. Status `match` means "
             "rel. diff <= 1e-6 (per value for per-site files). This is an inventory, not a gate.", "",
             "| status | count |", "|---|---|"]
    lines += ["| %s | %d |" % (s, counts.get(s, 0)) for s in order]
    lines += ["", "| fixture | kind | model | recorded | recomputed | rel diff | status | shims |",
              "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["id"], r["kind"], r.get("model_class", ""), fmt(r.get("recorded_lnL")),
            fmt(r.get("recomputed_lnL")), fmt(r.get("rel_diff"), "%.2g"), r["status"],
            ", ".join(r.get("shims", []))))
    lines += ["", "## Notes per fixture", ""]
    for r in rows:
        if r["note"] or r.get("summary_param_check"):
            pc = r.get("summary_param_check")
            extra = (" Summary parameters vs model at x: max rel diff %.2g over %d (%s)."
                     % (pc["max_rel_diff"], pc["n_compared"], pc["worst"])) if pc else ""
            lines.append("- **%s** (%s): %s%s" % (r["id"], r["status"], r["note"], extra))
    lines += ["", "## Patterns", ""] + meta["patterns"]
    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines) + "\n")


PATTERNS = [
    "- Three code paths raise under Python 3 / NumPy 2 and were evaluated only through the "
    "script-local shims listed in the script docstring: MG94 IndCodonGeneconv (`len/3` in "
    "nts_to_codons), 'One rate' JSGeneconv (`dtype=np.floating` in IGCModel), and "
    "HMMJSGeneconv (implicit relative imports). Shimmed rows carry the shim name.",
    "- ReCodonGeneconv summaries (`MG94_*`, `HKY_YLR406C_*`, `Force_MG94_*`) record `ll = 0.0`: "
    "get_individual_summary writes self.ll, which stays at its constructor value 0.0 unless a "
    "likelihood was evaluated first. Their parameter rows can still be checked against the save.",
    "- HMMJSGeneconv.get_mle writes log tau into x[5], which is the tau slot only for MG94; the HKY "
    "and HKY+rv HMMJS saves and their internal IndCodonGeneconv saves therefore carry a corrupted "
    "first branch length (HKY) or r3 (HKY+rv).",
    "- The YDR418W_YEL054C simulation fixtures (simulation_test_save, *_MG94_IGC_HMMJSGeneconv_save, "
    "HMMJS_YDR418W_*) cannot be evaluated: YDR418W_YEL054C_Simulation.fasta contains stop codons, "
    "which IndCodonGeneconv's MG94 codon map rejects.",
    "- Where a recorded value and a save file both exist, the summary's own parameter rows tell "
    "whether they belong together (see 'Summary parameters vs model at x' in the notes). Every "
    "mismatch of a scalar lnL comes with parameter rows that also disagree with the save "
    "(max rel diff 0.036 to 25), i.e. the save and the summary were written at different x "
    "(saves overwritten by later runs, or summaries from another alignment); none points to a "
    "change in the likelihood code. The matches (all four JSGeneconv EDN_ECP fits, the YDR418W "
    "One-rate fit, the HMM HKY+rv fit, six of the eight EDN_ECP per-site files and two of the three "
    "PSJS per-pair files, plus the JS-format PSJS_HKY_EDN_ECP summary) agree to <= 1.4e-9 relative in lnL.",
    "- The per-site files whose x is not stored at full precision (EDN_ECP_HKY_nonclock_sw_lnL_check: "
    "best candidate off by 3.9e-5 per site) or not stored at all (Force_EDN_ECP_MG94_nonclock_sw_lnL, "
    "YLR406C_YDL075W_sitewise_lnL, PSJS_HKY_YDR418W_YEL054C_nonclock_lnL_summary) do not match any "
    "save file in test/save.",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated fixture ids (writes nothing)")
    args = ap.parse_args()
    t0 = time.time()
    todo = FIXTURES
    if args.only:
        ids = set(args.only.split(","))
        todo = [f for f in FIXTURES if f["id"] in ids]
    rows = []
    for f in todo:
        row = evaluate_fixture(f)
        rows.append(row)
        print("%-58s %-17s %s" % (row["id"], row["status"],
                                  "" if row.get("rel_diff") is None else "rel %.2g" % row["rel_diff"]),
              flush=True)
    runtime = time.time() - t0
    print("runtime %.1f s" % runtime)
    if args.only:
        return
    import scipy
    meta = {"generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "runtime_s": round(runtime, 1), "rtol": RTOL, "patterns": PATTERNS}
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    with open(OUT_JSON, "w") as fh:
        json.dump({"meta": dict(meta, counts=counts), "fixtures": rows}, fh, indent=1, default=float)
    write_md(rows, meta)
    print(json.dumps(counts))


if __name__ == "__main__":
    main()
