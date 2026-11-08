"""IS-IGC fit (ReCodonGeneconv) with the BEAGLE engine and exact gradients.

Uses the EDN/ECP example data of ../IS_IGC_Old. Run from this folder:

    python Run_IS_IGC_BEAGLE.py                      # MG94, BEAGLE on 1 CPU thread
    python Run_IS_IGC_BEAGLE.py --threads 8          # 8 BEAGLE threads
    python Run_IS_IGC_BEAGLE.py --gpu                # CUDA build of BEAGLE
    python Run_IS_IGC_BEAGLE.py --compare            # also time one legacy (jsonctmctree) gradient

BEAGLE is the default IS-IGC engine. The settings below are environment variables, so they
can also be set in the shell instead (IGC_BEAGLE_PREFIX, IGC_BEAGLE_THREADS,
IGC_BEAGLE_RESOURCE, IGC_RECODON_ENGINE=legacy to use the old code). They must be set before
IGCexpansion is imported.
"""

import argparse
import os
import time

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="MG94", choices=["MG94", "HKY"])
    ap.add_argument("--threads", type=int, default=1, help="BEAGLE CPU threads")
    ap.add_argument("--gpu", action="store_true", help="run BEAGLE on the GPU")
    ap.add_argument("--compare", action="store_true", help="also time one legacy lnL + gradient")
    a = ap.parse_args()

    os.environ["IGC_BEAGLE_THREADS"] = str(a.threads)
    if a.gpu:
        os.environ["IGC_BEAGLE_RESOURCE"] = "gpu"

    from IGCexpansion.CodonGeneconv import ReCodonGeneconv
    from IGCexpansion.beagle import recodon

    if not recodon.available():
        raise SystemExit("BEAGLE was not found: set IGC_BEAGLE_PREFIX to its install prefix (see README).")

    data = os.path.join("..", "IS_IGC_Old")
    os.makedirs("save", exist_ok=True)
    os.makedirs("summary", exist_ok=True)
    save_name = os.path.join("save", a.model + "_EDN_ECP_nonclock_beagle_save.txt")
    obj = ReCodonGeneconv(os.path.join(data, "EDN_ECP_tree.newick"), os.path.join(data, "EDN_ECP_Cleaned.fasta"),
                          ["EDN", "ECP"], Model=a.model, Force=None, clock=None,
                          save_path="./save/", save_name=save_name)

    # One objective call: f = -lnL and its exact gradient (x is on the log scale).
    t0 = time.perf_counter()
    f, g = obj.loglikelihood_and_gradient()
    print("BEAGLE  lnL = %.6f  (%.2f s for lnL + exact gradient; first call includes set-up)"
          % (-f, time.perf_counter() - t0))
    if a.compare:
        os.environ["IGC_RECODON_ENGINE"] = "legacy"          # read at every call
        t0 = time.perf_counter()
        f0, g0 = obj.loglikelihood_and_gradient()
        print("legacy  lnL = %.6f  (%.2f s, finite differences)  max |grad diff| = %.2e"
              % (-f0, time.perf_counter() - t0, np.max(np.abs(np.asarray(g0) - g))))
        del os.environ["IGC_RECODON_ENGINE"]

    # Maximum likelihood fit (L-BFGS-B with exact gradients), then the usual summaries.
    t0 = time.perf_counter()
    result = obj.get_mle(display=False)
    print("fit: lnL = %.6f in %d evaluations, %.1f s" % (-result.fun, result.nfev, time.perf_counter() - t0))
    print("tau = %.4f, kappa = %.4f, omega = %.4f" % (obj.tau, obj.kappa, obj.omega))
    obj.get_ExpectedNumGeneconv()
    obj.get_individual_summary(summary_path="./summary/")


if __name__ == "__main__":
    main()
