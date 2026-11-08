# IGCexpansion
[![tests](https://github.com/xji3/IGCexpansion/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/xji3/IGCexpansion/actions/workflows/tests.yml)

IGC expansion development folder

##### Dependent software:

[jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation
engine by Alex Griffing)

[Biopython](http://biopython.org/wiki/Biopython)
[networkx](https://networkx.github.io/)
[numpy](https://numpy.org/)
[scipy](https://www.scipy.org/)

(you could install them by

`
pip install --user Biopython networkx numpy scipy
`)


#### Coding Language

Python 3.10 or higher


#### Preparation

*Mac OS / Linux*

1. To install python packages, you need to use [pip](https://pip.pypa.io/en/stable/installing/) (package management).

2. You might need to install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

*Windows*

1. I recommand using [Anaconda](https://www.anaconda.com/products/individual#windows) on Windows that ships with pip functionality and many more useful features.

2. You still need to install [git](https://git-scm.com/download/win)


#### Install python packages

1. Install jsonctmctree package by Alex Griffing (slightly modified for version updates):

`
pip install --user git+https://github.com/xji3/jsonctmctree.git
`

2. Install IGCexpansion package:

`
pip install --user git+https://github.com/xji3/IGCexpansion.git
`

3. Similarly install any other python packages (they should have been installed with IGCexpansion already)

`
pip install --user networkx
`

`
pip install --user Biopython
`

To uninstall:

`
pip uninstall IGCexpansion
`


#### Faster analyses with BEAGLE (optional, recommended)

IGCexpansion can compute likelihoods and gradients with the
[BEAGLE](https://github.com/beagle-dev/beagle-lib) library's *action* plugin, which applies
the sparse rate matrix directly instead of building matrix exponentials. For IS-IGC
analyses (`ReCodonGeneconv`, `JointAnalysis`), this gives:

* **exact gradients** for every parameter. Branch lengths come from BEAGLE's derivative
  kernel; the model parameters (kappa, omega, tau, ...) come from Gauss-Legendre Frechet
  integrals instead of finite differences. They agree with the legacy finite differences to
  about 1e-8 and with an exact reference to 1e-14.
* **speed, especially on a GPU.** On the 1,129-codon yeast pair YGL062W/YBR218C (MG94), one
  log-likelihood plus full gradient takes 36 s with the legacy code on 8 processes, 17 s with
  8 BEAGLE threads, 1.4 s on an RTX A6000 and 0.7 s on an A100.
* **the same answers.** An MG94 fit reaches the legacy optimum (difference 2e-8 in lnL) in
  the same number of iterations.

BEAGLE is used **by default whenever it is installed**. Without it, IGCexpansion prints one
warning and runs the original jsonctmctree code, so existing scripts keep working.

**Install BEAGLE** (the `action` branch of beagle-dev, commit `07fc8b6`; needs CMake and a C++17
compiler; Eigen 3.4 is downloaded automatically if it is not installed; add `-DBUILD_CUDA=ON`
and a CUDA toolkit for GPU support):

```bash
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib && git checkout 07fc8b6 && cd ..
cmake -S beagle-lib -B beagle-build -DBUILD_ACTION=ON -DBUILD_JNI=OFF -DBUILD_OPENCL=OFF -DBUILD_CUDA=OFF \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/opt/beagle-action-07fc8b6
cmake --build beagle-build -j 8 && cmake --install beagle-build
```

IGCexpansion looks for BEAGLE in `~/opt/beagle-action-07fc8b6`. For any other location, set
`IGC_BEAGLE_PREFIX` to the install prefix.

**Settings** (environment variables, set before importing IGCexpansion):

| Variable | Default | Meaning |
|---|---|---|
| `IGC_BEAGLE_PREFIX` | `~/opt/beagle-action-07fc8b6` | BEAGLE install prefix |
| `IGC_BEAGLE_THREADS` | `1` | CPU threads per BEAGLE instance (use threads before processes) |
| `IGC_BEAGLE_RESOURCE` | `cpu` | `gpu` runs on a CUDA GPU (needs a CUDA build of BEAGLE) |
| `IGC_RECODON_ENGINE` | `beagle` | `legacy` uses the original jsonctmctree code for IS-IGC |
| `IGC_RECODON_PROCS` | `1` | more than 1 splits one gene pair's gradient over worker processes |
| `IGC_PSJS_ENGINE` | `legacy` | `beagle` runs the PS-IGC (`PSJSGeneconv`) objective on BEAGLE workers |

For example, on a GPU machine:

```bash
IGC_BEAGLE_RESOURCE=gpu python my_analysis.py
```

`JointAnalysis` (several gene pairs at once) runs on a persistent pool of worker processes,
one per pair by default. Scripts that use worker processes (`JointAnalysis`,
`IGC_RECODON_PROCS`, or `IGC_PSJS_ENGINE=beagle`) must be run as a `.py` file with an
`if __name__ == "__main__":` guard. Expected-count summaries (`get_summary`) still use
jsonctmctree.

The [BEAGLE tutorial](tutorials/IS_IGC_BEAGLE/Run_IS_IGC_BEAGLE.py) fits IS-IGC to the EDN/ECP
example data and, with `--compare`, checks one gradient against the legacy code:

```bash
cd IGCexpansion/tutorials/IS_IGC_BEAGLE
python Run_IS_IGC_BEAGLE.py --compare              # MG94 on 1 CPU thread
python Run_IS_IGC_BEAGLE.py --threads 8            # 8 BEAGLE threads
python Run_IS_IGC_BEAGLE.py --gpu                  # CUDA build of BEAGLE
```

Expected result (MG94): one gradient check with `max |grad diff|` around 3e-7 (the legacy
finite-difference error), then `fit: lnL = -1700.518782` with tau = 0.6312, kappa = 2.0626,
omega = 0.8271. The fit takes about 90 s with 4 threads on an Apple M5 and 23 s on an RTX A6000.
Fitted parameters go to `save/` and the summary to `summary/`.

Development details, tests and benchmarks: `tests/` (`pytest`), `benchmarks/README.md`.

##### Getting a local copy of the package

`
git clone https://github.com/xji3/IGCexpansion.git
`

You can now run the tutorial file or edit it to perform analyses.

`
cd IGCexpansion/tutorials/IS_IGC
`

`
python Run_IS_IGC.py
`


##### Tutorials

For independent site IGC (IS-IGC) analyses, please refer to this [tutorial](tutorials/IndependentSiteTutorial.md).  
There are two versions of the software that implement the same IS-IGC approach.  The difference between them is the flexibility of considering different duplication/loss histories where the first version assumes one single duplication event without loss which the second version does not.

To run IS-IGC with BEAGLE (exact gradients, optional GPU), see the [BEAGLE tutorial](tutorials/IS_IGC_BEAGLE/Run_IS_IGC_BEAGLE.py) and the BEAGLE section above.
