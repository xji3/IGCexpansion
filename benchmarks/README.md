# Benchmarks and baselines (BEAGLE integration, M0)

These are the jsonctmctree baselines that the BEAGLE engine has to beat, from milestone M0
of the BEAGLE integration plan. Results are in `results/`; `results/m0_summary.md` is the
table.

## Setup

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv -e '.[test]'
```

`import IGCexpansion` installs two runtime shims (`IGCexpansion/_jct_compat.py`) that
jsonctmctree 0.2.0 needs on NumPy >= 2 and SciPy >= 1.14. Import IGCexpansion before
jsonctmctree.

## Runs

```bash
.venv/bin/python benchmarks/m0_baselines.py isigc --procs 8
.venv/bin/python benchmarks/m0_baselines.py psigc --procs 10
.venv/bin/python benchmarks/summarize_m0.py
```

- `isigc`: MG94 on YGL062W_YBR218C (3721 states). Needs the IGCTract alignments next to
  this repository, or `IGC_MAFFT_DIR`.
- `psigc`: the PS-IGC objective (PSJSGeneconv, HKY, rate variation) on YDR418W_YEL054C, 1,461 scenes.

Run them on an otherwise idle machine; they use one BLAS thread per process.

## What "fair" means here

- **Compression.** jsonctmctree's `w` observation reduction weights rows after evaluating
  all of them, so `jct_compressed.py` passes the unique rows and their counts
  (`wnnlogl`, `wdnderi`, `wdngrad`). Results agree with the uncompressed calls to
  about 1e-16 relative for lnL and 3e-13 absolute for edge derivatives.
- **O(E) edge pass.** `sdngrad`/`wdngrad` instead of IGCexpansion's per-edge `sdnderi`.
- **Processes.** A pool of single-threaded worker processes, each holding its own model;
  PS-IGC scenes are grouped by compressed size and dispatched largest first.
- **Gradients.** IS-IGC: edge pass plus central differences for the model parameters.
  PS-IGC: edge derivatives plus IGCexpansion's own per-parameter forward differences
  (step 1e-8), as in `loglikelihood_and_gradient_for_one_n`.

## PSJSGeneconv evaluation point

`PSJSGeneconv.__init__` leaves every branch length at 1.0; the rate part of `self.x` is
applied only by the first `unpack_x`. The package's own test values (and
`tests/golden/l2.json` case `psjs_ydr418w`) are at that constructor state.
The baselines call `unpack_x(self.x)` first, so lnL and gradient are at the same point
(golden case `psjs_ydr418w_at_x`). Fitting code always unpacks first, so fits are not
affected. The legacy gradient code, called before any `unpack_x`, differences two
different states and returns nonsense; that path is not used by `objective_and_gradient`.

## Results (2026-11-07, hmm.stat.iastate.edu: 2 x Xeon Gold 6438Y+, 64 cores; Python 3.9, numpy 2.0.2)

Full tables: `results/m0_summary.md`; raw: `results/m0_*.json`, `results/m0_hmm.log`.
The BEAGLE CPU/GPU numbers for the same machine are in the integration plan's `gpu-bench/`.

| Baseline | Wall |
|---|---|
| IS-IGC YGL062W lnL, 1 process (uncompressed / compressed) | 20.1 s / 16.1 s |
| IS-IGC YGL062W lnL + edge derivatives, 1 process (`sdnderi` / `sdngrad` / compressed `wdngrad`) | 69.9 / 40.5 / 32.6 s |
| **IS-IGC YGL062W lnL + full gradient, 8 processes (fair)** | **36.2 s** |
| PS-IGC objective lnL, status quo, 1 process | 104.9 s |
| PS-IGC objective lnL, compressed, 1 / 10 / 32 processes | 80.2 / 8.4 / 2.9 s |
| **PS-IGC lnL + gradient, compressed, 10 / 32 processes (fair)** | **79.2 / 29.8 s** |

- PS-IGC compression is only 1.7x in rows (119,316 -> 69,624) and saves about 25% of the
  objective, because building the 1,461 scenes in Python costs 27.5 s of the 105 s.
- The status-quo PS-IGC lnL + gradient (11,688 uncompressed calls on 1 process) was not
  rerun; the spike extrapolated about 1,175 s.
- On this machine one jsonctmctree MG94 lnL is about 5x slower than on an Apple M5
  (20.1 s vs 3.9 s); the PS-IGC objective is about the same (105 s vs 109 s).
