# Reproducing the PS-IGC tract-length results

This folder reproduces the simulation (arm L) and real-data results for the mean IGC tract length 1/p,
computed with the PS-IGC pair-site engine `IGCexpansion.psigc`. It contains every input, the run scripts and
the report generator.

```
psigc_reproduce/
├── README.md
├── constraints.txt            pinned versions used for the published runs
├── data/
│   ├── armL_bundle/           800 simulated alignments (2017 archive) + truth parameters (manifest.jsonl)
│   └── realdata_inputs/       17 real datasets (yeast pairs, EDN/ECP, HBG, Harpak introns)
├── scripts/
│   ├── armL.py                M2 fits of the simulated datasets
│   ├── uq_run.py              uncertainty (Godambe/CLR intervals) and parametric bootstrap, simulated data
│   ├── realdata.py            fits, intervals and bootstrap, real data
│   ├── recheck_pb_nd.py       optional: exact convergence re-check of bootstrap refits
│   └── report_data.py         collects all results into one JSON
└── report/
    ├── gen_report.py          writes the LaTeX report from that JSON
    └── PSIGC_RESULTS_body.tex report template
```

## 1. Hardware

The published runs used two NVIDIA A100 80 GB GPUs (CUDA 12.8, driver 595) on a 128-core Linux host.
The engine needs double-precision GPU throughput: consumer and workstation GPUs (e.g. RTX A6000) work but are
many times slower. Up to 55 GB of GPU memory was in use during the arm-L intervals, so an 80 GB card is the
tested configuration. The simulation of bootstrap replicates runs on CPU worker processes (`--sim-workers`).

Measured times (one A100 per process unless stated):

| Step | Time |
|---|---|
| Arm-L M2 fits, 800 datasets | 6.5 h on 2 GPUs |
| Arm-L intervals (UQ), 800 datasets | 59 h on 2 GPUs |
| Arm-L bootstrap, 16 datasets | 9 h |
| Real-data fits and intervals, 17 datasets | 6 h |
| Real-data bootstrap, 17 datasets | 38 h |

## 2. Install

Python 3.12 on Linux x86-64 (macOS works for the CPU tests only).

```bash
git clone -b psigc-engine https://github.com/xji3/IGCexpansion.git
cd IGCexpansion
python3.12 -m venv .venv
source .venv/bin/activate
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
pip install -e ".[gpu,test]" -c psigc_reproduce/constraints.txt
```

Check the install (CPU tests; the CUDA-only tests skip without a GPU):

```bash
python -m pytest -q tests/psigc
```

On a GPU host, `tests/psigc/test_gpu.py::test_graph_capture` checks CUDA-graph capture, which has broken
across torch versions; rerun it after any change of torch, CUDA or driver.

All later commands run from `psigc_reproduce/`, with these settings:

```bash
cd psigc_reproduce
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Every runner writes one JSON per task and skips tasks whose JSON exists, so a run can be stopped and restarted.
Failed tasks leave a traceback under `<out>/errors/`.

## 3. Simulation study (arm L)

The 800 datasets are the 2017 archive: 8 cells with true 1/p of 3, 10, 50, 100, 200, 300, 400 and 500,
100 replicates each, simulated under HKY + PS-IGC on the YDR418W/YEL054C yeast tree. `manifest.jsonl` holds the
truth parameters, file hashes and event counts. Replicate 88 is the same alignment at 1/p = 400 and 500; the 2017
simulator yields it in both cells (shared seeds), so it is a genuine draw and is included.

**Step 1: M2 fits** (joint composite-likelihood estimate of all parameters):

```bash
python scripts/armL.py run --bundle data/armL_bundle --out results/armL_M2 \
    --estimators M2 --backend gpu --devices 0,1 --workers-per-gpu 1
```

**Step 2: intervals** (Wald naive and sandwich, CLR with global and local λ, boundary test; J by parametric
bootstrap with R = 100, H_E from 25 replicates):

```bash
python scripts/uq_run.py run --bundle data/armL_bundle --m2 results/armL_M2 --out results/armL_UQ \
    --devices 0,1 --arm 0 --sim-workers 24 --local-he --h-reps 25
```

**Step 3: bootstrap percentile intervals** with full refits, for replicates 1–2 of each cell:

```bash
python scripts/uq_run.py pb --bundle data/armL_bundle --m2 results/armL_M2 --out results/armL_PB \
    --devices 0 --arm 0 --reps 1-2 --sim-workers 24
```

`--cells` and `--reps` (e.g. `--cells 50,500 --reps 1-10`) restrict any step to a subset.

## 4. Real data

The published runs cover 17 datasets: every dataset in `data/realdata_inputs` except HBG_genomic (one
4,277-site locus, far costlier, not analysed).

```bash
NAMES=$(python scripts/realdata.py list --root data/realdata_inputs | cut -d' ' -f1 | grep -v '^HBG_genomic$' | paste -sd, -)
python scripts/realdata.py fit --root data/realdata_inputs --names "$NAMES" --device cuda:0 \
    --out results/realdata_fits --sim-workers 24 --local-he --h-reps 25
python scripts/realdata.py pb --root data/realdata_inputs --names "$NAMES" --device cuda:0 \
    --fits results/realdata_fits --out results/realdata_pb --sim-workers 24
```

`python scripts/realdata.py list --root data/realdata_inputs` lists the datasets; `--names A,B` runs a subset.

Harpak_introns is the only non-coding dataset. The published run labelled its fits "optimizer failure" because
the codon-position rates r₂, r₃ were left free although the likelihood does not depend on them; an exact
Newton-decrement check (`scripts/recheck_pb_nd.py`) showed every fit converged. The code here holds r₂, r₃ fixed
for non-coding data, so a rerun does not need the re-check.

## 5. Tables and report

```bash
python scripts/report_data.py results report_data.json
python report/gen_report.py report_data.json report/PSIGC_RESULTS.tex
cd report && xelatex PSIGC_RESULTS.tex && xelatex PSIGC_RESULTS.tex
```

`report_data.json` holds every number in the report: per cell and method the coverage, width and interval
score; the bootstrap comparison; and the real-data intervals.

## 6. Reproducibility

- **Seeds:** every simulated replicate (bootstrap and J) has its own PCG64 stream from a fixed study entropy
  (`IGCexpansion.psigc.seeds`) and a key of (stage, arm, cell, replicate), so reruns use identical datasets.
- **GPU arithmetic:** reductions on the GPU are not bit-reproducible across hardware or library versions.
  Reruns can differ in the last digits; interval endpoints can move slightly where the profile is flat.
- **Provenance:** each result JSON records the code hash, package versions, GPU and seed key.
