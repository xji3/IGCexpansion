# M0 fair baselines (jsonctmctree)

## IS-IGC: MG94 YGL062W_YBR218C, 3721 states

Linux-5.14.0-687.52.1.el9_8.x86_64-x86_64-with-glibc2.34, 128 CPUs, Python 3.9.25, numpy 2.0.2; 2026-11-07 10:49. 1129 codons -> 1013 unique patterns; get_scene 0.00 s.

| Request (1 process) | Rows | Wall (best of repeats) |
|---|---|---|
| lnL | 1129 | 20.1 s |
| lnL_compressed | 1013 | 16.1 s |
| lnL+sdnderi | 1129 | 69.9 s |
| lnL+sdngrad | 1129 | 40.5 s |
| lnL+wdngrad_compressed | 1013 | 32.6 s |

Compressed lnL relative difference: 0.0e+00.

**Fair full gradient** (compressed lnL + O(E) edge pass + 12 central differences, 8 processes): **36.2 s**.

## PS-IGC: PSJSGeneconv HKY One rate, rate variation, YDR418W_YEL054C

Linux-5.14.0-687.52.1.el9_8.x86_64-x86_64-with-glibc2.34, 128 CPUs, Python 3.9.25, numpy 2.0.2; 2026-11-07 10:56. 1461 scenes; 119316 observation rows -> 69624 unique (1.7x).

| Objective | Processes | Wall | lnL |
|---|---|---|---|
| lnL, status quo (uncompressed) | 1 | 104.9 s | -941165.536943 |
| lnL, compressed | 1 | 80.2 s | -941165.536943 |
| lnL, compressed | 10 | 8.40 s | -941165.536943 |
| lnL + gradient, compressed (11688 calls) | 10 | 79.2 s | -941165.536943 |

Building all scenes in Python (get_scene, 1 process): 27.5 s.

Gradient check against `loglikelihood_and_gradient_for_one_n` (compressed vs status quo):

- n = 1: lnL diff 0.0e+00, max gradient diff 4.5e-05 (max |g| 157)
- n = 488: lnL diff 0.0e+00, max gradient diff 0.0e+00 (max |g| 2.95)

## PS-IGC: PSJSGeneconv HKY One rate, rate variation, YDR418W_YEL054C

Linux-5.14.0-687.52.1.el9_8.x86_64-x86_64-with-glibc2.34, 128 CPUs, Python 3.9.25, numpy 2.0.2; 2026-11-07 11:03. 1461 scenes; 119316 observation rows -> 69624 unique (1.7x).

| Objective | Processes | Wall | lnL |
|---|---|---|---|
| lnL, compressed | 1 | 80.9 s | -941165.536943 |
| lnL, compressed | 32 | 2.93 s | -941165.536943 |
| lnL + gradient, compressed (11688 calls) | 32 | 29.8 s | -941165.536943 |

Building all scenes in Python (get_scene, 1 process): 28.1 s.

Gradient check against `loglikelihood_and_gradient_for_one_n` (compressed vs status quo):

- n = 1: lnL diff 0.0e+00, max gradient diff 4.5e-05 (max |g| 157)
- n = 488: lnL diff 0.0e+00, max gradient diff 0.0e+00 (max |g| 2.95)
