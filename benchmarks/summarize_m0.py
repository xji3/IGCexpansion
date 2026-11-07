"""Build benchmarks/results/m0_summary.md (the fair-baseline table) from the m0_*.json files."""

import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def f(s):
    return "%.1f s" % s if s >= 10 else "%.2f s" % s


def main():
    out = ["# M0 fair baselines (jsonctmctree)", ""]
    for path in sorted(glob.glob(os.path.join(RES, "m0_isigc_*.json"))):
        R = json.load(open(path))
        e, o, p = R["env"], R["one_process"], R["full_gradient_pool"]
        out += ["## IS-IGC: %s" % R["case"], "",
                "%s, %d CPUs, Python %s, numpy %s; %s. %d codons -> %d unique patterns; get_scene %s."
                % (e["platform"], e["cpu_count"], e["python"], e["numpy"], e["date"], R["codons"], R["patterns"],
                   f(R["sec_get_scene"])), "",
                "| Request (1 process) | Rows | Wall (best of repeats) |", "|---|---|---|"]
        for k, v in o.items():
            out.append("| %s | %d | %s |" % (k, v["rows"], f(v["sec"])))
        out += ["", "Compressed lnL relative difference: %.1e." % R["lnL_rel_diff_compressed"], "",
                "**Fair full gradient** (compressed lnL + O(E) edge pass + %d central differences, %d processes): **%s**."
                % (p["tasks"] - 1, p["procs"], f(p["sec"])), ""]
    for path in sorted(glob.glob(os.path.join(RES, "m0_psigc_*.json"))):
        R = json.load(open(path))
        e = R["env"]
        out += ["## PS-IGC: %s" % R["case"], "",
                "%s, %d CPUs, Python %s, numpy %s; %s. %d scenes; %d observation rows -> %d unique (%.1fx)."
                % (e["platform"], e["cpu_count"], e["python"], e["numpy"], e["date"], R["scenes"], R["rows_total"],
                   R["rows_compressed_total"], R["rows_total"] / R["rows_compressed_total"]), "",
                "| Objective | Processes | Wall | lnL |", "|---|---|---|---|"]
        if "status_quo_lnL_1proc" in R:
            v = R["status_quo_lnL_1proc"]
            out.append("| lnL, status quo (uncompressed) | 1 | %s | %.6f |" % (f(v["sec"]), v["lnL"]))
        v = R["compressed_lnL_1proc"]
        out.append("| lnL, compressed | 1 | %s | %.6f |" % (f(v["sec"]), v["lnL"]))
        v = R["compressed_lnL_pool"]
        out.append("| lnL, compressed | %d | %s | %.6f |" % (e["procs"], f(v["sec"]), v["lnL"]))
        v = R["compressed_lnL_and_gradient_pool"]
        out.append("| lnL + gradient, compressed (%d calls) | %d | %s | %.6f |"
                   % (v["calls"], e["procs"], f(v["sec"]), v["lnL"]))
        out += ["", "Building all scenes in Python (get_scene, 1 process): %s." % f(R["sec_get_scene_all_1proc"]), ""]
        out.append("Gradient check against `loglikelihood_and_gradient_for_one_n` (compressed vs status quo):")
        out.append("")
        for n, c in R["gradient_check_vs_status_quo"].items():
            out.append("- n = %s: lnL diff %.1e, max gradient diff %.1e (max |g| %.3g)"
                       % (n, c["lnL_abs_diff"], c["grad_max_abs_diff"], c["grad_max_abs"]))
        out.append("")
    open(os.path.join(RES, "m0_summary.md"), "w").write("\n".join(out))
    print("\n".join(out))


if __name__ == "__main__":
    main()
