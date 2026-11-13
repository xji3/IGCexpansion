"""Study-wide random-number keys (SIMULATION_REDO_PLAN.md sections 4 and 7): every simulated dataset gets
its own PCG64 stream from SeedSequence(STUDY_ENTROPY, spawn_key=(stage, ...)). Stages are never pooled.
STUDY_ENTROPY was drawn once with secrets.randbits(128) on 2026-11-08 and must not change."""
STUDY_ENTROPY = 32077950205718396361442678092716565789

STAGE = {
    'arm_g': 1,          # arm G datasets: (stage, cell index, m index, rep)
    'pilot': 2,          # G5 pilot datasets: (stage, cell index, rep)
    'g3': 3,             # G3 simulator checks
    'j_pb': 10,          # J_PB replicates: (stage, arm, cell index, rep, [restricted-fit flag], r)
    'boot_p_only': 11,   # p-only bootstrap negative control
    'warp': 12,          # warp-speed bootstrap (one replicate per dataset)
    'nested': 13,        # nested bootstrap check (B = 50)
    'j0': 14,            # oracle J0 (scores at theta_0 across replicates)
    'diag_local': 15,    # local-lambda diagnostic: J_PB at restricted fits (stage, cell index, rep, seed, r)
}
