"""
phase4_gcgp_sweep.py -- Phase 4, GCGP parameter set only: run the ideal
vapor-compression cycle (cubic PR EoS) across the full ambient sweep
(10/15/20/25 C) and write phase4_cop_GCGP.csv.

See phase4_common.py for the shared sweep procedure/spec and
phase4_report.py for the combined comparison table (all methods vs the
Phase 3a Helmholtz baseline).

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24

Debug mode (2026-07-27, added while investigating the T_amb=20 anomaly,
task #37): run with `--diagnose` to print the compressor's real-outlet
vs isentropic T/h (and evap/cond superheat/subcool) at EVERY ambient in
this same sweep, e.g.:
    python3 phase4_gcgp_sweep.py --diagnose
This reuses the exact same case-building code as the real sweep (via
phase4_common.py's diagnose= option) instead of re-deriving the T_amb=20
case in a separate script -- so what you see is guaranteed to be the
same run that produces phase4_cop_GCGP.csv, just with the internal
compressor state exposed for all 4 ambients at once (not just the one
point compressor_fix_regression_check.py has been checking).
"""
import os
import sys
from phase4_common import run_sweep, write_csv

METHOD = "GCGP"
DIAGNOSE = "--diagnose" in sys.argv

rows = run_sweep(METHOD, diagnose=DIAGNOSE)
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        f"phase4_cop_{METHOD}.csv")
write_csv(rows, out_path)
