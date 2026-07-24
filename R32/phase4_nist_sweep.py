"""
phase4_nist_sweep.py -- Phase 4, NIST parameter set only: run the ideal
vapor-compression cycle (cubic PR EoS) across the full ambient sweep
(10/15/20/25 C) and write phase4_cop_NIST.csv.

See phase4_common.py for the shared sweep procedure/spec and
phase4_report.py for the combined comparison table (all methods vs the
Phase 3a Helmholtz baseline).

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import os
from phase4_common import run_sweep, write_csv

METHOD = "NIST"

rows = run_sweep(METHOD)
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        f"phase4_cop_{METHOD}.csv")
write_csv(rows, out_path)
