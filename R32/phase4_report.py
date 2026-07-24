"""
phase4_report.py -- read the three per-method CSVs written by
phase4_nist_sweep.py / phase4_gcgp_sweep.py / phase4_spgp_sweep.py and build
the combined Phase 4 comparison table: COP vs ambient, all three Colon-group
R-32 parameter sets (NIST, GCGP, SPGP) vs the Phase 3a Helmholtz baseline.

Run the three sweep scripts FIRST (each writes its own
phase4_cop_<METHOD>.csv next to this file), then run this script.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import os
import csv
from phase4_common import AMBIENTS, HELMHOLTZ_REF

METHODS = ["NIST", "GCGP", "SPGP"]
HERE = os.path.dirname(os.path.abspath(__file__))


def load_method_csv(method):
    """Return {ambient_C: row_dict} for one method's CSV, or {} if the
    file doesn't exist yet (so the report can still run partially and
    tell the user which sweep(s) still need to be run)."""
    path = os.path.join(HERE, f"phase4_cop_{method}.csv")
    if not os.path.exists(path):
        print(f"  (missing: {path} -- run phase4_{method.lower()}_sweep.py first)")
        return {}
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            out[int(float(row["ambient_C"]))] = row
    return out


data = {m: load_method_csv(m) for m in METHODS}

print(f"\n{'='*100}\n  PHASE 4: COP vs ambient, all methods vs Helmholtz\n{'='*100}")
header = f"{'T_amb':>7}{'Helmholtz':>11}"
for m in METHODS:
    header += f"{m:>12}"
print(header)

for Tamb in AMBIENTS:
    line = f"{Tamb:>7}{HELMHOLTZ_REF[Tamb]:>11.2f}"
    for m in METHODS:
        row = data[m].get(Tamb)
        if row is None:
            line += f"{'(n/a)':>12}"
        elif row["converged"] == "True":
            flag = "*" if row.get("third_attempt") == "True" else ""
            line += f"{float(row['cop']):>11.4f}{flag:>1}"
        else:
            resid = row.get("max_residual", "")
            line += f"{('fail:' + resid[:6]) if resid else '--':>12}"
    print(line)

print("\n(* = only reached via SPGP's third, more permissive solver-tolerance "
      "attempt -- low confidence, see phase4_spgp_sweep.py's note and the "
      "max_residual column in phase4_cop_SPGP.csv)")

# --- % vs Helmholtz table ---
print(f"\n{'='*100}\n  % difference vs Helmholtz reference\n{'='*100}")
header = f"{'T_amb':>7}"
for m in METHODS:
    header += f"{m:>12}"
print(header)
for Tamb in AMBIENTS:
    line = f"{Tamb:>7}"
    for m in METHODS:
        row = data[m].get(Tamb)
        if row is not None and row["converged"] == "True" and row.get("pct_vs_helmholtz"):
            line += f"{row['pct_vs_helmholtz']+'%':>12}"
        else:
            line += f"{'--':>12}"
    print(line)

# --- write a single combined CSV too ---
out_path = os.path.join(HERE, "phase4_combined_report.csv")
fieldnames = ["ambient_C", "helmholtz_cop"] + \
             [f"{m}_cop" for m in METHODS] + \
             [f"{m}_converged" for m in METHODS] + \
             [f"{m}_pct_vs_helmholtz" for m in METHODS] + \
             [f"{m}_max_residual" for m in METHODS]
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for Tamb in AMBIENTS:
        out_row = {"ambient_C": Tamb, "helmholtz_cop": HELMHOLTZ_REF[Tamb]}
        for m in METHODS:
            row = data[m].get(Tamb, {})
            out_row[f"{m}_cop"] = row.get("cop", "")
            out_row[f"{m}_converged"] = row.get("converged", "")
            out_row[f"{m}_pct_vs_helmholtz"] = row.get("pct_vs_helmholtz", "")
            out_row[f"{m}_max_residual"] = row.get("max_residual", "")
        writer.writerow(out_row)
print(f"\nwrote {out_path}")
