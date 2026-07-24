"""
phase4_property_comparison.py -- Phase 4 (extended): compare cubic-PR COP
vs ambient across all three Colon-group R-32 parameter sets (NIST, GCGP,
SPGP), against the Phase 3a Helmholtz baseline sweep (same 4 ambient
points, same ideal-cycle spec).

Sweep (matches phase_3a_helmholtz_cop.py exactly):
  T_amb in [10, 15, 20, 25] C, T_cond_sat = T_amb + 9 (condenser_approach),
  evap_sat_temperature = -29 C, superheating = subcooling = 0 (ideal cycle,
  Shridhar 2016), max_pressure_ratio = 10.

Helmholtz reference (Phase 3a, PHASE3_NOTES.md Section 1):
  T_amb  10    15    20    25
  COP   3.95  3.53  3.19  2.91

Single-point (T_amb=20) cubic-PR result already confirmed this session:
  NIST COP=3.1624 (converged), GCGP COP=3.1899 (converged, essentially
  exact vs Helmholtz), SPGP failed to initialize ("locally infeasible").

SPGP's failure is a real candidate STRUCTURAL issue, not just a numerics/
init bug: its own fitted Pitzer acentric factor is negative (omega=
-0.2741), which makes the PR alpha-function's kappa parameter negative
too (kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2 ~= -0.068) --
unphysical for a real substance (kappa > 0 always, otherwise alpha(T) can
misbehave/go non-monotonic). So SPGP is EXPECTED to fail across some or
all of this sweep too, consistent with Phase 0-2 flagging it as the
"problem child." This script still attempts every (method, ambient) pair
independently (fresh model each time) so a SPGP failure at one ambient
doesn't block NIST/GCGP or other ambients.

FIX (2026-07-24): the first sweep run showed NIST itself (not just SPGP)
crashing HARD at T_amb=15 -- an uncaught InitializationError inside
vc.initialize()'s compressor step (idaes/models/unit_models/pressure_
changer.py's init_isentropic -> properties_out.initialize()), well before
set_specifications()/optimize_COP() are even reached. This is a different,
more fundamental failure point than anything fixed earlier this session
(those were all downstream of a successful vc.initialize()), and it's a
reminder that the earlier fixes were only ever validated at T_amb=20 --
nothing guaranteed they generalize to every ambient in this sweep. Because
this script's whole POINT is to characterize which (method, ambient) cells
converge and which don't (mirroring phase_3a_helmholtz_cop.py's own
fallback-retry, which already treats non-convergence as a data point, not
a fatal error), an uncaught exception killing the ENTIRE remaining sweep
defeats the purpose -- so run_one() is now wrapped in a per-cell
try/except that records the failure (converged=False, cop=None, the
exception type/message) and moves on to the next cell. This is a
deliberate exception to the "no try/except" SOP for exactly this survey
script (not for the core class methods, which still fail loudly).

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import os
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

METHODS_TO_RUN = ["NIST", "GCGP", "SPGP"]
AMBIENTS = [10, 15, 20, 25]  # deg C -- matches phase_3a_helmholtz_cop.py

# Phase 3a Helmholtz baseline, PHASE3_NOTES.md Section 1
HELMHOLTZ_REF = {10: 3.95, 15: 3.53, 20: 3.19, 25: 2.91}


def run_one(method, Tamb):
    Tcond_sat = Tamb + 9
    vc = SimpleVaporCompressionCycle(
        "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX, method=method
    )
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)

    def apply_specs(disable_arc_p):
        vc.set_specifications(
            ambient_temperature=Tamb,
            condenser_approach=9,
            evap_sat_temperature=-29,
            superheating=0,
            subcooling=0,
            max_pressure_ratio=10,
            debug_disable_arc_pressure_eq=disable_arc_p,
        )

    apply_specs(False)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    fallback = False
    if not converged:
        fallback = True
        apply_specs(True)
        cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    return cop, converged, fallback


rows = []
results = {m: {} for m in METHODS_TO_RUN}
for method in METHODS_TO_RUN:
    for Tamb in AMBIENTS:
        print(f"\n{'='*70}\n  METHOD = {method}, T_amb = {Tamb} C\n{'='*70}")
        error_note = ""
        try:
            cop, converged, fallback = run_one(method, Tamb)
        except Exception as e:
            cop, converged, fallback = None, False, None
            error_note = f"{type(e).__name__}: {e}"
            print(f"  CRASHED: {error_note}")
        results[method][Tamb] = {"cop": cop, "converged": converged, "fallback": fallback,
                                  "error": error_note}
        if error_note == "":
            cop_str = f"{cop:.4f}" if converged else "--"
            print(f"  COP = {cop_str}, converged = {converged}, fallback_used = {fallback}")
        rows.append({
            "method": method, "ambient_C": Tamb, "T_cond_sat_C": Tamb + 9,
            "cop": cop if converged else "", "converged": converged,
            "fallback_arc_p_disabled": fallback,
            "helmholtz_ref_cop": HELMHOLTZ_REF[Tamb],
            "pct_vs_helmholtz": round(100.0 * (cop - HELMHOLTZ_REF[Tamb]) / HELMHOLTZ_REF[Tamb], 2)
                                if converged else "",
            "error": error_note,
        })

# --- summary table: ambient down the rows, method+Helmholtz across columns ---
print(f"\n{'='*90}\n  SUMMARY: COP vs ambient, all methods vs Helmholtz\n{'='*90}")
header = f"{'T_amb':>7}{'Helmholtz':>11}"
for m in METHODS_TO_RUN:
    header += f"{m:>10}"
print(header)
for Tamb in AMBIENTS:
    line = f"{Tamb:>7}{HELMHOLTZ_REF[Tamb]:>11.2f}"
    for m in METHODS_TO_RUN:
        r = results[m][Tamb]
        line += f"{r['cop']:>10.4f}" if r["converged"] else f"{'--':>10}"
    print(line)

# --- write CSV next to this script ---
import csv
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "phase4_cop_vs_ambient_all_methods.csv")
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"\nwrote {out_path}")
