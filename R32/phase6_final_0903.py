"""
phase6_final_GCN.py -- full COP sweep, ALL FIVE methods (Helmholtz, NIST,
GCGP, first_principle, gcn -- SPGP replaced per the 08/26 note below),
each run through phase5_monotonic_warmstart.py's run-1/run-2/... chase:

  Attempt 1: solve the method normally (no warmstart at all -- generic
             initial guess) at ALL FOUR ambients (10, 15, 20, 25).
  Check: did every ambient converge, AND is COP monotonically DECREASING
         as Tamb increases (10->15->20->25)? If yes, that method is done
         -- no more attempts needed.
  If no (any ambient failed, or the trend isn't strictly decreasing):
         Attempt 2 re-solves EACH ambient warmstarted from Attempt 1's
         OWN converged points at that SAME ambient (via
         warmstart_from_state). Any ambient that didn't converge in
         Attempt 1 falls back to a plain, unwarmstarted solve for
         Attempt 2, since there's nothing to warmstart from there.
  Repeat: Attempt 3 seeded by Attempt 2's own results, and so on, up to
          MAX_ATTEMPTS = 10, stopping the moment an attempt is fully
          converged AND monotonically decreasing.

This extends phase5_monotonic_warmstart.py -- which only ran this chase
for GCGP and SPGP, treating Helmholtz and NIST as a fixed, one-shot
"context" reference solved once with no iteration -- to run ALL FOUR
methods through the identical chase. Helmholtz and NIST are expected to
land on Attempt 1 every time (neither has ever shown the branch-jump
instability GCGP/SPGP show), so for them this is mainly a sanity check,
not a fix -- but it means every method in this report went through
exactly the same procedure, with no method treated as special.

Every attempt's full COP-vs-Tamb sweep is saved (not just the final
one), same as phase5, so the progression toward (or failure to reach) a
monotonic trend is visible. On top of that, this script also tracks how
many attempts each method actually needed (attempts_used), since a
method landing on Attempt 1 vs. never converging in 10 is itself a
meaningful result worth reporting, not just the final COP numbers.

Final report format matches phase4_report.py's combined table (COP vs
ambient, with fail:<residual> for non-converged cells, plus a %-vs-
Helmholtz table) -- with the attempts-used line/column added.

Design choice made here, worth flagging: the %-vs-Helmholtz table below
uses THIS script's own live, freshly-computed Helmholtz COP as the
reference (falling back to the hardcoded Phase 3a value only if
Helmholtz itself failed to converge in this run) -- rather than always
using the hardcoded Phase 3a reference the way phase4_common.py does.
This keeps the whole report internally self-consistent (everything
computed fresh, same run), but means the %-vs-Helmholtz numbers here
could differ very slightly from phase4's if Helmholtz's live COP isn't
bit-for-bit the same as the hardcoded 3.95/3.53/3.19/2.91. Flag if you'd
rather pin it to the hardcoded reference instead.

REVERTED (2026-08-11): the vapor_compression_unified.py / molar-basis-
unification approach was set aside -- Shilpa determined the pasted-code
identity check didn't hold up (see BREADCRUMB_07-20.md) and asked to
revert. Back to the original two-file split: NIST/GCGP/SPGP via
vapor_compression_cubic_refstate.py, Helmholtz via a dedicated cycle
file, exactly as this script had it before that detour.

UPDATED (2026-08-11): Helmholtz now imports from vapor_compression.py
(R32/) instead of vapor_compression_plr.py -- this is the file with the
full ambient/approach/evap_sat_temperature spec interface matching
cubic-refstate's, converted to AmountBasis.MOLE so both paths share the
same basis. relax_enth_bounds and _helm_point were updated accordingly
(enth_mol, with the H_MAX*M_R32 basis split).

COPY (2026-08-26): this file is the _GCN working copy of phase6_final.py,
for running the Colon group's two new Shomate-fit methods through the
same monotonic-chase COP sweep as Helmholtz/NIST/GCGP. METHODS (below)
is ["Helmholtz", "NIST", "GCGP", "first_principle", "gcn"], replacing
"SPGP". F/G for every method feeding this sweep are 0.0 (see
phase_1_cubic_eos_validation_refstate_GCN.py's own note) -- confirmed
algebraically to leave every COP number unchanged from what non-zero
F/G would have given, since COP only depends on h/s differences and the
compressor's isentropic equality, both of which cancel any constant
per-method F/G offset.

Bug found and fixed (2026-08-26): this file's CubicCycle import was
still pointing at vapor_compression_cubic_refstate (the production
file), not vapor_compression_cubic_refstate_GCN. The production file's
property package has no "first_principle"/"gcn" keys, so every call
with those method names hit vapor_compression_cubic_refstate.py's own
`assert method in METHODS` and raised immediately -- before any IPOPT
solve even started. run_cubic()'s try/except caught that and reported
it as converged=False, which is why a full run showed "FAILED" at
every ambient, on every one of the 10 monotonic-chase attempts,
identically -- not a numerical non-convergence at all, just wrong
wiring. Fixed by importing from vapor_compression_cubic_refstate_GCN
instead. Not yet re-run after this fix -- first_principle/gcn's earlier
poor fit to Linde (65.99%/46.72% pressure MAPE, see
compare_cp_methods_GCN.py) may still cause genuine IPOPT convergence
trouble even with the import corrected; that would show up now as an
actual max_residual or a different exception, not this assertion.

COPY (2026-09-03): this file is the _0903 working copy of
phase6_final_GCN.py, updated to run first_principle's post-
thermoreconciliation data (spgp_r32 (4).xlsx) through the real cycle.
CubicCycle import repointed to vapor_compression_cubic_refstate_0903
(which itself imports from phase_1_cubic_eos_validation_refstate_0903 --
the file holding today's reconciled Pc=49.62 bar, Tc=392.094 K,
omega=-0.090036 for first_principle, already verified against the
IDAES-vs-vanilla cubic-EOS gate, |dZ|=4.10e-07).

"gcn" ADDED BACK (09/03/2026, later same day): now present in
phase_1_cubic_eos_validation_refstate_0903.py's METHODS dict using the
(4).xlsx thermoreconciled Pc/Tc/Shomate, with omega computed via the
SAME mmHg-reading convention as first_principle (procedural
consistency -- same spreadsheet, same column). Result: omega=-0.88,
well outside any physically normal range. Per Shilpa's explicit
direction, today's run uses ONLY the thermoreconciled numbers, so gcn
is included despite this -- but its result should be treated as
provisional/exploratory, not a validated COP, given the Colon group's
own "4.5% feasible within CI" flag on gcn's data.

Author: Shilpa Narasimhan
Support: Claude AI

Date created: 08/11/2026 (original phase6_final.py)
This _GCN copy created: 08/26/2026
This _0903 copy created: 09/03/2026
"""
from pyomo.environ import value, Var, Constraint
from vapor_compression_cubic_refstate_0903 import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
import pandas as pd
import os

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3      # J/kg
MAX_ATTEMPTS = 10
METHODS = ["Helmholtz", "NIST", "GCGP", "first_principle", "gcn"]

# Phase 3a baseline (PHASE3_NOTES.md Section 1) -- used only as a fallback
# reference if Helmholtz itself fails to converge in this run; otherwise
# this script's own live Helmholtz COP is the reference (see docstring).
HELMHOLTZ_REF = {10: 3.95, 15: 3.53, 20: 3.19, 25: 2.91}

STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]


def relax_enth_bounds(vc, hmax=H_MAX):
    """Widens the enth_mass/enth_mol upper bound wherever it's tighter
    than hmax. Needed on the Helmholtz path specifically -- its default
    bounds can be too tight for some ambients, causing spurious
    infeasibilities unrelated to the actual physics.

    H_MAX (700e3) is a mass-basis value (J/kg). vapor_compression.py
    (R32/) is now AmountBasis.MOLE (enth_mol), so the same numeric
    ceiling is wrong there: 700e3 J/mol is ~19x too large relative to
    R32's real molar enthalpy. The physically equivalent molar ceiling
    is hmax * M_R32 (700e3 * 0.052024 = ~36.4e3 J/mol)."""
    hmax_mol = hmax * M_R32
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        name = v.parent_component().local_name
        if name == "enth_mass" and v.ub is not None and v.ub < hmax:
            v.setub(hmax)
        elif name == "enth_mol" and v.ub is not None and v.ub < hmax_mol:
            v.setub(hmax_mol)


def _cubic_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mol) / M_R32,
        "s": value(state.entr_mol) / M_R32,
        "x": value(state.phase_frac["Vap"]),
    }


def _helm_point(state):
    """Helmholtz counterpart to _cubic_point() -- vapor_compression.py
    (R32/) is now molar basis (AmountBasis.MOLE), same as the cubic-PR
    path, so this now divides by M_R32 identically to _cubic_point()
    to report h/s in mass-basis units (kJ/kg-equivalent) for the CSV/
    report tables. (Previously read enth_mass/entr_mass directly when
    this cycle was mass-basis -- no longer applicable.)"""
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mol) / M_R32,
        "s": value(state.entr_mol) / M_R32,
        "x": value(state.phase_frac["Vap"]),
    }


def _max_constraint_residual(model):
    """Same diagnostic as phase4_common.py -- scans every active
    Constraint for the largest |body - target| residual. Lets a
    non-converged cell report fail:<residual> in the final table instead
    of a bare pass/fail flag, matching phase4_report.py's table style."""
    worst = 0.0
    for con in model.component_data_objects(Constraint, active=True, descend_into=True):
        try:
            body_val = value(con.body, exception=False)
            lower = value(con.lower, exception=False) if con.lower is not None else 0.0
            upper = value(con.upper, exception=False) if con.upper is not None else 0.0
            if body_val is None:
                continue
            target = lower if lower == upper else 0.0
            resid = abs(body_val - target)
        except Exception:
            continue
        if resid is not None and resid > worst:
            worst = resid
    return worst


def warmstart_from_state(model, points):
    """Sets temperature/pressure initial values on a not-yet-solved
    cycle's stream outlets directly from a previous solve's own converged
    `points` dict. Valid under IMPROVED_TPX mode, where temperature and
    pressure are the actual primary state variables, not derived
    expressions."""
    unit_map = {
        "evap_out": model.fs.evaporator,
        "comp_out": model.fs.compressor,
        "cond_out": model.fs.condenser,
        "valve_out": model.fs.expansion_valve,
    }
    for stream, unit in unit_map.items():
        state = unit.control_volume.properties_out[0]
        state.temperature.set_value(points[stream]["T"])
        state.pressure.set_value(points[stream]["P"])


def run_cubic(method, Tamb, warmstart_points=None):
    """One (method, Tamb) solve for a cubic-PR property fit (NIST/GCGP/
    SPGP). warmstart_points, if given, is a previous attempt's own
    converged `points` dict at this SAME Tamb -- applied after
    vc.initialize() and before vc.set_specifications(), so it only
    changes the starting guess, never the target condition."""
    Tcond_sat = Tamb + 9
    vc = CubicCycle(FLUID, compressor_efficiency=0.9999, mode=CubicMode.IMPROVED_TPX, method=method)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
    if warmstart_points is not None:
        warmstart_from_state(vc.model, warmstart_points)
    vc.set_specifications(ambient_temperature=Tamb, condenser_approach=9,
                           evap_sat_temperature=-29, superheating=0, subcooling=0,
                           max_pressure_ratio=10)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    m = vc.model
    points = {
        "evap_out": _cubic_point(m.fs.evaporator.control_volume.properties_out[0]),
        "comp_out": _cubic_point(m.fs.compressor.control_volume.properties_out[0]),
        "cond_out": _cubic_point(m.fs.condenser.control_volume.properties_out[0]),
        "valve_out": _cubic_point(m.fs.expansion_valve.control_volume.properties_out[0]),
    }
    result = {"cop": cop, "converged": converged, "points": points}
    if not converged:
        result["max_residual"] = _max_constraint_residual(m)
    return result


def run_helm(Tamb, warmstart_points=None):
    """Helmholtz counterpart to run_cubic() -- no `method` argument,
    since there's only one Helmholtz package. Accepts warmstart_points,
    mirroring run_cubic exactly, so Helmholtz goes through the same
    monotonic-chase loop as the three cubic-PR methods."""
    Tcond_sat = Tamb + 9
    vc = HelmCycle(FLUID, compressor_efficiency=0.9999, mode=HelmMode.IMPROVED_TPX)
    relax_enth_bounds(vc)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
    if warmstart_points is not None:
        warmstart_from_state(vc.model, warmstart_points)
    vc.set_specifications(ambient_temperature=Tamb, condenser_approach=9,
                           evap_sat_temperature=-29, superheating=0, subcooling=0,
                           max_pressure_ratio=10)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    m = vc.model
    points = {
        "evap_out": _helm_point(m.fs.evaporator.control_volume.properties_out[0]),
        "comp_out": _helm_point(m.fs.compressor.control_volume.properties_out[0]),
        "cond_out": _helm_point(m.fs.condenser.control_volume.properties_out[0]),
        "valve_out": _helm_point(m.fs.expansion_valve.control_volume.properties_out[0]),
    }
    result = {"cop": cop, "converged": converged, "points": points}
    if not converged:
        result["max_residual"] = _max_constraint_residual(m)
    return result


def solve_one(method, Tamb, warmstart_points=None):
    """Dispatches to run_helm or run_cubic depending on method -- the one
    place the call-signature difference (run_helm has no `method` arg)
    is handled, so the chase loop below doesn't need to branch on
    method type itself."""
    if method == "Helmholtz":
        return run_helm(Tamb, warmstart_points=warmstart_points)
    return run_cubic(method, Tamb, warmstart_points=warmstart_points)


def is_monotonic_decreasing(cop_by_amb):
    """cop_by_amb: dict {Tamb: cop or None}. True only if EVERY ambient
    converged (no None values) AND COP strictly decreases as Tamb
    increases, matching the physically-expected trend."""
    cops = [cop_by_amb[Tamb] for Tamb in AMBIENTS]
    if any(c is None for c in cops):
        return False
    return all(cops[i] > cops[i + 1] for i in range(len(cops) - 1))


# --- Run every method through the identical run-1/run-2/... chase. ---
all_attempt_records = []       # one row per (method, attempt, Tamb) -- full history, all methods
final_state_by_method = {}     # method -> {Tamb: result dict}, from whichever attempt stopped the loop
attempts_used = {}             # method -> attempt number the chase stopped/gave up on
reached_monotonic = {}         # method -> True/False

for method in METHODS:
    print(f"\n{'='*90}\n  {method}: iterative monotonic-trend chase (cap = {MAX_ATTEMPTS} attempts)\n{'='*90}")

    prior_points_by_amb = {Tamb: None for Tamb in AMBIENTS}  # warmstart seed per ambient, None = no seed yet
    reached_monotonic[method] = False

    for attempt in range(1, MAX_ATTEMPTS + 1):
        attempt_results = {}
        for Tamb in AMBIENTS:
            seed = prior_points_by_amb[Tamb]
            try:
                r = solve_one(method, Tamb, warmstart_points=seed)
            except Exception as e:
                r = {"converged": False, "error": f"{type(e).__name__}: {e}", "points": None, "cop": None}
            attempt_results[Tamb] = r
            all_attempt_records.append({
                "method": method, "attempt": attempt, "Tamb_C": Tamb,
                "seeded_by": ("none" if seed is None else f"attempt {attempt - 1} own result"),
                "converged": r.get("converged", False),
                "COP": r.get("cop") if r.get("converged") else None,
                "max_residual": r.get("max_residual"),
                "error": r.get("error"),
            })

        cop_by_amb = {Tamb: (attempt_results[Tamb]["cop"] if attempt_results[Tamb].get("converged") else None)
                      for Tamb in AMBIENTS}
        cop_str = ", ".join(f"{Tamb}C={cop_by_amb[Tamb]:.4f}" if cop_by_amb[Tamb] is not None else f"{Tamb}C=FAILED"
                             for Tamb in AMBIENTS)
        monotonic = is_monotonic_decreasing(cop_by_amb)
        print(f"  Attempt {attempt}: {cop_str}  -- {'MONOTONIC (done)' if monotonic else 'not monotonic / has failures'}")

        # Next attempt's warmstart seed, per ambient: this attempt's own
        # converged points where available, otherwise no seed (None) for
        # that ambient next time.
        for Tamb in AMBIENTS:
            r = attempt_results[Tamb]
            prior_points_by_amb[Tamb] = r["points"] if r.get("converged") else None

        if monotonic:
            reached_monotonic[method] = True
            final_state_by_method[method] = attempt_results
            attempts_used[method] = attempt
            print(f"  -> Reached a fully-converged, monotonically decreasing trend on attempt {attempt}.")
            break
    else:
        final_state_by_method[method] = attempt_results
        attempts_used[method] = MAX_ATTEMPTS
        print(f"  -> Did NOT reach a monotonic trend within {MAX_ATTEMPTS} attempts.")

# --- Save every attempt, every ambient, every method, to CSV. ---
HERE = os.path.dirname(os.path.abspath(__file__))
attempts_df = pd.DataFrame(all_attempt_records)
attempts_df.to_csv(os.path.join(HERE, "phase6_monotonic_warmstart_results.csv"), index=False)
print(f"\nSaved {len(all_attempt_records)} rows to phase6_monotonic_warmstart_results.csv")

# --- Combined report, matching phase4_report.py's table format. ---
print(f"\n{'='*100}\n  PHASE 6: COP vs ambient, all methods, after monotonic-chase warmstart\n{'='*100}")
header = f"{'T_amb':>7}"
for m in METHODS:
    header += f"{m:>14}"
print(header)
for Tamb in AMBIENTS:
    line = f"{Tamb:>7}"
    for m in METHODS:
        r = final_state_by_method[m][Tamb]
        if r.get("converged"):
            line += f"{r['cop']:>14.4f}"
        else:
            resid = r.get("max_residual")
            cell = f"fail:{resid:.0f}" if resid is not None else "--"
            line += f"{cell:>14}"
    print(line)

print(f"\n{'='*100}\n  % difference vs Helmholtz reference (this run's own live Helmholtz COP)\n{'='*100}")
header = f"{'T_amb':>7}"
for m in METHODS[1:]:  # skip Helmholtz itself -- it's the reference column
    header += f"{m:>12}"
print(header)
for Tamb in AMBIENTS:
    line = f"{Tamb:>7}"
    helm_r = final_state_by_method["Helmholtz"][Tamb]
    helm_ref = helm_r["cop"] if helm_r.get("converged") else HELMHOLTZ_REF[Tamb]
    for m in METHODS[1:]:
        r = final_state_by_method[m][Tamb]
        if r.get("converged"):
            pct = 100.0 * (r["cop"] - helm_ref) / helm_ref
            cell = f"{pct:.2f}%"
            line += f"{cell:>12}"
        else:
            line += f"{'--':>12}"
    print(line)

print("\nAttempts used to reach a monotonic trend (cap = 10):")
for m in METHODS:
    status = f"attempt {attempts_used[m]}" if reached_monotonic[m] else f"{MAX_ATTEMPTS} (not reached)"
    print(f"  {m:<10} {status}")

# --- Combined summary CSV, one row per ambient, all methods + attempts_used. ---
summary_rows = []
for Tamb in AMBIENTS:
    row = {"ambient_C": Tamb}
    for m in METHODS:
        r = final_state_by_method[m][Tamb]
        row[f"{m}_cop"] = r["cop"] if r.get("converged") else None
        row[f"{m}_converged"] = r.get("converged", False)
        row[f"{m}_max_residual"] = r.get("max_residual")
        row[f"{m}_attempts_used"] = attempts_used[m]
        row[f"{m}_reached_monotonic"] = reached_monotonic[m]
    summary_rows.append(row)
summary_df = pd.DataFrame(summary_rows)
summary_out_path = os.path.join(HERE, "phase6_combined_report.csv")
summary_df.to_csv(summary_out_path, index=False)
print(f"\nwrote {summary_out_path}")
