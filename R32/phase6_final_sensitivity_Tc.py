"""
phase6_final_sensitivity_Tc.py -- NIST critical-constant sensitivity
analysis, built on top of phase6_final.py's monotonic-chase warmstart
logic.

*** STATUS (2026-08-11): scaffolding only. NIST_BASE and DEV_PCTS below
are defined but NOT YET WIRED UP -- the loop starting at "Run every
method through the identical run-1/run-2/... chase" below is still the
ORIGINAL phase6 4-method (Helmholtz/NIST/GCGP/SPGP) loop, unchanged from
phase6_final.py. Still to add: a register_case(tc_dev, pc_dev) helper
that builds a uniquely-named entry in p1cev.METHODS from NIST_BASE with
Tc and/or Pc swapped in, three sweep loops (Tc_only, Pc_only,
Tc_Pc_grid) that call the existing solve_one()/run_cubic() chase logic
with that generated name in place of "NIST", and an .xlsx export with
one sheet per sweep (header row per Tc/Pc combination + 4 ambient rows,
COP cell showing the literal string "not feasible" when an ambient
doesn't converge). See BREADCRUMB_07-20.md, "sensitivity analysis scope
finalized" section, for the full spec. ***

Intended sensitivity analysis, once wired up: NIST's critical
temperature (Tc) and/or critical pressure (Pc) are perturbed by
+/-{0.01%, 0.1%, 0.5%, 1%} (plus a 0% baseline) relative to NIST_BASE's
values (Tc=351.3 K, Pc=57.82e5 Pa), in three separate sweeps -- Tc only
(Pc held at baseline), Pc only (Tc held at baseline), and the full 9x9
grid of both varied together. omega and the Shomate/refstate (F/G)
coefficients in NIST_BASE never change in any of the three sweeps. Each
perturbed combination is registered as a new entry in
phase_1_cubic_eos_validation_refstate.METHODS (imported here as p1cev)
and run through the SAME run-1/run-2/... monotonic-trend chase already
implemented below for the original 4 methods:

  Attempt 1: solve normally (no warmstart -- generic initial guess) at
             ALL FOUR ambients (10, 15, 20, 25).
  Check: did every ambient converge, AND is COP monotonically DECREASING
         as Tamb increases (10->15->20->25)? If yes, done -- no more
         attempts needed for that combination.
  If no: Attempt 2 re-solves EACH ambient warmstarted from Attempt 1's
         OWN converged points at that SAME ambient (via
         warmstart_from_state). Any ambient that didn't converge in
         Attempt 1 falls back to a plain, unwarmstarted solve for
         Attempt 2, since there's nothing to warmstart from there.
  Repeat up to MAX_ATTEMPTS = 10, stopping the moment an attempt is
  fully converged AND monotonically decreasing.

This chase logic, and everything below it in this file for now, is
copied unchanged from phase6_final.py (itself a locked copy of
phase6_monotonic_warmstart.py) -- including its Helmholtz-vs-cubic-PR
history (the vapor_compression_unified.py revert, and the later switch
of Helmholtz to vapor_compression.py on AmountBasis.MOLE). See
phase6_final.py's own history in BREADCRUMB_07-20.md for that context;
it is not repeated here since it predates this file's actual purpose.

Author: Shilpa Narasimhan
Support: Claude AI

Date created: 08/11/2026
"""
from pyomo.environ import value, Var, Constraint
from vapor_compression_cubic_refstate import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
import pandas as pd
import os
import phase_1_cubic_eos_validation_refstate as p1cev

# NIST's full parameter set, hardcoded here as literals rather than
# imported live from phase_1_cubic_eos_validation_refstate.METHODS["NIST"]
# -- pins this script's baseline independent of that file, so it can't
# silently drift if the source dict ever changes later. Pc/Tc/omega feed
# the PR-EoS alpha function directly; A-E are the Shomate ideal-gas Cp
# coefficients; F/G are the IIR-calibrated reference-state offsets (see
# phase_1_cubic_eos_validation_refstate.py's own comment on how F/G were
# derived). omega and A-G are held fixed across every sweep below --
# only Tc and/or Pc ever get perturbed.
NIST_BASE = {
    "Pc": 57.82e5, "Tc": 351.3, "omega": 0.2769,
    "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361,
    "F": 25.9022, "G": 84.9031,
}

# Relative deviations (percent, NOT fraction -- confirmed with Shilpa,
# since the fraction reading would set Tc/Pc to 0 at the +/-1 case,
# breaking the PR EoS). Applied to NIST_BASE["Tc"] and/or
# NIST_BASE["Pc"] as base_value * (1 + dev/100). 0 is an explicit
# baseline case included in all three sweeps (Tc_only, Pc_only,
# Tc_Pc_grid), not just implied by omission.
DEV_PCTS = [-1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1]

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3      # J/kg
MAX_ATTEMPTS = 10
METHODS = ["Helmholtz", "NIST", "GCGP", "SPGP"]

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


def register_case(tc_dev, pc_dev):
    """Registers one Tc/Pc sensitivity case as a new entry in
    p1cev.METHODS (the SAME dict object vapor_compression_cubic_refstate.py
    reads from at CubicCycle construction time -- see BREADCRUMB_07-20.md).

    tc_dev, pc_dev: percent deviations (e.g. 1 = +1%, -0.5 = -0.5%)
    applied to NIST_BASE["Tc"]/["Pc"] respectively. Pass 0 for whichever
    one should stay at its baseline value (e.g. register_case(1, 0) for
    a Tc-only sweep point). omega and the Shomate/refstate (A-G)
    coefficients always come straight from NIST_BASE, unperturbed.

    Each call adds a brand-new dict entry rather than mutating NIST_BASE
    or overwriting an existing METHODS entry -- so nothing needs to be
    restored afterward, and no case can accidentally clobber another.

    Returns (name, Tc, Pc): `name` is the unique key to pass as `method`
    into solve_one()/run_cubic(); Tc/Pc are the actual perturbed values
    (in K, Pa) for logging into the report rows."""
    Tc = NIST_BASE["Tc"] * (1 + tc_dev / 100)
    Pc = NIST_BASE["Pc"] * (1 + pc_dev / 100)
    name = f"NIST_Tc{tc_dev:+.2f}pct_Pc{pc_dev:+.2f}pct"
    p1cev.METHODS[name] = {**NIST_BASE, "Tc": Tc, "Pc": Pc}
    return name, Tc, Pc


def run_chase(method_name):
    """Runs the run-1/run-2/... monotonic-trend chase for one method
    name across all 4 ambients. Returns (final_state, attempts_used,
    reached_monotonic)."""
    prior_points_by_amb = {Tamb: None for Tamb in AMBIENTS}
    attempt_results = {}
    for attempt in range(1, MAX_ATTEMPTS + 1):
        attempt_results = {}
        for Tamb in AMBIENTS:
            seed = prior_points_by_amb[Tamb]
            try:
                r = solve_one(method_name, Tamb, warmstart_points=seed)
            except Exception as e:
                r = {"converged": False, "error": f"{type(e).__name__}: {e}", "points": None, "cop": None}
            attempt_results[Tamb] = r

        cop_by_amb = {Tamb: (attempt_results[Tamb]["cop"] if attempt_results[Tamb].get("converged") else None)
                      for Tamb in AMBIENTS}
        monotonic = is_monotonic_decreasing(cop_by_amb)

        for Tamb in AMBIENTS:
            r = attempt_results[Tamb]
            prior_points_by_amb[Tamb] = r["points"] if r.get("converged") else None

        if monotonic:
            return attempt_results, attempt, True
    return attempt_results, MAX_ATTEMPTS, False


def rows_for_case(sheet, tc_dev, pc_dev, Tc, Pc, final_state, attempts_used, reached_monotonic):
    """Header row + 4 ambient rows for one Tc/Pc combination."""
    rows = [{
        "Tc_dev_pct": tc_dev, "Pc_dev_pct": pc_dev, "Tc_K": Tc, "Pc_Pa": Pc,
        "description": f"Tc {tc_dev:+.2f}%, Pc {pc_dev:+.2f}%",
        "attempts_used": attempts_used, "reached_monotonic": reached_monotonic,
        "Tamb_C": None, "COP": None, "converged": None,
    }]
    for Tamb in AMBIENTS:
        r = final_state[Tamb]
        rows.append({
            "Tc_dev_pct": tc_dev, "Pc_dev_pct": pc_dev, "Tc_K": Tc, "Pc_Pa": Pc,
            "description": None, "attempts_used": None, "reached_monotonic": None,
            "Tamb_C": Tamb,
            "COP": r["cop"] if r.get("converged") else "not feasible",
            "converged": r.get("converged", False),
        })
    return rows
# --- Run all three sensitivity sweeps: Tc only, Pc only, full grid. ---
# Each sweep reuses the exact same run_chase() (the run-1/run-2/...
# monotonic-trend chase, unchanged from the original 4-method loop)
# and rows_for_case() (header row + 4 ambient rows) -- only what varies
# on the outside (which Tc/Pc combination gets registered) changes.
tc_only_rows, pc_only_rows, grid_rows = [], [], []

for tc_dev in DEV_PCTS:
    name, Tc, Pc = register_case(tc_dev, 0)
    print(f"[Tc_only] Tc {tc_dev:+.2f}% (Pc baseline) -- chasing...")
    final_state, attempts, reached = run_chase(name)
    tc_only_rows += rows_for_case("Tc_only", tc_dev, 0, Tc, Pc, final_state, attempts, reached)

for pc_dev in DEV_PCTS:
    name, Tc, Pc = register_case(0, pc_dev)
    print(f"[Pc_only] Pc {pc_dev:+.2f}% (Tc baseline) -- chasing...")
    final_state, attempts, reached = run_chase(name)
    pc_only_rows += rows_for_case("Pc_only", 0, pc_dev, Tc, Pc, final_state, attempts, reached)

for tc_dev in DEV_PCTS:
    for pc_dev in DEV_PCTS:
        name, Tc, Pc = register_case(tc_dev, pc_dev)
        print(f"[Tc_Pc_grid] Tc {tc_dev:+.2f}%, Pc {pc_dev:+.2f}% -- chasing...")
        final_state, attempts, reached = run_chase(name)
        grid_rows += rows_for_case("Tc_Pc_grid", tc_dev, pc_dev, Tc, Pc, final_state, attempts, reached)

# --- Export: one .xlsx, three sheets (Tc_only, Pc_only, Tc_Pc_grid). ---
HERE = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(HERE, "phase6_sensitivity_NIST.xlsx")
with pd.ExcelWriter(out_path) as writer:
    pd.DataFrame(tc_only_rows).to_excel(writer, sheet_name="Tc_only", index=False)
    pd.DataFrame(pc_only_rows).to_excel(writer, sheet_name="Pc_only", index=False)
    pd.DataFrame(grid_rows).to_excel(writer, sheet_name="Tc_Pc_grid", index=False)
print(f"\nwrote {out_path}")
