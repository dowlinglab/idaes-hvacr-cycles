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

COPY (2026-09-09): this file is meant to become the _0909 working copy
of phase6_final_0903.py, to run the (5).xlsx delivery's data through
the real cycle -- the Sep-6-revised first_principle/gcn reconciliation,
plus two brand-new methods, first_principle_multi and gcn_multi, from
the new "multi after recon" section ("multi-output" not yet defined by
the collaborators; question asked, reply pending).

FIXED (2026-09-09): this file had one real bug, now corrected -- the
module docstring was not properly closed (ended in two double-quote
characters instead of three), an actual SyntaxError that meant the file
couldn't run at all.

2026-09-09 (later): the import/METHODS/column-width changes below were
applied -- CubicCycle import repointed to vapor_compression_cubic_
refstate_0909, METHODS extended to all 7 entries (Helmholtz, NIST,
GCGP, first_principle, gcn, first_principle_multi, gcn_multi), and the
report table's column widths widened from 14/12/10 to 24 chars
throughout (main COP table, %-vs-Helmholtz table, attempts-used line)
to fit "first_principle_multi" (22 characters) without adjacent columns
running together.

--- 09/10/2026 (this _0910 copy) ---
New spreadsheet delivery (spgp_r32 (7).xlsx) -- first_principle, gcn,
and first_principle_multi all got revised Shomate a-e fits this round;
gcn_multi's fit failed and falls back to the 09/06 data (full details
in phase_1_cubic_eos_validation_refstate_0910.py's own 09/10 note).
CubicCycle import repointed to vapor_compression_cubic_refstate_0910
(which itself imports from phase_1_cubic_eos_validation_refstate_0910).
METHODS list unchanged (already had all 7 entries from 09/09).

Updated omega values (all via the mmHg-reading omega_from_psat
convention): first_principle -0.235212 (was -0.209932), gcn +0.450573
(was -0.426564 -- sign flip, now closest to real R32's literature omega
~0.277 of any gcn revision so far), first_principle_multi 0.157137
(was 0.060099, still comfortably physical), gcn_multi unchanged at
1.254451 (still on the 09/06 fallback data, still unphysical). Treat
gcn and gcn_multi results as provisional/exploratory, not validated COP
numbers, pending the Colon group's own confirmation of their
reconciled data.

Not yet run with real IDAES as of this writing -- needs Shilpa's own
machine (myidaesenv). p-h/T-s diagrams for this same (7).xlsx data were
generated first, via compare_cp_methods_0910.py (standalone, no IDAES
needed), per the requested order: diagrams before COP.

Author: Shilpa Narasimhan
Support: Claude AI

Date created: 08/11/2026 (original phase6_final.py)
This _GCN copy created: 08/26/2026
This _0903 copy created: 09/03/2026
This _0909 copy created: 09/09/2026
This _0910 copy created: 09/10/2026
"""
from pyomo.environ import value, Var, Constraint
from vapor_compression_cubic_refstate_R1234zeZ_1008 import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
# HELMHOLTZ IS NOT AVAILABLE FOR THIS FLUID.
# IDAES's general_helmholtz package registers 'r1234ze' -- but that is
# R-1234ze(E), the TRANS isomer, Tc = 382.513 K. Our fluid is the CIS
# isomer, R-1234ze(Z), Tc = 423.27 K. Both have molar mass 114.0416 g/mol,
# so a wrong-isomer run would look entirely plausible and be entirely wrong.
# Checked 10/08/2026: registered_components() returns
#   butane, co2, h2o, isobutane, nh3, propane, r1234ze, r125, r134a,
#   r227ea, r32
# -- R-1234ze(Z) is not among them.
#
# The Helmholtz reference for this fluid therefore comes from CoolProp
# (Akasaka & Lemmon), not from IDAES. See compare_cp_methods_R1234zeZ_1008.py
# and verify_akasaka_lemmon_SI.py. run_helm() below is left in place but is
# never reached, because "Helmholtz" is not in METHODS.
try:
    from vapor_compression import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
except Exception as _e:          # pragma: no cover
    HelmCycle = HelmMode = None
    print(f"note: IDAES Helmholtz cycle unavailable ({_e}); not needed for "
          f"R-1234ze(Z).")
import pandas as pd
import os

FLUID = "R1234ze(Z)"
# =========================================================================
# OPERATING CONDITIONS -- the only place the duty is defined.
# =========================================================================
# R-32 ran T_evap = -29 C with T_cond = T_amb + 9 and ambients 10/15/20/25.
# That point is UNUSABLE for R-1234ze(Z). Its saturation pressure at -29 C
# is 0.169 bar -- about 0.17 atm -- and the chain failed to converge there
# at every ambient and every solver setting tried (16 attempts, 10/08/2026:
# default, tol=1e-6, tol=1e-5 + constr_viol_tol, gradient-based scaling).
#
# R-1234ze(Z) is a HIGH-TEMPERATURE fluid: NBP 282.878 K (9.7 C), Tc 423.27 K
# (150.1 C). Petr & Raabe, Energy 93 (2015) 266-274, run it between a
# near-ambient condenser and a 100-250 C heat source with a 10 K pinch; the
# window below sits inside that. The chain converges 4/4 here.
#
# AMBIENTS keeps its name for a minimal diff, but for a high-temperature heat
# pump it is the heat-SINK temperature, not the weather.
AMBIENTS = [60, 70, 80, 90]   # heat-sink temperature [C]   (was 10/15/20/25)
T_EVAP_C = 30.0               # evaporator sat. temperature [C]   (was -29)
COND_APPROACH = 10.0          # condenser approach [K]  (was 9; 10 K matches
                              #   the Petr & Raabe pinch)
SUPERHEAT = 0                 # K   (unchanged)
SUBCOOL = 0                   # K   (unchanged)
# => T_cond = 70 / 80 / 90 / 100 C
M_R32 = 0.1140416  # kg/mol, R-1234ze(Z) (name kept for minimal diff)
H_MAX = 700e3      # J/kg
MAX_ATTEMPTS = 10
METHODS = ["Helmholtz", "PR"]
# Two methods. "Helmholtz" MUST stay first -- the report tables below use
# METHODS[0] as the reference column and METHODS[1:] as the comparison
# columns. Helmholtz here is CoolProp, not IDAES (see run_helm and the note
# at the import above). "PR" is the cubic, solved in IDAES.
#
# When the Colon group delivers R-1234ze(Z) parameter sets, add their names
# to this list AND their rows to METHODS_RAW in
# phase_1_cubic_eos_validation_refstate_R1234zeZ_1008.py.

# Phase 3a baseline (PHASE3_NOTES.md Section 1) -- used only as a fallback
# reference if Helmholtz itself fails to converge in this run; otherwise
# this script's own live Helmholtz COP is the reference (see docstring).
HELMHOLTZ_REF = {60: 6.2386, 70: 4.6539, 80: 3.5623, 90: 2.7473}   # CoolProp, T_evap = 30 C

STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]


def relax_enth_bounds(vc, hmax=H_MAX):
    """ Widens the enth_mass/enth_mol upper bound wherever it's tighter
    than hmax. Needed on the Helmholtz path specifically -- its default
    bounds can be too tight for some ambients, causing spurious
    infeasibilities unrelated to the actual physics.

    H_MAX (700e3) is a mass-basis value (J/kg). vapor_compression.py
    (R32/) is now AmountBasis.MOLE (enth_mol), so the same numeric
    ceiling is wrong there: 700e3 J/mol is ~19x too large relative to
    R32's real molar enthalpy. The physically equivalent molar ceiling
    is hmax * M_R32 (700e3 * 0.052024 = ~36.4e3 J/mol). """
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
    Tcond_sat = Tamb + COND_APPROACH
    vc = CubicCycle(FLUID, compressor_efficiency=0.9999, mode=CubicMode.IMPROVED_TPX, method=method)
    vc.specify_initial_conditions(low_side_temperature=T_EVAP_C, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
    if warmstart_points is not None:
        warmstart_from_state(vc.model, warmstart_points)
    vc.set_specifications(ambient_temperature=Tamb, condenser_approach=COND_APPROACH,
                           evap_sat_temperature=T_EVAP_C, superheating=SUPERHEAT, subcooling=SUBCOOL,
                           max_pressure_ratio=1e6)
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
    """Helmholtz reference for R-1234ze(Z) -- via CoolProp, NOT via IDAES.

    IDAES's general_helmholtz has no R-1234ze(Z): its 'r1234ze' is the E
    (trans) isomer, Tc = 382.513 K, against our 423.27 K. Same molar mass,
    so a wrong-isomer run would look plausible and be wrong. CoolProp
    carries the Akasaka & Lemmon equation, checked against that paper's own
    Supporting Information by verify_akasaka_lemmon_SI.py.

    This is therefore CLOSED FORM, not a solve. Consequences, all
    deliberate:
      - it cannot fail, so converged is always True;
      - there is nothing to warmstart, so warmstart_points is accepted and
        ignored (keeps solve_one()'s dispatch free of a special case);
      - the monotonic-chase loop will finish it on attempt 1.

    Same cycle definition as run_cubic(), so the two columns are comparable:
    saturated vapour at the evaporator outlet, isentropic compression to the
    condensing pressure at the stated efficiency, saturated liquid at the
    condenser outlet, isenthalpic valve. Superheat and subcool are read from
    SUPERHEAT / SUBCOOL and must stay 0 for this closed form to be correct;
    if they are ever made non-zero, this function needs reworking and will
    say so rather than quietly return a wrong number.
    """
    import CoolProp.CoolProp as CP

    if SUPERHEAT != 0 or SUBCOOL != 0:
        raise NotImplementedError(
            f"run_helm() assumes saturated states (SUPERHEAT=SUBCOOL=0); "
            f"got SUPERHEAT={SUPERHEAT}, SUBCOOL={SUBCOOL}. Rework the "
            f"closed form before trusting this column.")

    eta = 0.9999                      # same as run_cubic's compressor_efficiency
    Te = T_EVAP_C + 273.15
    Tc = (Tamb + COND_APPROACH) + 273.15

    p_e = CP.PropsSI("P", "T", Te, "Q", 1, FLUID)
    p_c = CP.PropsSI("P", "T", Tc, "Q", 0, FLUID)

    h1 = CP.PropsSI("H", "T", Te, "Q", 1, FLUID)          # 1 evaporator outlet
    s1 = CP.PropsSI("S", "T", Te, "Q", 1, FLUID)
    h2s = CP.PropsSI("H", "P", p_c, "S", s1, FLUID)       # isentropic discharge
    h2 = h1 + (h2s - h1) / eta                            # 2 compressor outlet
    T2 = CP.PropsSI("T", "P", p_c, "H", h2, FLUID)
    s2 = CP.PropsSI("S", "P", p_c, "H", h2, FLUID)
    h3 = CP.PropsSI("H", "T", Tc, "Q", 0, FLUID)          # 3 condenser outlet
    s3 = CP.PropsSI("S", "T", Tc, "Q", 0, FLUID)
    h4 = h3                                               # 4 valve: isenthalpic
    s4 = CP.PropsSI("S", "P", p_e, "H", h4, FLUID)
    x4 = CP.PropsSI("Q", "P", p_e, "H", h4, FLUID)

    cop = (h1 - h4) / (h2 - h1)

    # h and s in J/kg and J/kg/K, matching _cubic_point (which divides the
    # molar values by M_R32).
    points = {
        "evap_out":  {"T": Te, "P": p_e, "h": h1, "s": s1, "x": 1.0},
        "comp_out":  {"T": T2, "P": p_c, "h": h2, "s": s2, "x": 1.0},
        "cond_out":  {"T": Tc, "P": p_c, "h": h3, "s": s3, "x": 0.0},
        "valve_out": {"T": Te, "P": p_e, "h": h4, "s": s4, "x": x4},
    }
    return {"cop": cop, "converged": True, "points": points}


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
attempts_df.to_csv(os.path.join(HERE, "phase6_monotonic_warmstart_results_1001.csv"), index=False)
print(f"\nSaved {len(all_attempt_records)} rows to phase6_monotonic_warmstart_results_1001.csv")

# --- Combined report, matching phase4_report.py's table format. ---
print(f"\n{'='*100}\n  PHASE 6: COP vs ambient, all methods, after monotonic-chase warmstart\n{'='*100}")
# Column width widened to 24 (from 14/12/10) -- "first_principle_multi" is
# 22 characters, wider than the old fields, which caused adjacent columns
# to run together with no separating space (the same "GCGPfirst_principle"
# merge bug seen in an earlier run's output).
header = f"{'T_amb':>7}"
for m in METHODS:
    header += f"{m:>24}"
print(header)
for Tamb in AMBIENTS:
    line = f"{Tamb:>7}"
    for m in METHODS:
        r = final_state_by_method[m][Tamb]
        if r.get("converged"):
            line += f"{r['cop']:>24.4f}"
        else:
            resid = r.get("max_residual")
            cell = f"fail:{resid:.0f}" if resid is not None else "--"
            line += f"{cell:>24}"
    print(line)

print(f"\n{'='*100}\n  % difference vs Helmholtz reference (this run's own live Helmholtz COP)\n{'='*100}")
header = f"{'T_amb':>7}"
for m in METHODS[1:]:  # skip Helmholtz itself -- it's the reference column
    header += f"{m:>24}"
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
            line += f"{cell:>24}"
        else:
            line += f"{'--':>24}"
    print(line)

print("\nAttempts used to reach a monotonic trend (cap = 10):")
for m in METHODS:
    status = f"attempt {attempts_used[m]}" if reached_monotonic[m] else f"{MAX_ATTEMPTS} (not reached)"
    print(f"  {m:<24} {status}")

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
summary_out_path = os.path.join(HERE, "phase6_combined_report_1001.csv")
summary_df.to_csv(summary_out_path, index=False)
print(f"\nwrote {summary_out_path}")


# =============================================================================
# Excel export (added 09/24/2026)
# =============================================================================
# The two CSVs above are kept as-is so anything that already reads them keeps
# working. This bundles the same data into one multi-sheet workbook, same
# pd.ExcelWriter convention as phase6_final_sensitivity_Tc_0917.py, plus a
# tidy long-form COP table (one row per method/ambient) which is easier to
# pivot or plot than the wide summary above.
#
#   cop_long        -- method, ambient, COP, converged, max_residual, attempts
#   summary_wide    -- exactly the phase6_combined_report_1001.csv layout
#   all_attempts    -- every attempt of the warmstart chase, every ambient
#   method_status   -- one row per method: attempts used, monotonic reached

cop_long_rows = []
for m in METHODS:
    for Tamb in AMBIENTS:
        r = final_state_by_method[m][Tamb]
        cop_long_rows.append({
            "method": m,
            "ambient_C": Tamb,
            "COP": r["cop"] if r.get("converged") else None,
            "converged": bool(r.get("converged", False)),
            "max_residual": r.get("max_residual"),
            "attempts_used": attempts_used[m],
            "reached_monotonic": bool(reached_monotonic[m]),
        })

status_rows = [{"method": m,
                "attempts_used": attempts_used[m],
                "reached_monotonic": bool(reached_monotonic[m]),
                "cap": MAX_ATTEMPTS,
                "n_converged_ambients": sum(
                    1 for Tamb in AMBIENTS
                    if final_state_by_method[m][Tamb].get("converged")),
                "n_ambients": len(AMBIENTS)}
               for m in METHODS]

xlsx_path = os.path.join(HERE, "phase6_final_1001.xlsx")
with pd.ExcelWriter(xlsx_path) as w:
    pd.DataFrame(cop_long_rows).to_excel(w, sheet_name="cop_long", index=False)
    summary_df.to_excel(w, sheet_name="summary_wide", index=False)
    attempts_df.to_excel(w, sheet_name="all_attempts", index=False)
    pd.DataFrame(status_rows).to_excel(w, sheet_name="method_status", index=False)
print(f"wrote {xlsx_path}")
