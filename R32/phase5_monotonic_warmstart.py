"""
phase5_monotonic_warmstart.py -- copy of phase5_gcgp_warmstart.py, with a
different warmstart strategy. See BREADCRUMB_07-20.md for full context.
Original phase5_gcgp_warmstart.py (and phase4 files) untouched.

Difference from phase5_gcgp_warmstart.py: that file chained GCGP re-runs
WITHIN one ambient (Run 1 seeded by NIST, Run 2 by Run 1, ... up to 10,
all at the same Tamb). This file instead chains at the SWEEP level, for
BOTH GCGP and SPGP:

  Attempt 1: solve the method normally (no NIST seed, generic initial
             guess) at ALL FOUR ambients (10, 15, 20, 25), regardless of
             what the results look like.
  Check: is COP monotonically DECREASING as Tamb increases (10->15->20
         ->25), and did all four ambients converge? If yes, this
         method is done -- no "failure", nothing more to do.
  If no (non-monotonic trend, or any ambient failed to converge) --
         that attempt counts as "failed". Attempt 2: re-solve at EACH
         ambient, warmstarted from Attempt 1's OWN converged result at
         that SAME ambient (not from NIST). If Attempt 1 didn't
         converge at some particular ambient, that ambient falls back
         to a normal (non-warmstarted) solve for Attempt 2, since
         there's nothing to warmstart from there.
  Repeat: Attempt 3 seeded by Attempt 2, and so on, up to
          MAX_ATTEMPTS = 10 total attempts, stopping early the moment
          an attempt's COP-vs-Tamb sweep is fully converged AND
          monotonically decreasing.

Every attempt's full COP-vs-Tamb sweep is saved (not just the final
one), so the progression toward (or failure to reach) a monotonic trend
is visible, not just the end state.

Author: Shilpa Narasimhan
Support: Claude AI

Date created: 08/11/2026
"""
from pyomo.environ import value, Var
from vapor_compression_cubic_refstate import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression_plr import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
import pandas as pd

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3      # J/kg
MAX_ATTEMPTS = 10

STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]


def relax_enth_bounds(vc, hmax=H_MAX):
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        if v.parent_component().local_name == "enth_mass":
            if v.ub is not None and v.ub < hmax:
                v.setub(hmax)


def _cubic_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mol) / M_R32,
        "s": value(state.entr_mol) / M_R32,
        "x": value(state.phase_frac["Vap"]),
    }


def _helm_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mass),
        "s": value(state.entr_mass),
        "x": value(state.phase_frac["Vap"]),
    }


def warmstart_from_state(model, points):
    """Sets temperature/pressure initial values on a not-yet-solved
    cycle's stream outlets directly from a previous solve's converged
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
    return {"cop": cop, "converged": converged, "points": points}


def run_helm(Tamb):
    Tcond_sat = Tamb + 9
    vc = HelmCycle(FLUID, compressor_efficiency=0.9999, mode=HelmMode.IMPROVED_TPX)
    relax_enth_bounds(vc)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
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
    return {"cop": cop, "converged": converged, "points": points}


def is_monotonic_decreasing(cop_by_amb):
    """cop_by_amb: dict {Tamb: cop or None}. True only if EVERY ambient
    converged (no None values) AND COP strictly decreases as Tamb
    increases, matching the physically-expected trend."""
    cops = [cop_by_amb[Tamb] for Tamb in AMBIENTS]
    if any(c is None for c in cops):
        return False
    return all(cops[i] > cops[i + 1] for i in range(len(cops) - 1))


# --- Helmholtz and NIST: solved once per ambient, normally, as context
# (no iteration needed -- neither has shown this kind of trend problem). ---
helm_results = {}
nist_results = {}
for Tamb in AMBIENTS:
    try:
        helm_results[Tamb] = run_helm(Tamb)
    except Exception as e:
        helm_results[Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
    try:
        nist_results[Tamb] = run_cubic("NIST", Tamb)
    except Exception as e:
        nist_results[Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}

print(f"{'T_amb':>7}{'Helmholtz COP':>15}{'NIST COP':>10}")
for Tamb in AMBIENTS:
    h = helm_results[Tamb]
    n = nist_results[Tamb]
    h_cop = f"{h['cop']:.4f}" if h.get("converged") else "FAILED"
    n_cop = f"{n['cop']:.4f}" if n.get("converged") else "FAILED"
    print(f"{Tamb:>7}{h_cop:>15}{n_cop:>10}")

# --- GCGP and SPGP: the iterative monotonic-trend chase. ---
all_attempt_records = []  # one row per (method, attempt, Tamb) -- saved to CSV at the end
final_state_by_method = {}  # method -> {Tamb: result-or-None}, from whichever attempt stopped the loop

for method in ("GCGP", "SPGP"):
    print(f"\n{'='*90}\n  {method}: iterative monotonic-trend chase (cap = {MAX_ATTEMPTS} attempts)\n{'='*90}")

    prior_points_by_amb = {Tamb: None for Tamb in AMBIENTS}  # warmstart seed per ambient, None = no seed yet
    reached_monotonic = False

    for attempt in range(1, MAX_ATTEMPTS + 1):
        attempt_results = {}
        for Tamb in AMBIENTS:
            seed = prior_points_by_amb[Tamb]
            try:
                r = run_cubic(method, Tamb, warmstart_points=seed)
            except Exception as e:
                r = {"converged": False, "error": f"{type(e).__name__}: {e}", "points": None, "cop": None}
            attempt_results[Tamb] = r
            all_attempt_records.append({
                "method": method, "attempt": attempt, "Tamb_C": Tamb,
                "seeded_by": ("none" if seed is None else f"attempt {attempt - 1} own result"),
                "converged": r.get("converged", False),
                "COP": r.get("cop") if r.get("converged") else None,
                "error": r.get("error"),
            })

        cop_by_amb = {Tamb: (attempt_results[Tamb]["cop"] if attempt_results[Tamb].get("converged") else None)
                      for Tamb in AMBIENTS}
        cop_str = ", ".join(f"{Tamb}C={cop_by_amb[Tamb]:.4f}" if cop_by_amb[Tamb] is not None else f"{Tamb}C=FAILED"
                             for Tamb in AMBIENTS)
        monotonic = is_monotonic_decreasing(cop_by_amb)
        print(f"  Attempt {attempt}: {cop_str}  -- {'MONOTONIC (done)' if monotonic else 'not monotonic / has failures'}")

        # Next attempt's warmstart seed, per ambient: this attempt's own
        # converged points where available, otherwise fall back to no
        # seed (None) for that ambient next time.
        for Tamb in AMBIENTS:
            r = attempt_results[Tamb]
            prior_points_by_amb[Tamb] = r["points"] if r.get("converged") else None

        if monotonic:
            reached_monotonic = True
            final_state_by_method[method] = attempt_results
            print(f"  -> Reached a fully-converged, monotonically decreasing trend on attempt {attempt}.")
            break
    else:
        final_state_by_method[method] = attempt_results
        print(f"  -> Did NOT reach a monotonic trend within {MAX_ATTEMPTS} attempts.")

# --- Save every attempt, every ambient, both methods, to CSV. ---
attempts_df = pd.DataFrame(all_attempt_records)
attempts_df.to_csv("phase5_monotonic_warmstart_results.csv", index=False)
print(f"\nSaved {len(all_attempt_records)} rows to phase5_monotonic_warmstart_results.csv")
