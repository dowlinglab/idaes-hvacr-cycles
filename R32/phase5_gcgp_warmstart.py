"""
phase5_gcgp_warmstart.py -- copy of phase4_full_pfd_comparison_refstate.py,
extended to test whether GCGP's deviation from Helmholtz is genuine
prediction error (solution is unique regardless of starting point) or
partly a solver-path artifact (different starting points reach different
solutions). See BREADCRUMB_07-20.md for the full context this grew out
of. Original phase4 files are untouched.

New capability this file adds: warmstart_from_state() sets a NEW, not-
yet-solved cycle's temperature/pressure variables directly from a
PREVIOUS cycle's actual converged values (not just the two generic
low/high-side temperature anchors specify_initial_conditions uses).
This is valid because IMPROVED_TPX mode makes temperature and pressure
the real primary Pyomo variables being solved for, not derived
expressions.

Per ambient, for GCGP specifically:
  Run 1: GCGP warmstarted from NIST's converged state (NIST solved
         first, normally, no warmstart of its own).
  Run 2..10: GCGP warmstarted from the PREVIOUS GCGP run's own
         converged state -- each run seeded by the one before it,
         capped at 10 total runs.
The full sequence of 10 COP values is printed per ambient. If the
sequence is flat (no real trend), GCGP's solution is robust regardless
of starting point -- supporting "deviation = prediction error, not
solver path". If the sequence decreases (or otherwise moves), that's
evidence of real solution sensitivity to the starting point.

Helmholtz, NIST, and SPGP are still run each ambient exactly as in
phase4, for context -- only GCGP gets the warmstart-chain treatment.

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
MAX_GCGP_RERUNS = 10

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
    `points` dict (the same structure run_cubic/run_helm return). Valid
    under IMPROVED_TPX mode, where temperature/pressure are the actual
    primary state variables, not derived expressions."""
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
    """Same as phase4's run_cubic, with one addition: if warmstart_points
    is given (a previous solve's `points` dict), it overrides the
    generic initial guess AFTER vc.initialize() but BEFORE
    vc.set_specifications() -- so IDAES's own initialization still runs
    first, and our values are the last thing written before the real
    solve."""
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
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.ratioP[0]),
        "comp_isen_T": value(m.fs.compressor.properties_isentropic[0].temperature),
        "comp_isen_h": value(m.fs.compressor.properties_isentropic[0].enth_mol) / M_R32,
        "points": points,
    }


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
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.outlet.pressure[0]) / value(m.fs.compressor.inlet.pressure[0]),
        "comp_isen_T": value(m.fs.compressor.properties_isentropic[0].temperature),
        "comp_isen_h": value(m.fs.compressor.properties_isentropic[0].enth_mass),
        "points": points,
    }


METHODS = ["Helmholtz", "NIST", "GCGP", "SPGP"]
results = {m: {} for m in METHODS}
gcgp_warmstart_chains = {}  # Tamb -> list of COP values, one per rerun
warmstart_records = []      # one row per GCGP warmstart-chain attempt, saved to CSV at the end

for Tamb in AMBIENTS:
    print(f"\n{'='*70}\n  Running T_amb = {Tamb} C\n{'='*70}")

    # Helmholtz and NIST run exactly as in phase4 -- no warmstart chain.
    try:
        results["Helmholtz"][Tamb] = run_helm(Tamb)
        print(f"  Helmholtz: converged={results['Helmholtz'][Tamb]['converged']}, COP={results['Helmholtz'][Tamb]['cop']:.4f}")
    except Exception as e:
        results["Helmholtz"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
        print(f"  Helmholtz: FAILED -- {type(e).__name__}: {e}")

    try:
        results["NIST"][Tamb] = run_cubic("NIST", Tamb)
        print(f"  NIST: converged={results['NIST'][Tamb]['converged']}, COP={results['NIST'][Tamb]['cop']:.4f}")
    except Exception as e:
        results["NIST"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
        print(f"  NIST: FAILED -- {type(e).__name__}: {e}")

    # SPGP runs exactly as in phase4 -- expected to fail, included for context.
    try:
        results["SPGP"][Tamb] = run_cubic("SPGP", Tamb)
        print(f"  SPGP: converged={results['SPGP'][Tamb]['converged']}, COP={results['SPGP'][Tamb]['cop']:.4f}")
    except Exception as e:
        results["SPGP"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
        print(f"  SPGP: FAILED -- {type(e).__name__}: {e}")

    # GCGP: the actual warmstart-chain experiment.
    cop_chain = []
    nist_result = results["NIST"].get(Tamb, {})
    if not nist_result.get("converged"):
        print("  GCGP warmstart chain: SKIPPED -- NIST did not converge at this ambient, nothing to seed from.")
        results["GCGP"][Tamb] = {"converged": False, "error": "NIST did not converge; no warmstart seed"}
        gcgp_warmstart_chains[Tamb] = cop_chain
        continue

    warmstart_points = nist_result["points"]
    last_result = None
    for run_idx in range(1, MAX_GCGP_RERUNS + 1):
        try:
            r = run_cubic("GCGP", Tamb, warmstart_points=warmstart_points)
        except Exception as e:
            print(f"  GCGP run {run_idx}: FAILED -- {type(e).__name__}: {e}")
            warmstart_records.append({
                "Tamb_C": Tamb, "run_idx": run_idx, "seeded_by": ("NIST" if run_idx == 1 else f"GCGP run {run_idx-1}"),
                "converged": False, "COP": None, "error": f"{type(e).__name__}: {e}",
            })
            break
        seed_desc = "NIST" if run_idx == 1 else f"GCGP run {run_idx - 1}"
        print(f"  GCGP run {run_idx} (seeded by {seed_desc}): converged={r['converged']}, COP={r['cop']:.4f}")
        # Save this run's result regardless of outcome -- one row per
        # attempt, including the full per-stream state, not just COP.
        row = {
            "Tamb_C": Tamb, "run_idx": run_idx, "seeded_by": seed_desc,
            "converged": r["converged"], "COP": r["cop"] if r["converged"] else None,
            "error": None,
        }
        for sp in STREAM_POINTS:
            pt = r["points"][sp]
            row[f"{sp}_T"] = pt["T"]
            row[f"{sp}_P"] = pt["P"]
            row[f"{sp}_h"] = pt["h"]
            row[f"{sp}_s"] = pt["s"]
            row[f"{sp}_x"] = pt["x"]
        warmstart_records.append(row)
        if not r["converged"]:
            break
        cop_chain.append(r["cop"])
        last_result = r
        warmstart_points = r["points"]  # next run seeded by THIS run's own result

    gcgp_warmstart_chains[Tamb] = cop_chain
    results["GCGP"][Tamb] = last_result if last_result is not None else {"converged": False, "error": "no successful run in chain"}

# --- Summary: COP, pressure ratio, isentropic reference block ---
print(f"\n{'='*100}\n  Summary: COP, pressure ratio, isentropic reference (T,h)\n{'='*100}")
print(f"{'T_amb':>7}{'method':>12}{'COP':>9}{'ratioP':>9}{'T_isen(K)':>11}{'h_isen(kJ/kg)':>15}")
for Tamb in AMBIENTS:
    for method in METHODS:
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            err = r.get("error", "not converged")
            print(f"{Tamb:>7}{method:>12}  FAILED: {err}")
            continue
        print(f"{Tamb:>7}{method:>12}{r['cop']:>9.4f}{r['ratioP']:>9.4f}"
              f"{r['comp_isen_T']:>11.3f}{r['comp_isen_h']/1e3:>15.2f}")

# --- GCGP warmstart-chain COP sequences, one row per ambient ---
print(f"\n{'='*100}\n  GCGP warmstart-chain COP sequence per ambient (Run 1 seeded by NIST,\n"
      f"  each later run seeded by the previous GCGP run; capped at {MAX_GCGP_RERUNS})\n{'='*100}")
for Tamb in AMBIENTS:
    chain = gcgp_warmstart_chains.get(Tamb, [])
    if not chain:
        print(f"  T_amb={Tamb}C: no successful runs")
        continue
    chain_str = " -> ".join(f"{c:.4f}" for c in chain)
    print(f"  T_amb={Tamb}C ({len(chain)} runs): {chain_str}")

# --- Save every GCGP warmstart-chain attempt to CSV -- one row per run,
# per ambient, including the full per-stream state, not just COP, so
# nothing has to be re-run to inspect what each step in the chain
# actually converged to. ---
warmstart_df = pd.DataFrame(warmstart_records)
warmstart_df.to_csv("phase5_gcgp_warmstart_chain.csv", index=False)
print(f"\nSaved {len(warmstart_records)} GCGP warmstart-chain rows to phase5_gcgp_warmstart_chain.csv")

# --- Full PFD table: one block per (T_amb, method), all 4 stream points ---
print(f"\n{'='*100}\n  Full PFD state table: T(K), P(bar), h(kJ/kg), s(kJ/kg-K), quality\n{'='*100}")
for Tamb in AMBIENTS:
    for method in METHODS:
        r = results[method].get(Tamb, {})
        print(f"\n--- T_amb={Tamb}C, {method} ---")
        if not r.get("converged"):
            print(f"  NOT CONVERGED: {r.get('error', 'unknown')}")
            continue
        print(f"  {'stream':>10}{'T(K)':>10}{'P(bar)':>10}{'h(kJ/kg)':>11}{'s(kJ/kg-K)':>12}{'quality':>9}")
        for sp in STREAM_POINTS:
            pt = r["points"][sp]
            print(f"  {sp:>10}{pt['T']:>10.3f}{pt['P']/1e5:>10.4f}{pt['h']/1e3:>11.3f}{pt['s']/1e3:>12.4f}{pt['x']:>9.4f}")
