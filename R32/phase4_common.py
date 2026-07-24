"""
phase4_common.py -- shared sweep logic for the Phase 4 per-method scripts
(phase4_nist_sweep.py, phase4_gcgp_sweep.py, phase4_spgp_sweep.py).

Each method gets its OWN dedicated driver script (per project direction:
"write a dedicated code for each method"), but all three need the exact
same sweep procedure and the exact same spec, so that logic lives here
once and each driver just calls run_sweep(method) and writes its own CSV.
A separate phase4_report.py reads all three CSVs and builds the combined
comparison table.

Sweep (matches phase_3a_helmholtz_cop.py exactly):
  T_amb in [10, 15, 20, 25] C, T_cond_sat = T_amb + 9 (condenser_approach),
  evap_sat_temperature = -29 C, superheating = subcooling = 0 (ideal cycle,
  Shridhar 2016), max_pressure_ratio = 10.

Helmholtz reference (Phase 3a, PHASE3_NOTES.md Section 1):
  T_amb  10    15    20    25
  COP   3.95  3.53  3.19  2.91

For each (method, ambient) cell we record, in order of what's actually
knowable:
  - cop, converged: the normal optimize_COP() outcome.
  - fallback: whether the debug_disable_arc_pressure_eq retry was needed.
  - max_residual: if NOT converged (but the model didn't crash), the
    largest constraint residual found by DiagnosticsToolbox -- a direct,
    quantitative answer to "how bad is it", not just a pass/fail flag.
  - error: if run_one() raised an exception (e.g. an InitializationError
    from vc.initialize() itself), the exception type/message. This is a
    deliberate, narrow exception to the "no try/except" SOP for exactly
    this reason -- a survey script's whole point is to characterize which
    cells work, so one bad cell must not kill the rest of the grid.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

AMBIENTS = [10, 15, 20, 25]  # deg C -- matches phase_3a_helmholtz_cop.py

# Phase 3a Helmholtz baseline, PHASE3_NOTES.md Section 1
HELMHOLTZ_REF = {10: 3.95, 15: 3.53, 20: 3.19, 25: 2.91}


def _max_constraint_residual(model):
    """Directly scan every active Constraint for |body - 0| (equality
    constraints only, which is everything in this model) and return the
    largest absolute residual. More robust than parsing
    DiagnosticsToolbox's printed report -- gives us a single float we can
    put straight into a CSV column."""
    from pyomo.environ import Constraint
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


def _compressor_diagnostics(vc):
    """Pull the same compressor-branch numbers
    compressor_fix_regression_check.py prints, for use INSIDE the actual
    sweep pipeline (added 2026-07-24 to debug phase4_gcgp_sweep.py
    directly, instead of re-deriving the same case in a separate one-off
    script -- see BREADCRUMB_07-20.md's pass 4/5/6 entries for why the
    isentropic-vs-real-outlet T/h relationship is the thing to watch)."""
    m = vc.model
    comp_out = m.fs.compressor.control_volume.properties_out[0]
    comp_isen = m.fs.compressor.properties_isentropic[0]
    e_out = m.fs.evaporator.control_volume.properties_out[0]
    c_out = m.fs.condenser.control_volume.properties_out[0]
    e_tsat = value(e_out.temperature_bubble["Vap", "Liq"])
    c_tsat = value(c_out.temperature_bubble["Vap", "Liq"])
    return {
        "comp_T_out": value(comp_out.temperature),
        "comp_h_out": value(comp_out.enth_mol),
        "comp_T_isen": value(comp_isen.temperature),
        "comp_h_isen": value(comp_isen.enth_mol),
        "evap_superheat": value(e_out.temperature) - e_tsat,
        "cond_subcool": c_tsat - value(c_out.temperature),
    }


def run_one(method, Tamb, extra_solver_options=None, diagnose=False):
    """Run the ideal-cycle spec for one (method, ambient) cell. Returns a
    dict with cop, converged, fallback, max_residual, error -- never
    raises (all exceptions are caught and recorded), so a caller can loop
    over the whole grid without any one cell killing the rest.

    extra_solver_options: if given, tried as a THIRD attempt (after the
    normal spec and the debug_disable_arc_pressure_eq fallback both fail),
    passed straight through to optimize_COP()'s solver_options override.
    Used by the SPGP sweep to make one more good-faith attempt at a
    tougher-converging parameter set before giving up on that cell --
    NOT applied to NIST/GCGP, which already converge with the standard
    two-attempt sequence.

    diagnose: if True, capture the compressor's real-outlet vs isentropic
    T/h (and evap/cond superheat/subcool) at the FINAL converged state
    and fold them into the returned dict -- lets a sweep expose the same
    "is this branch actually consistent" numbers used throughout this
    session's compressor debugging, for every ambient in one run, without
    duplicating the case-building logic in a separate script."""
    result = {"method": method, "ambient_C": Tamb, "T_cond_sat_C": Tamb + 9,
              "cop": None, "converged": False, "fallback": None,
              "max_residual": None, "third_attempt": None, "error": ""}
    try:
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

        if not converged and extra_solver_options:
            result["third_attempt"] = True
            cop, converged = vc.optimize_COP(
                verbose=False, initialize=True, optimize=False,
                solver_options=extra_solver_options,
            )

        result["cop"] = cop
        result["converged"] = converged
        result["fallback"] = fallback
        if not converged:
            result["max_residual"] = _max_constraint_residual(vc.model)
        if diagnose:
            result.update(_compressor_diagnostics(vc))
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"
    return result


def run_sweep(method, extra_solver_options=None, diagnose=False):
    """Run every ambient point in AMBIENTS for one method. Returns a list
    of per-cell result dicts (see run_one)."""
    rows = []
    for Tamb in AMBIENTS:
        print(f"\n{'='*70}\n  METHOD = {method}, T_amb = {Tamb} C\n{'='*70}")
        r = run_one(method, Tamb, extra_solver_options=extra_solver_options, diagnose=diagnose)
        r["helmholtz_ref_cop"] = HELMHOLTZ_REF[Tamb]
        r["pct_vs_helmholtz"] = (
            round(100.0 * (r["cop"] - HELMHOLTZ_REF[Tamb]) / HELMHOLTZ_REF[Tamb], 2)
            if r["converged"] else None
        )
        if r["error"]:
            print(f"  CRASHED: {r['error']}")
        else:
            cop_str = f"{r['cop']:.4f}" if r["converged"] else "--"
            resid_str = (f", max_residual={r['max_residual']:.4g}"
                         if r["max_residual"] is not None else "")
            print(f"  COP = {cop_str}, converged = {r['converged']}, "
                  f"fallback_used = {r['fallback']}{resid_str}")
            if diagnose and r["converged"]:
                print(f"    compressor REAL outlet:  T={r['comp_T_out']:.3f} K, h={r['comp_h_out']:.2f} J/mol")
                print(f"    compressor ISENTROPIC:   T={r['comp_T_isen']:.3f} K, h={r['comp_h_isen']:.2f} J/mol")
                print(f"    T_out - T_isen = {r['comp_T_out']-r['comp_T_isen']:+.3f} K, "
                      f"h_out - h_isen = {r['comp_h_out']-r['comp_h_isen']:+.3f} J/mol")
                print(f"    evap superheat = {r['evap_superheat']:+.3f} K, "
                      f"cond subcool = {r['cond_subcool']:+.3f} K")
        rows.append(r)
    return rows


def write_csv(rows, out_path):
    import csv
    fieldnames = ["method", "ambient_C", "T_cond_sat_C", "cop", "converged",
                  "fallback", "third_attempt", "max_residual", "helmholtz_ref_cop",
                  "pct_vs_helmholtz", "error"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"\nwrote {out_path}")
