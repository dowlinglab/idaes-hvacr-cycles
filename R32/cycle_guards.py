"""
cycle_guards.py -- physical checks that a solved vapour-compression cycle has
to pass before its COP is reported. 10/08/2026.

WHY THIS EXISTS
    On 10/08/2026 the flowsheet was found to be producing converged,
    monotonic, plausible-looking COPs from a loop that does not conserve
    energy. The arcs enforce flow_mol, mole_frac_comp, pressure and
    temperature, but NOT enthalpy or phase fraction. For a pure fluid at
    saturation, (T, P) are not independent -- they lie on the saturation
    line -- so the vapour fraction is left free and two blocks can satisfy
    every arc equality while being in different thermodynamic states.

    Measured on the R-32 poster configuration (NIST, optimize=False):
    the stream leaving the expansion valve arrived at the evaporator inlet
    with an enthalpy 26-30 % lower, and the loop energy balance closed to
    -1280 W at T_amb = 15 C and +117 W at 25 C -- opposite signs.

    Every one of those runs reported converged = True. Nothing in the
    workflow asked whether the answer was physically possible.

WHAT IT CHECKS
    1. Carnot ceiling          COP <= T_evap / (T_cond - T_evap)
    2. Loop energy closure     Q_evap + W_comp + Q_cond == 0
    3. Arc enthalpy continuity h_out of one unit == h_in of the next
    4. Mass continuity         same flow_mol everywhere
    5. Spec compliance         superheat / subcool match what was requested.
                               Pass the SAME values you gave
                               set_specifications via superheat_K / subcool_K.
                               NOTE: set_specifications DEACTIVATES these
                               constraints when passed 0, so 0 does not mean
                               "saturated" -- it means "unconstrained". Any
                               non-zero value above 0.1 does activate them.

HOW TO USE IT
    from cycle_guards import check_cycle

    cop, converged = vc.optimize_COP(verbose=False, initialize=True,
                                     optimize=False)
    rpt = check_cycle(vc.model, cop, T_evap_C=-29.0, T_cond_C=Tamb+9.0,
                      fluid="R32", strict=False)
    if not rpt["ok"]:
        print(rpt["summary"])

    strict=True raises AssertionError instead of returning. Use strict=True
    in anything whose output gets reported; strict=False when you are
    deliberately studying the broken behaviour.

Works with either fluid -- it reads enth_mol and flow_mol, never a component
name. `fluid` is only used for the saturation lookups in check 5 and may be
left as None to skip that check.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""

from pyomo.environ import value

# Connections, in flow order. (upstream unit, downstream unit)
ARCS = [
    ("evaporator", "compressor"),
    ("compressor", "condenser"),
    ("condenser", "expansion_valve"),
    ("expansion_valve", "evaporator"),
]

UNITS = ["evaporator", "compressor", "condenser", "expansion_valve"]


def _cv(model, name):
    return getattr(model.fs, name).control_volume


def check_cycle(model, cop, T_evap_C, T_cond_C, fluid=None,
                superheat_K=0.0, subcool_K=0.0,
                tol_W=1.0, tol_arc_J=1.0, tol_spec_K=0.05, strict=True):
    """Run every physical check on a solved cycle.

    Returns a dict with keys: ok, failures (list of strings), summary (str),
    and the measured quantities. Raises AssertionError on any failure if
    strict is True.
    """
    fail = []
    out = {}

    Te = T_evap_C + 273.15
    Tc = T_cond_C + 273.15

    # ---- 1. Carnot ceiling ---------------------------------------------
    carnot = Te / (Tc - Te)
    out["cop"] = cop
    out["carnot"] = carnot
    if cop is not None and cop > carnot:
        fail.append(f"COP {cop:.4f} EXCEEDS the Carnot ceiling {carnot:.4f} "
                    f"for T_evap={T_evap_C} C, T_cond={T_cond_C} C "
                    f"(ratio {cop/carnot:.2f})")

    # ---- 2. Loop energy closure ----------------------------------------
    Q = value(model.fs.evaporator.heat_duty[0])
    W = value(model.fs.compressor.work_mechanical[0])
    K = value(model.fs.condenser.heat_duty[0])
    close = Q + W + K
    out.update(Q_evap=Q, W_comp=W, Q_cond=K, loop_closure=close)
    if abs(close) > tol_W:
        fail.append(f"loop does not close: Q_evap + W_comp + Q_cond = "
                    f"{close:.2f} W  ({100*abs(close)/abs(Q):.1f} % of Q_evap)")

    # ---- 3. Arc enthalpy continuity ------------------------------------
    # THE one that was broken. A stream leaving a unit must arrive at the
    # next unit with the same specific enthalpy -- they are joined by a pipe.
    out["arc_breaks"] = {}
    for up, down in ARCS:
        h_out = value(_cv(model, up).properties_out[0].enth_mol)
        h_in = value(_cv(model, down).properties_in[0].enth_mol)
        d = h_in - h_out
        out["arc_breaks"][f"{up}->{down}"] = d
        if abs(d) > tol_arc_J:
            pct = 100 * abs(d) / abs(Q / value(_cv(model, up).properties_out[0].flow_mol)) \
                  if Q else float("nan")
            fail.append(f"arc {up} -> {down}: enthalpy jumps {d:.1f} J/mol "
                        f"({pct:.1f} % of the evaporator duty per mol)")

    # ---- 4. Mass continuity --------------------------------------------
    flows = {}
    for u in UNITS:
        flows[u] = (value(_cv(model, u).properties_in[0].flow_mol),
                    value(_cv(model, u).properties_out[0].flow_mol))
    out["flows"] = flows
    allf = [x for pair in flows.values() for x in pair]
    if max(allf) - min(allf) > 1e-4:
        fail.append(f"mass flow is not constant around the loop: "
                    f"{min(allf):.6f} to {max(allf):.6f} mol/s")

    # ---- 5. Spec compliance (optional, needs CoolProp + fluid) ---------
    if fluid is not None:
        try:
            from CoolProp.CoolProp import PropsSI
            eo = _cv(model, "evaporator").properties_out[0]
            ko = _cv(model, "condenser").properties_out[0]
            sh = value(eo.temperature) - PropsSI("T", "P", value(eo.pressure),
                                                 "Q", 1, fluid)
            sc = PropsSI("T", "P", value(ko.pressure), "Q", 0, fluid) \
                - value(ko.temperature)
            out.update(superheat_K=sh, subcool_K=sc)
            # NOTE the constraints in set_specifications are INEQUALITIES
            # (>= for superheat, <= for subcool), so the solved value is only
            # required to be on one side of the request. Reported either way.
            if abs(sh - superheat_K) > tol_spec_K:
                fail.append(f"evaporator outlet superheat is {sh:+.2f} K, "
                            f"not the {superheat_K:g} K requested")
            if abs(sc - subcool_K) > tol_spec_K:
                fail.append(f"condenser outlet subcool is {sc:+.2f} K, "
                            f"not the {subcool_K:g} K requested"
                            + ("  (NEGATIVE: outlet is ABOVE saturation while "
                               "flagged as liquid)" if sc < -tol_spec_K else ""))
        except Exception as e:                      # pragma: no cover
            out["spec_check_error"] = f"{type(e).__name__}: {e}"

    # ---- report ---------------------------------------------------------
    out["ok"] = not fail
    out["failures"] = fail
    lines = [f"  COP {cop:.4f}   Carnot {carnot:.4f}   "
             f"loop closure {close:+.2f} W" if cop is not None else
             f"  (no COP)   Carnot {carnot:.4f}   loop closure {close:+.2f} W"]
    if fail:
        lines.append(f"  {len(fail)} PHYSICAL CHECK(S) FAILED:")
        lines += [f"    - {f}" for f in fail]
    else:
        lines.append("  all physical checks passed")
    out["summary"] = "\n".join(lines)

    if strict and fail:
        raise AssertionError("cycle failed physical checks:\n" + out["summary"])
    return out


if __name__ == "__main__":
    # Demonstration on the R-32 poster configuration, which is known to fail.
    import logging, warnings, io, contextlib
    logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
    from vapor_compression_cubic_refstate_1001 import (
        SimpleVaporCompressionCycle as C, Mode)

    T_EVAP, APPROACH, TAMB = -29.0, 9.0, 20
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        vc = C("R32", compressor_efficiency=0.9999,
               mode=Mode.IMPROVED_TPX, method="NIST")
        vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                      high_side_temperature=TAMB + APPROACH)
        vc.initialize(verbose=False)
        vc.set_specifications(ambient_temperature=TAMB,
                              condenser_approach=APPROACH,
                              evap_sat_temperature=T_EVAP, superheating=0,
                              subcooling=0, max_pressure_ratio=10)
        cop, conv = vc.optimize_COP(verbose=False, initialize=True,
                                    optimize=False)

    print(f"solver reported: COP = {cop:.4f}, converged = {conv}\n")
    rpt = check_cycle(vc.model, cop, T_evap_C=T_EVAP,
                      T_cond_C=TAMB + APPROACH, fluid="R32", strict=False)
    print(rpt["summary"])
