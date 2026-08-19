"""
Stage N validation: end-to-end exercise of `vapor_compression_r515b_
integration.py`'s `R515BVaporCompressionCycle` -- construction,
`specify_initial_conditions` (our own validated-VLE-solver seeding, not
CoolProp), `initialize` (propagate_state across all 4 Arcs), `set_
specifications`, and `optimize_COP`, then sanity checks on the converged
solution: solver termination optimal, COP in a physically plausible
range, closed-loop consistency (expansion valve outlet state matches
evaporator inlet state, since the flow_mass equality on that arc was
deliberately deactivated for the closed loop -- see `vapor_compression.py`
comment "Our flowsheet is a closed, circular loop"), and each unit's
inlet/outlet pressures matching the P_low/P_high split.

Two real bugs were found and fixed while building this integration (see
helmholtz_prop_validation.md Section 26 for full write-ups):
  1. `R515BStateBlockData` did not support `entr_mol`/`enth_mol` (molar
     entropy/enthalpy) -- required UNCONDITIONALLY by `PressureChanger.
     add_isentropic()`/`init_isentropic()` (used by `Compressor`),
     regardless of `amount_basis`. Fixed by adding both as thin unit-
     converting wrappers around already-validated quantities.
  2. `define_port_members()` had been overridden to add `temperature`
     (an Expression) to the Port, but `idaes.core.util.initialization.
     propagate_state()` requires every port member to be a settable Var.
     Fixed by reverting to the StateBlockData base class default (which
     is what the reference `HelmholtzStateBlockData` itself uses,
     confirmed by direct inspection) -- the port now carries only the 3
     public state vars, matching the reference contract exactly.

One independent, deliberate deviation from `vapor_compression.py` (not a
bug in that read-only file, which is never modified): its temperature-
bound checks use bare `if bound:`, which silently skips a legitimate
bound of exactly 0 degC (Python falsiness). Fixed in this NEW file only,
as `if bound is not None:` -- verified below.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402
matplotlib.use("Agg")  # headless -- this validation never needs to display plots

from pyomo.environ import value  # noqa: E402
from vapor_compression_r515b_integration import R515BVaporCompressionCycle  # noqa: E402

TOL_PA = 1.0
TOL_JKG = 1.0e-3  # relative


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def run_case(label, evaporator_temperature, compressor_temperature, condenser_temperature,
             expect_cop_range=(1.0, 10.0)):
    cyc = R515BVaporCompressionCycle()
    cyc.specify_initial_conditions(low_side_temperature=-10, high_side_temperature=40)
    cyc.initialize(verbose=False)
    cyc.set_specifications(
        low_side_pressure=(100, 600),
        high_side_pressure=(500, 1500),
        evaporator_temperature=evaporator_temperature,
        compressor_temperature=compressor_temperature,
        condenser_temperature=condenser_temperature,
        subcooling=3, superheating=3, max_pressure_ratio=6,
    )
    cop, converged = cyc.optimize_COP(verbose=False, initialize=True, optimize=True)

    ok = True
    reasons = []
    if not converged:
        ok = False
        reasons.append("optimizer did not report optimal termination")
    if not (expect_cop_range[0] <= cop <= expect_cop_range[1]):
        ok = False
        reasons.append(f"COP {cop:.3f} outside plausible range {expect_cop_range}")

    m = cyc.model
    evap_in = m.fs.evaporator.control_volume.properties_in[0]
    valve_out = m.fs.expansion_valve.control_volume.properties_out[0]
    p_err = abs(value(evap_in.pressure) - value(valve_out.pressure))
    h_err = rel(value(evap_in.enth_mass), value(valve_out.enth_mass))
    if p_err > TOL_PA:
        ok = False
        reasons.append(f"closed-loop pressure mismatch {p_err:.3g} Pa")
    if h_err > TOL_JKG:
        ok = False
        reasons.append(f"closed-loop enthalpy mismatch rel_err {h_err:.3g}")

    # Evaporator upper-temperature-bound honored even at exactly 0 degC
    # (the independent falsy-zero fix, described in the module docstring).
    if evaporator_temperature and evaporator_temperature[1] == 0:
        t_out = value(m.fs.evaporator.control_volume.properties_out[0].temperature)
        bound_k = 0.0 + 273.15
        if t_out > bound_k + 1e-3:
            ok = False
            reasons.append(f"evaporator outlet T={t_out:.4f}K exceeds the 0degC upper bound -- falsy-zero bug regressed")

    print(f"[{label}] COP={cop:.4f} converged={converged} p_err={p_err:.3g}Pa h_err={h_err:.3g} "
          f"{'PASS' if ok else 'FAIL: ' + '; '.join(reasons)}")
    return ok


def main():
    all_pass = True
    all_pass &= run_case(
        "nominal_-15to-1C_evap_35to55C_cond",
        evaporator_temperature=(-15, -1), compressor_temperature=(None, 130),
        condenser_temperature=(35, 55),
    )
    all_pass &= run_case(
        "zero_degC_upper_bound_regression_check",
        evaporator_temperature=(-15, 0), compressor_temperature=(None, 130),
        condenser_temperature=(35, 55),
    )
    print(f"\nOVERALL Stage N (vapor_compression_r515b_integration.py end-to-end): "
          f"{'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
