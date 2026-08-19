"""
Validates the 2026-08-18 additive extension to `mixture_ph_flash_
residuals_expr` (s_liq_jmolk/s_vap_jmolk/S_actual_jmolk) and the new
`R515BStateBlockData.entr_mass`/`temperature_sat` Expressions (added for
Stage N/vapor_compression.py parity) against the SciPy reference:
`r515b_helmholtz_core.mix_entropy_direct` for entropy, and `flash_ph_
pseudopure`'s own `T_sat_K` for saturation temperature. Same real-IPOPT-
solve pattern as `validate_state_block_construction.py`, at the same 4
representative test points.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel, value  # noqa: E402
from idaes.core import FlowsheetBlock  # noqa: E402

import r515b_helmholtz_core as core  # noqa: E402
from r515b_property_package import R515BParameterBlock  # noqa: E402

TOL_REL = 1e-6


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def build_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = R515BParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
    return m


def main():
    all_pass = True

    d1 = core.load_idaes_helmholtz_json(core.FLUID1)
    d2 = core.load_idaes_helmholtz_json(core.FLUID2)
    x1_val = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1_val * rhoc1 + (1.0 - x1_val) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0
    mw_mix = x1_val * mw1 + (1.0 - x1_val) * mw2

    sat300 = core.solve_pseudopure_saturation_at_t(d1, d2, x1_val, 300.0, rho_l_seed0, rho_v_seed0)
    P_test = sat300["P_Pa"]
    h_l_sat, h_v_sat = sat300["h_l_Jmol"], sat300["h_v_Jmol"]

    test_points = [
        ("subcooled_liquid", h_l_sat - 3000.0),
        ("two_phase_q30", h_l_sat + 0.30 * (h_v_sat - h_l_sat)),
        ("two_phase_q70", h_l_sat + 0.70 * (h_v_sat - h_l_sat)),
        ("superheated_vapor", h_v_sat + 3000.0),
    ]

    for label, H_molar_test in test_points:
        ref = core.flash_ph_pseudopure(d1, d2, x1_val, P_test, H_molar_test, 300.0, rho_l_seed0, rho_v_seed0)
        if ref["status"] != "CONVERGED":
            print(f"[{label}] SciPy reference flash DIVERGED -- skipping")
            all_pass = False
            continue

        # Reference entropy: direct evaluation at the reference's own
        # converged (T,rho) state(s) -- lever-rule blended in the two-phase
        # region, exactly matching the pyomo side's S_actual_jmolk blend.
        if ref["region"] == "two_phase":
            s_l_ref = core.mix_entropy_direct(d1, d2, ref["T_sat_K"], ref["rho_l_sat_molm3"], x1_val)
            s_v_ref = core.mix_entropy_direct(d1, d2, ref["T_sat_K"], ref["rho_v_sat_molm3"], x1_val)
            s_ref_jmolk = (1.0 - ref["vapor_frac"]) * s_l_ref + ref["vapor_frac"] * s_v_ref
        else:
            s_ref_jmolk = core.mix_entropy_direct(d1, d2, ref["T_K"], ref["rho_molm3"], x1_val)
        entr_mass_ref = s_ref_jmolk / mw_mix
        t_sat_ref = ref["T_sat_K"]

        h_mass_test = H_molar_test / mw_mix

        m2 = build_flowsheet()
        blk2 = m2.fs.state[0]
        blk2.flow_mass.fix(1.0)
        blk2.pressure.fix(P_test)
        blk2.enth_mass.fix(h_mass_test)
        m2.fs.state.initialize(outlvl=0)

        entr_mass_pyomo = value(blk2.entr_mass)
        t_sat_pyomo = value(blk2.temperature_sat)

        e_s = rel(entr_mass_pyomo, entr_mass_ref)
        e_tsat = rel(t_sat_pyomo, t_sat_ref)
        ok = e_s <= TOL_REL and e_tsat <= TOL_REL
        all_pass = all_pass and ok
        print(f"[{label}] entr_mass: pyomo={entr_mass_pyomo:.6f} scipy={entr_mass_ref:.6f} rel_err={e_s:.3g}  "
              f"temperature_sat: pyomo={t_sat_pyomo:.6f} scipy={t_sat_ref:.6f} rel_err={e_tsat:.3g}  "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\nOVERALL entropy/temperature_sat validation: {'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
