"""
Stage L (part 3, final): validate that `r515b_property_package.py`'s
`R515BParameterBlock`/`R515BStateBlock` actually CONSTRUCT correctly inside
a real IDAES FlowsheetBlock, have the expected degrees of freedom before/
after fixing state vars, and that `.initialize()` (SciPy-seeded, per the
rule-44 carve-out) followed by a real IPOPT solve of the embedded 8-equation
smooth PH-flash system reproduces `r515b_helmholtz_core.flash_ph_pseudopure`
(the validated SciPy reference) at representative subcooled/two-phase/
superheated test points.

This directly tests the "plumbing" added in r515b_property_package.py -- the
underlying thermodynamic math (mixture_ph_flash_residuals_expr) was already
independently validated in validate_pyomo_flash_vs_core.py.

Per the 2026-08-18 clarification (R-515B ancillary correlations do NOT
exist / are NOT to be silently invented -- Option 1 "preserve the validated
VLE equations algebraically in the IDAES formulation" is the only approach
in use): this script exercises exactly that pathway. There is no fitted
P_sat=f(T)/rho_sat=g(T) anywhere in the property package; the StateBlockData
solves the full 8-equation binary-Helmholtz-derived VLE system as ordinary
Pyomo Constraints, with SciPy used only as an initialize()-time seed
generator (rule-44 carve-out), never inside the active NLP's constraint
bodies.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel, value  # noqa: E402
from idaes.core import FlowsheetBlock  # noqa: E402
from idaes.core.util.model_statistics import degrees_of_freedom  # noqa: E402

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

    # ---- Step 1: construction + DOF check ----
    # NOTE (self-caught test-script bug, fixed here -- see helmholtz_prop_
    # validation.md Section 22 2026-08-18 entry "DOF-before expectation was
    # wrong, not the package"): a first draft of this script asserted
    # dof_before == 3, reasoning naively from total_vars(11) -
    # total_constraints(8). That is NOT how IDAES's degrees_of_freedom()
    # actually computes -- it is number_unfixed_variables_in_activated_
    # equalities(block) - number_activated_equalities(block), i.e. it only
    # counts variables that actually APPEAR in an activated equality. Our
    # `flow_mass` never appears in any of the 8 embedded VLE-flash
    # constraints (correctly -- it is a pure extensive throughput variable,
    # structurally independent of the intensive T/P/h/rho system), so it is
    # excluded from the count. Direct check against the reference pure-
    # fluid `HelmholtzParameterBlock` (StateVars.PH, AmountBasis.MASS)
    # confirms this is the SAME convention already used by the package we
    # are matching: that reference block reports degrees_of_freedom()==0
    # both before and after fixing state vars (it has ZERO Pyomo-visible
    # constraints at all -- temperature/vapor_frac come from black-box
    # external functions there, per the general_helmholtz research already
    # recorded). Our block has 10 vars-in-constraints and 8 constraints
    # (pressure, enth_mass are both counted; flow_mass is not), so the
    # correct expectation is dof_before == 2, dof_after == 0 -- both are
    # exactly what fixing the 3 public state vars must deliver. This was a
    # test-expectation bug, not a package bug; documented per spec rule
    # 74/89 (never hide a failure, including self-inflicted test bugs).
    m = build_flowsheet()
    blk = m.fs.state[0]
    dof_before = degrees_of_freedom(blk)
    print(f"[construction] DOF before fixing state vars: {dof_before} (expect 2 -- see note in source)")
    ok_dof_before = dof_before == 2
    all_pass = all_pass and ok_dof_before

    blk.flow_mass.fix(1.0)
    blk.pressure.fix(5.0e5)
    blk.enth_mass.fix(3.0e5)
    dof_after = degrees_of_freedom(blk)
    print(f"[construction] DOF after fixing flow_mass/pressure/enth_mass: {dof_after} (expect 0)")
    ok_dof_after = dof_after == 0
    all_pass = all_pass and ok_dof_after

    # ---- Step 2: initialize() + solve, compared against SciPy reference ----
    d1 = core.load_idaes_helmholtz_json(core.FLUID1)
    d2 = core.load_idaes_helmholtz_json(core.FLUID2)
    x1_val = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1_val * rhoc1 + (1.0 - x1_val) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0
    mw_mix = value(m.fs.properties.mw)

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

        h_mass_test = H_molar_test / mw_mix

        m2 = build_flowsheet()
        blk2 = m2.fs.state[0]
        blk2.flow_mass.fix(1.0)
        blk2.pressure.fix(P_test)
        blk2.enth_mass.fix(h_mass_test)

        m2.fs.state.initialize(outlvl=0)

        T_pyomo = value(blk2.temperature)
        vf_pyomo = value(blk2.vapor_frac)
        dens_mass_pyomo = value(blk2.dens_mass)

        rho_ref = ref["rho_molm3"] if ref["region"] != "two_phase" else (
            1.0 / ((1.0 - ref["vapor_frac"]) / ref["rho_l_sat_molm3"] + ref["vapor_frac"] / ref["rho_v_sat_molm3"])
        )
        dens_mass_ref = rho_ref * mw_mix

        e_T = rel(T_pyomo, ref["T_K"])
        e_vf = abs(vf_pyomo - ref["vapor_frac"])
        e_dens = rel(dens_mass_pyomo, dens_mass_ref)
        ok = e_T <= TOL_REL and e_vf <= 1e-5 and e_dens <= TOL_REL
        all_pass = all_pass and ok
        print(f"[{label}] T: pyomo={T_pyomo:.6f} scipy={ref['T_K']:.6f} rel_err={e_T:.3g}  "
              f"vf: pyomo={vf_pyomo:.6f} scipy={ref['vapor_frac']:.6f} abs_err={e_vf:.3g}  "
              f"dens_mass: pyomo={dens_mass_pyomo:.6f} scipy={dens_mass_ref:.6f} rel_err={e_dens:.3g}  "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\nOVERALL Stage L (part 3, StateBlockData construction+DOF+initialize/solve): "
          f"{'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
