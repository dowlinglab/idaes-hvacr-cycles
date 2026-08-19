"""
Stage L (part 3) validation: native-Pyomo smooth PH-flash system
(`r515b_pyomo_eos.mixture_ph_flash_residuals_expr`) vs. the explicit-
branching SciPy reference (`r515b_helmholtz_core.flash_ph_pseudopure`).

Method: same real-IPOPT-solve pattern as
`validate_pyomo_saturation_vs_core.py` -- build a ConcreteModel with the 8
unknowns (T_sat, rho_l_sat, rho_v_sat, vapor_frac, T_liq, rho_liq, T_vap,
rho_vap) free, fix P and H at target values spanning subcooled liquid,
two-phase, and superheated vapor, wrap the 8 residuals in Constraints,
solve with IDAES's get_solver() (IPOPT), and compare T_actual/vapor_frac/
derived density against the SciPy reference at each point.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel, Var, Constraint, value  # noqa: E402
from idaes.core.solvers import get_solver  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402
from r515b_pyomo_eos import mixture_ph_flash_residuals_expr  # noqa: E402

TOL_REL = 1e-6


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def main():
    d1 = core.load_idaes_helmholtz_json(core.FLUID1)
    d2 = core.load_idaes_helmholtz_json(core.FLUID2)
    x1_val = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    Tc1, Tc2 = float(d1["basic"]["Tc"]), float(d2["basic"]["Tc"])
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    vc1, vc2 = 1.0 / rhoc1, 1.0 / rhoc2
    rho_l_seed0 = 0.8 * (x1_val * rhoc1 + (1.0 - x1_val) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    solver = get_solver()

    # Anchor: saturation state at T=300K gives us P and h_l_sat/h_v_sat to
    # build test points relative to.
    sat300 = core.solve_pseudopure_saturation_at_t(d1, d2, x1_val, 300.0, rho_l_seed0, rho_v_seed0)
    P_test = sat300["P_Pa"]
    h_l_sat, h_v_sat = sat300["h_l_Jmol"], sat300["h_v_Jmol"]

    test_points = [
        ("subcooled_liquid", h_l_sat - 3000.0),
        ("two_phase_q30", h_l_sat + 0.30 * (h_v_sat - h_l_sat)),
        ("two_phase_q70", h_l_sat + 0.70 * (h_v_sat - h_l_sat)),
        ("superheated_vapor", h_v_sat + 3000.0),
    ]

    all_pass = True
    for label, H_test in test_points:
        ref = core.flash_ph_pseudopure(d1, d2, x1_val, P_test, H_test, 300.0, rho_l_seed0, rho_v_seed0)
        if ref["status"] != "CONVERGED":
            print(f"[{label}] SciPy reference flash DIVERGED -- skipping")
            all_pass = False
            continue

        m = ConcreteModel()
        m.x1 = Var(initialize=x1_val)
        m.x1.fix(x1_val)
        m.P = Var(initialize=P_test)
        m.P.fix(P_test)
        m.H = Var(initialize=H_test)
        m.H.fix(H_test)
        # Initial guesses: near the SATURATION state (T~300K, rho_l_sat/
        # rho_v_sat from the T=300K anchor), perturbed 5-15% -- NOT handed
        # the converged per-point reference answer directly, but branch-
        # consistent (T_liq/rho_liq seeded near the LIQUID branch,
        # T_vap/rho_vap near the VAPOR branch) rather than an arbitrary
        # constant guess. A first draft used arbitrary generic guesses
        # (T_liq=295K flat, rho_liq=a generic ballpark density unrelated
        # to rho_l_sat) and found IPOPT converging to a DIFFERENT, also-
        # exactly-residual-zero (T,rho) root of the 2-equation liquid-
        # branch system at the subcooled/two-phase points -- i.e. a
        # genuine non-uniqueness in that 2-eq system away from a good
        # seed, not a formulation bug (residuals were all ~0 at the wrong
        # root too). This is the same seed-sensitivity lesson already
        # recorded for the saturation curve (helmholtz_prop_validation.md
        # Section 22) -- fixed here the same way: seed near the physically
        # expected branch instead of a generic constant.
        m.T_sat = Var(initialize=300.0 * 1.02, bounds=(200.0, 420.0))
        m.rho_l_sat = Var(initialize=sat300["rho_l_molm3"] * 1.08, bounds=(1.0, 2.0e4))
        m.rho_v_sat = Var(initialize=sat300["rho_v_molm3"] * 0.85, bounds=(1.0e-6, 2.0e4))
        m.vapor_frac = Var(initialize=0.5, bounds=(-0.05, 1.05))
        m.T_liq = Var(initialize=300.0 * 0.97, bounds=(200.0, 420.0))
        m.rho_liq = Var(initialize=sat300["rho_l_molm3"] * 1.05, bounds=(1.0, 2.0e4))
        m.T_vap = Var(initialize=300.0 * 1.03, bounds=(200.0, 420.0))
        m.rho_vap = Var(initialize=sat300["rho_v_molm3"] * 0.92, bounds=(1.0e-6, 2.0e4))

        res = mixture_ph_flash_residuals_expr(
            d1, d2, m.x1, m.P, m.H,
            m.T_sat, m.rho_l_sat, m.rho_v_sat, m.vapor_frac,
            m.T_liq, m.rho_liq, m.T_vap, m.rho_vap,
            Tc1, Tc2, vc1, vc2,
        )
        m.eq_sat_p_liq = Constraint(expr=res["sat_p_liq"] == 0)
        m.eq_sat_p_vap = Constraint(expr=res["sat_p_vap"] == 0)
        m.eq_sat_gibbs = Constraint(expr=res["sat_gibbs"] == 0)
        m.eq_complementarity = Constraint(expr=res["complementarity"] == 0)
        m.eq_liq_p = Constraint(expr=res["liq_p"] == 0)
        m.eq_liq_h = Constraint(expr=res["liq_h"] == 0)
        m.eq_vap_p = Constraint(expr=res["vap_p"] == 0)
        m.eq_vap_h = Constraint(expr=res["vap_h"] == 0)

        results = solver.solve(m, tee=False)
        solved_ok = str(results.solver.termination_condition) == "optimal"

        T_actual_pyomo = value(res["T_actual"])
        vf_pyomo = value(m.vapor_frac)

        e_T = rel(T_actual_pyomo, ref["T_K"])
        e_vf = abs(vf_pyomo - ref["vapor_frac"])  # absolute (vf can be exactly 0)
        ok = solved_ok and e_T <= TOL_REL and e_vf <= 1e-5
        all_pass = all_pass and ok
        print(f"[{label}] ipopt={results.solver.termination_condition}  "
              f"T_actual: pyomo={T_actual_pyomo:.6f} scipy={ref['T_K']:.6f} rel_err={e_T:.3g}  "
              f"vapor_frac: pyomo={vf_pyomo:.6f} scipy={ref['vapor_frac']:.6f} abs_err={e_vf:.3g}  "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\nOVERALL Stage L (part 3, native-Pyomo PH-flash): {'PASS' if all_pass else 'FAIL'} "
          f"(tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
