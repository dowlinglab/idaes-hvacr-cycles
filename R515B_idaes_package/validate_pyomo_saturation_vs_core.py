"""
Stage L (part 3) validation: native Pyomo saturation-curve system
(`r515b_pyomo_eos.saturation_residuals_expr`) vs. the SciPy pseudo-pure
saturation solve (`r515b_helmholtz_core.solve_pseudopure_saturation_at_t`).

Method (deliberately different from the "fixed-Var, evaluate value()"
pattern used by Stage L parts 1-2's validators): THIS is the first Stage L
component that must actually be SOLVED as a genuine implicit NLP (3
unknowns T_sat/rho_l_sat/rho_v_sat, 3 equations), not just evaluated. Uses
IDAES's own `get_solver()` (IPOPT) via a small standalone ConcreteModel:
fixes pressure at several target values (taken from the already-validated
SciPy pseudo-pure solve's own converged output, so both sides are being
asked "what is the saturation state at THIS pressure"), frees T_sat/
rho_l_sat/rho_v_sat, and checks that IPOPT converges to the same state the
SciPy solver found -- i.e. this validates that the native-Pyomo
formulation is not just algebraically plausible but ACTUALLY SOLVABLE by
the real NLP solver this property package will run under.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from pyomo.environ import ConcreteModel, Var, Constraint, value  # noqa: E402
from idaes.core.solvers import get_solver  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402
from r515b_pyomo_eos import saturation_residuals_expr  # noqa: E402

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

    # NOTE: a first draft of this script used a coarse target list
    # (260/290/320/350/370K, i.e. up-to-30K seed jumps between points) and
    # found a spurious FAIL at 350K -- traced (see helmholtz_prop_
    # validation.md Section 22's matching entry) to the SciPy reference
    # solver itself converging to a non-physical (non-monotonic-in-T) root
    # when seeded from 30K away, NOT a Pyomo/IPOPT problem -- IPOPT's own
    # answer was actually CLOSER to the true continuation-consistent value
    # than the badly-seeded SciPy reference was. Fixed by walking a fine
    # (5K-step) continuation grid, exactly as `solve_bubble_at_t`/
    # `solve_pseudopure_saturation_at_t` are meant to be used in practice
    # (each temperature reseeded from the PREVIOUS converged solution, per
    # this project's own established convention -- see
    # `validate_pseudopure_saturation.py`'s working 15-point grid), so the
    # SciPy reference itself is trustworthy at every comparison point.
    t_targets = list(np.arange(255.0, 376.0, 5.0))
    solver = get_solver()

    all_pass = True
    n_checked = 0
    rl_seed, rv_seed = rho_l_seed0, rho_v_seed0
    for T_target in t_targets:
        sp = core.solve_pseudopure_saturation_at_t(d1, d2, x1_val, float(T_target), rl_seed, rv_seed)
        if sp["status"] != "CONVERGED":
            print(f"[T={T_target}] SciPy pseudo-pure solve DIVERGED -- skipping")
            all_pass = False
            continue
        rl_seed, rv_seed = sp["rho_l_molm3"], sp["rho_v_molm3"]
        P_target = sp["P_Pa"]
        n_checked += 1

        # Deliberately OFFSET initial guesses from the SciPy answer (not
        # handed the converged solution directly) so this is a genuine
        # solver-convergence check, not just a value() echo.
        m = ConcreteModel()
        m.x1 = Var(initialize=x1_val)
        m.x1.fix(x1_val)
        m.P = Var(initialize=P_target)
        m.P.fix(P_target)
        m.T_sat = Var(initialize=T_target * 1.05, bounds=(200.0, 420.0))
        m.rho_l_sat = Var(initialize=sp["rho_l_molm3"] * 1.15, bounds=(1.0, 2.0e4))
        m.rho_v_sat = Var(initialize=sp["rho_v_molm3"] * 0.7, bounds=(1.0e-6, 2.0e4))

        res_p_liq, res_p_vap, res_g = saturation_residuals_expr(
            d1, d2, m.x1, m.T_sat, m.rho_l_sat, m.rho_v_sat, m.P, Tc1, Tc2, vc1, vc2
        )
        m.eq_p_liq = Constraint(expr=res_p_liq == 0)
        m.eq_p_vap = Constraint(expr=res_p_vap == 0)
        m.eq_gibbs = Constraint(expr=res_g == 0)

        results = solver.solve(m, tee=False)
        solved_ok = str(results.solver.termination_condition) == "optimal"

        T_pyomo = value(m.T_sat)
        rl_pyomo = value(m.rho_l_sat)
        rv_pyomo = value(m.rho_v_sat)

        e_T = rel(T_pyomo, sp["T_K"])
        e_rl = rel(rl_pyomo, sp["rho_l_molm3"])
        e_rv = rel(rv_pyomo, sp["rho_v_molm3"])
        ok = solved_ok and max(e_T, e_rl, e_rv) <= TOL_REL
        all_pass = all_pass and ok
        print(f"[T_target={T_target}]  ipopt={results.solver.termination_condition}  "
              f"T_sat: pyomo={T_pyomo:.6f} scipy={sp['T_K']:.6f} rel_err={e_T:.3g} | "
              f"rho_l: pyomo={rl_pyomo:.4f} scipy={sp['rho_l_molm3']:.4f} rel_err={e_rl:.3g} | "
              f"rho_v: pyomo={rv_pyomo:.4f} scipy={sp['rho_v_molm3']:.4f} rel_err={e_rv:.3g}  "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\nChecked {n_checked}/{len(t_targets)} temperature points.")
    print(f"OVERALL Stage L (part 3, native-Pyomo saturation curve): {'PASS' if all_pass else 'FAIL'} "
          f"(tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
