"""
Stage L (part 1) validation: native Pyomo EOS-kernel expressions
(`r515b_pyomo_eos.py`) vs. the already-validated NumPy core
(`r515b_helmholtz_core.py`, itself validated exact vs. the oracle in
Stages D-K).

Method: build a tiny Pyomo ConcreteModel with Var T, rho, x1 FIXED (not
solved) at the same representative states used throughout this task,
evaluate `mixture_alpha_and_derivs_expr`'s returned `alpha_mix` expression
via `pyomo.environ.value()`, and compare against
`r515b_helmholtz_core._mix_alpha_and_derivs`'s own `alpha` output at the
identical (T, rho, x1). This isolates whether the Pyomo transcription
itself is correct -- no solving happens here, just expression evaluation.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel, Var, value  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402
from r515b_pyomo_eos import mixture_alpha_and_derivs_expr  # noqa: E402

TOL_REL = 1e-10


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

    states = [
        ("liquid_direct", 300.0, 10773.988537),
        ("vapor_direct", 300.0, 124.343493),
        ("supercritical_direct", 389.5, 3000.0),
    ]

    all_pass = True
    for label, T_val, rho_val in states:
        m = ConcreteModel()
        m.T = Var(initialize=T_val)
        m.T.fix(T_val)
        m.rho = Var(initialize=rho_val)
        m.rho.fix(rho_val)
        m.x1 = Var(initialize=x1_val, bounds=(1e-9, 1 - 1e-9))
        m.x1.fix(x1_val)

        tau_e, delta_e, alpha_mix_e, Tred_e, vred_e = mixture_alpha_and_derivs_expr(
            d1, d2, m.x1, m.T, m.rho, Tc1, Tc2, vc1, vc2
        )
        alpha_pyomo = value(alpha_mix_e)
        tau_pyomo = value(tau_e)
        delta_pyomo = value(delta_e)

        alpha_ref, alpha_tau_ref, ar_del_ref = core._mix_alpha_and_derivs(d1, d2, T_val, rho_val, x1_val)

        e_alpha = rel(alpha_pyomo, alpha_ref)
        ok = e_alpha <= TOL_REL
        all_pass = all_pass and ok
        print(f"[{label}] T={T_val} rho={rho_val}  status={'PASS' if ok else 'FAIL'}")
        print(f"    alpha_mix: pyomo={alpha_pyomo:.12g}  numpy_core={alpha_ref:.12g}  rel_err={e_alpha:.3g}")
        print(f"    tau: pyomo={tau_pyomo:.10g}  delta: pyomo={delta_pyomo:.10g}")

    print(f"\nOVERALL Stage L (part 1, Pyomo EOS kernel): {'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
