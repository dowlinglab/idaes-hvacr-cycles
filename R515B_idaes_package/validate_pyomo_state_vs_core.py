"""
Stage L (part 2) validation: native Pyomo P/h/s/g/Z and mu1/mu2 expressions
(`r515b_pyomo_eos.mixture_state_expr` / `mixture_chemical_potentials_expr`)
vs. the already-validated NumPy core (`r515b_helmholtz_core.py`'s
`mix_state` / `mix_entropy_direct` / `chemical_potentials_analytic`,
themselves validated exact vs. the oracle in Sections 5-7).

Method: identical pattern to `validate_pyomo_eos_vs_core.py` (Stage L part
1) -- build a Pyomo ConcreteModel with Var T, rho, x1 FIXED (not solved) at
the same 3 representative states, evaluate the new Pyomo expressions via
`value()`, and compare against the NumPy core's own output at the
identical (T, rho, x1). No solving happens -- this isolates whether the
Pyomo transcription (here: whether Pyomo's symbolic differentiation +
the chain-rule identities used to convert d/dT, d/drho, d/dx1 into
alpha_tau_mix/ar_del_mix/dar_dx1) is correct.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel, Var, value  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402
from r515b_pyomo_eos import mixture_state_expr, mixture_chemical_potentials_expr  # noqa: E402

TOL_REL = 1e-8


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

        state_e = mixture_state_expr(d1, d2, m.x1, m.T, m.rho, Tc1, Tc2, vc1, vc2)
        mu1_e, mu2_e = mixture_chemical_potentials_expr(d1, d2, m.x1, m.T, m.rho, Tc1, Tc2, vc1, vc2)

        P_pyomo = value(state_e["P_pa"])
        h_pyomo = value(state_e["h_jmol"])
        g_pyomo = value(state_e["g_jmol"])
        s_pyomo = value(state_e["s_jmolk"])  # raw, no Honeywell rebase (apply_entropy_offset=False)
        mu1_pyomo = value(mu1_e)
        mu2_pyomo = value(mu2_e)

        ref_state = core.mix_state(d1, d2, T_val, rho_val, x1_val)
        s_ref_rebased = core.mix_entropy_direct(d1, d2, T_val, rho_val, x1_val)
        s_ref_raw = s_ref_rebased - core.ENTROPY_REFERENCE_OFFSET_JMOLK
        mu1_ref, mu2_ref = core.chemical_potentials_analytic(d1, d2, T_val, rho_val, x1_val)

        checks = [
            ("P_pa", P_pyomo, ref_state.p_pa),
            ("h_jmol", h_pyomo, ref_state.h_jmol),
            ("g_jmol", g_pyomo, ref_state.g_jmol),
            ("s_jmolk (raw)", s_pyomo, s_ref_raw),
            ("mu1_jmol", mu1_pyomo, mu1_ref),
            ("mu2_jmol", mu2_pyomo, mu2_ref),
        ]
        state_pass = True
        print(f"[{label}] T={T_val} rho={rho_val}")
        for name, pyomo_val, ref_val in checks:
            e = rel(pyomo_val, ref_val)
            ok = e <= TOL_REL
            state_pass = state_pass and ok
            print(f"    {name}: pyomo={pyomo_val:.10g}  numpy_core={ref_val:.10g}  rel_err={e:.3g}  "
                  f"{'PASS' if ok else 'FAIL'}")
        all_pass = all_pass and state_pass
        print(f"    -> {label}: {'PASS' if state_pass else 'FAIL'}")

    print(f"\nOVERALL Stage L (part 2, Pyomo P/h/s/g/Z + mu1/mu2): {'PASS' if all_pass else 'FAIL'} "
          f"(tolerance {TOL_REL:.0e})")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
