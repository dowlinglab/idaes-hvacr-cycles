"""
Stage H (partial) validation: confirm that `linear_model_codex.py`'s
`compute_table1_properties` fugacity/chemical-potential outputs (f1, f2, mu1,
mu2) match the oracle's own `chemical_potentials_analytic` (mu1, mu2) at
representative single-phase states.

Background: the oracle's `chemical_potentials_analytic` uses a LOCAL
re-implementation `_bell2023_reducing_derivs_binary_local` (defined inside
mixture_fully_validated.py itself, NOT imported) of the same composition-
derivative formulas that `linear_model_codex.py`'s own
`_bell2023_reducing_derivs_binary` implements -- confirmed by reading both
function bodies side by side: same algebraic formula, same variable-by-
variable structure, just independently written out. This is exactly the
kind of "mathematically equivalent, not automatically numerically
identical" situation the MASTER TASK spec warns about (rule 10), so unlike
Sections 2-4 (genuinely shared code paths), this one is NOT guaranteed
identical by construction and must be checked numerically -- which is what
this script does.

Everything else feeding chemical_potentials_analytic (mixture_alpha0_
alphar_derivs, alphar_idaes_with_derivs, bell2023_departure_alphar/base,
bell2023_Tred_vred) IS a shared code path (Sections 2-4), so any mismatch
found here would isolate specifically to the reducing-function composition-
derivative piece.

Read-only reference use of the oracle (spec rule 9 carve-out) for
comparison only.
"""

import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))

from mixture_fully_validated import (  # noqa: E402
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    chemical_potentials_analytic,
    solve_bubble_at_t,
    solve_dew_at_t,
)
from linear_model_codex import compute_table1_properties  # noqa: E402

FLUID1, FLUID2 = "r1234ze", "r227ea"
W1 = 0.911
TOL_REL = 1e-6  # composition-derivative piece is an independent re-implementation,
                 # not a shared code path -- use a looser (still tight) check than
                 # the 1e-8 tier-1 floor reserved for genuinely shared-code-path comparisons.


def compare_state(label, T, rho_mol, x1, mw1, mw2):
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mu1_ref, mu2_ref = chemical_potentials_analytic(d1, d2, T, rho_mol, x1)

    rho_mass = rho_mol * (x1 * mw1 + (1.0 - x1) * mw2)
    new = compute_table1_properties(FLUID1, FLUID2, W1, T, rho_mass, fd_rel=1e-6)
    mu1_new, mu2_new = new["mu1_Jmol"], new["mu2_Jmol"]

    mu1_rel = abs(mu1_new - mu1_ref) / abs(mu1_ref)
    mu2_rel = abs(mu2_new - mu2_ref) / abs(mu2_ref)
    status = "PASS" if (mu1_rel <= TOL_REL and mu2_rel <= TOL_REL) else "FAIL"
    print(f"[{label}] T={T:.3f}K rho_mol={rho_mol:.3f} mol/m3  status={status}")
    print(f"    mu1: ref={mu1_ref:.10g} J/mol   new={mu1_new:.10g} J/mol   rel_err={mu1_rel:.3g}")
    print(f"    mu2: ref={mu2_ref:.10g} J/mol   new={mu2_new:.10g} J/mol   rel_err={mu2_rel:.3g}")
    return status == "PASS", dict(label=label, mu1_rel=mu1_rel, mu2_rel=mu2_rel)


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1, mw2 = mw_from_json(d1), mw_from_json(d2)
    x1 = w1_to_x1(W1, mw1, mw2)
    print(f"x1={x1:.10f}\n")

    T_SAT = 300.0
    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0
    bubble_ref = solve_bubble_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    dew_ref = solve_dew_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    RHO_LIQ = 1.10 * bubble_ref["rho_l_molm3"]
    RHO_VAP = 0.50 * dew_ref["rho_v_molm3"]

    all_pass = True
    for label, T, rho in [
        ("liquid_direct_subcooled", T_SAT, RHO_LIQ),
        ("vapor_direct_superheated", T_SAT, RHO_VAP),
        ("supercritical_direct", 389.5, 3000.0),
    ]:
        ok, _ = compare_state(label, T, rho, x1, mw1, mw2)
        all_pass = all_pass and ok
        print()

    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e} relative)")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
