"""
Stage D/F/G (partial) validation: confirm that `linear_model_codex.py`'s
`compute_table1_properties` -- an EXISTING, unmodified, legitimate
production dependency (explicitly permitted by MASTER TASK spec rule 6's
"legitimate production dependency" carve-out, NOT the oracle) -- reproduces
the oracle's own `mix_state` (P, h) and entropy (`_mix_entropy_direct`,
before the Honeywell rebase offset) at representative direct
("no-solver-either-side") states.

Why this matters
-----------------
Both `mixture_fully_validated.py` (the oracle) AND `linear_model_codex.py`
call the SAME shared, unmodified functions for the reducing functions
(`bell2023_Tred_vred`) and the combined ideal+residual+departure Helmholtz
kernel (`mixture_alpha0_alphar_derivs`) -- confirmed by reading the oracle's
own `from linear_model_codex import (...)` statement and its
`_mix_alpha_and_derivs`/`mix_state` bodies. This means Stage D (data),
Stage E (reducing functions), and the core Helmholtz-kernel part of Stage F
are numerically identical BY CONSTRUCTION, not merely "should match" --
there is a single shared code path, not two independent ports that could
drift apart. What this script actually checks is the DOWNSTREAM formulas:
does `compute_table1_properties`'s own P/h/s identities (written
independently of `mix_state`/`_mix_entropy_direct`, even though both pull
from the same alpha/derivative values) produce the same numbers.

Entropy convention note: `compute_table1_properties`'s `s_over_r` uses the
RAW (un-rebased) entropy convention (s/R = tau*alpha_tau - alpha, same
identity as the oracle's `_mix_entropy_direct` before its Honeywell
rebase). The oracle then adds `ENTROPY_REFERENCE_OFFSET_JMOLK` to match
Honeywell's stated reference state. This script compares the RAW (pre-
offset) entropy from both sources, then separately confirms that adding
the oracle's own offset constant reproduces the oracle's final rebased
value -- so the new production package knows exactly where that offset
needs to be applied.

This is READ-ONLY reference use of the oracle for comparison purposes only
(spec rule 9 carve-out). `linear_model_codex.py` (not the oracle) is what
the eventual production package will actually import at runtime.
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
    mix_state,
    _mix_entropy_direct,
    ENTROPY_REFERENCE_OFFSET_JMOLK,
    solve_bubble_at_t,
    solve_dew_at_t,
)
from linear_model_codex import compute_table1_properties  # noqa: E402

FLUID1, FLUID2 = "r1234ze", "r227ea"
W1 = 0.911
TOL_REL = 1e-8  # tier-1 "direct-equation" frozen tolerance (Section 1)


def compare_state(label, T, rho_mol, x1, mw1, mw2):
    ref = mix_state(load_idaes_helmholtz_json(FLUID1), load_idaes_helmholtz_json(FLUID2), T, rho_mol, x1)
    # IMPORTANT: `_mix_entropy_direct` already returns the Honeywell-REBASED
    # value internally (it adds ENTROPY_REFERENCE_OFFSET_JMOLK before
    # returning -- confirmed by reading its body directly, line 2569-2570:
    # `s_raw = R_u*(tau*alpha_tau-alpha); return s_raw + ENTROPY_REFERENCE_
    # OFFSET_JMOLK`). An earlier draft of this script incorrectly treated
    # its return value as pre-offset "raw" entropy and then added the offset
    # a SECOND time, producing a false ~1.4595 J/(mol*K) "mismatch" against
    # `compute_table1_properties` that was actually just double-counting the
    # rebase -- not a real physics discrepancy. Fixed here: back out the
    # offset explicitly to get the model's true raw (pre-Honeywell-rebase)
    # entropy for the apples-to-apples comparison against
    # `compute_table1_properties`, which has never been rebased.
    s_ref_rebased = _mix_entropy_direct(load_idaes_helmholtz_json(FLUID1), load_idaes_helmholtz_json(FLUID2), T, rho_mol, x1)
    s_ref_raw = s_ref_rebased - ENTROPY_REFERENCE_OFFSET_JMOLK

    rho_mass = rho_mol * (x1 * mw1 + (1.0 - x1) * mw2)
    new = compute_table1_properties(FLUID1, FLUID2, W1, T, rho_mass, fd_rel=1e-6)

    p_rel = abs(new["p_Pa"] - ref.p_pa) / abs(ref.p_pa)
    h_rel = abs(new["h_molar_Jmol"] - ref.h_jmol) / abs(ref.h_jmol)
    s_rel_raw = abs(new["s_molar_JmolK"] - s_ref_raw) / abs(s_ref_raw)

    status = "PASS" if (p_rel <= TOL_REL and h_rel <= TOL_REL and s_rel_raw <= TOL_REL) else "FAIL"
    print(f"[{label}] T={T:.3f}K rho_mol={rho_mol:.3f} mol/m3  status={status}")
    print(f"    P:  ref={ref.p_pa:.10g} Pa   new={new['p_Pa']:.10g} Pa   rel_err={p_rel:.3g}")
    print(f"    h:  ref={ref.h_jmol:.10g} J/mol   new={new['h_molar_Jmol']:.10g} J/mol   rel_err={h_rel:.3g}")
    print(f"    s(raw, pre-Honeywell-rebase):  ref={s_ref_raw:.10g} J/molK   new={new['s_molar_JmolK']:.10g} J/molK   rel_err={s_rel_raw:.3g}")
    print(f"    s(Honeywell-rebased, oracle's actual _mix_entropy_direct output): {s_ref_rebased:.10g} J/molK "
          f"(= raw + {ENTROPY_REFERENCE_OFFSET_JMOLK})")
    return status == "PASS", dict(label=label, T=T, rho_mol=rho_mol, p_rel=p_rel, h_rel=h_rel, s_rel_raw=s_rel_raw)


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
    results = []
    for label, T, rho in [
        ("liquid_direct_subcooled", T_SAT, RHO_LIQ),
        ("vapor_direct_superheated", T_SAT, RHO_VAP),
        ("supercritical_direct", 389.5, 3000.0),
    ]:
        ok, rec = compare_state(label, T, rho, x1, mw1, mw2)
        all_pass = all_pass and ok
        results.append(rec)
        print()

    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e} relative)")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
