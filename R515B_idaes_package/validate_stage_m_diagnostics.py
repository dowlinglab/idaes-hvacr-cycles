"""
Stage M: IDAES structural/numerical diagnostics on the R-515B property
package, using `idaes.core.util.diagnostics_tools.diagnostics_toolbox.
DiagnosticsToolbox` (the modern 2.12.0 location -- the older
`idaes.core.util.model_diagnostics.DiagnosticsToolbox` import path is
deprecated and warns on import).

Per spec: Stage M is a structural-tests gate that must pass (or have any
findings explicitly triaged/documented) before Stage N (vapor_compression.py
integration). This script runs the toolbox against:
  1. A single R515BStateBlock, unfixed (structural-only checks -- the
     toolbox's structural methods are designed to run before a model is
     square/solved).
  2. The same StateBlock after fixing the 3 public state vars and solving
     (structural + numerical checks on a converged square problem), at one
     representative two-phase point.

`assert_no_structural_warnings()`/`assert_no_numerical_warnings()` raise on
any WARNING-level finding; this script catches and reports rather than
letting a bare exception propagate, per the "never hide a failure, but
also never silently swallow one" spirt of rule 74/89 -- any finding is
printed in full and recorded as a FAIL for this stage, not hidden.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import ConcreteModel  # noqa: E402
from idaes.core import FlowsheetBlock  # noqa: E402
from idaes.core.util.diagnostics_tools.diagnostics_toolbox import DiagnosticsToolbox  # noqa: E402

import r515b_helmholtz_core as core  # noqa: E402
from r515b_property_package import R515BParameterBlock  # noqa: E402


def build_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = R515BParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
    return m


def main():
    all_pass = True

    # ---- Part 1: structural diagnostics on the freshly-built, SQUARE block ----
    # NOTE (self-caught test-design bug, fixed here, not the package): a
    # first draft ran `assert_no_structural_warnings()` on the block with
    # its 3 public state vars still UNFIXED (DOF=2, matching Section 24c's
    # already-documented convention). That correctly reported "2 Degrees
    # of Freedom" / "Structural singularity found" WARNINGS -- but those
    # warnings are the EXPECTED, correct description of a deliberately
    # under-determined block, not a defect. `DiagnosticsToolbox`'s
    # structural-warning assertion is meant to run on a SQUARE (DOF=0)
    # problem, matching how it's actually used downstream (a StateBlock
    # inside a solved flowsheet always has its state vars fixed). Fixed by
    # fixing the 3 public state vars first, exactly as Stage L part 3's own
    # validation already established is the correct precondition.
    print("=" * 70)
    print("Part 1: structural diagnostics, state vars FIXED (square, DOF=0), not yet solved")
    print("=" * 70)
    m1 = build_flowsheet()
    blk1 = m1.fs.state[0]
    blk1.flow_mass.fix(1.0)
    blk1.pressure.fix(5.0e5)
    blk1.enth_mass.fix(3.0e5)
    dt1 = DiagnosticsToolbox(model=blk1)
    dt1.report_structural_issues()
    try:
        dt1.assert_no_structural_warnings()
        print("[structural, square/unsolved] PASS -- no structural warnings")
    except AssertionError as e:
        print(f"[structural, square/unsolved] FAIL -- structural warnings present:\n{e}")
        all_pass = False

    # ---- Part 2: structural + numerical diagnostics on a solved, square block ----
    print("\n" + "=" * 70)
    print("Part 2: structural + numerical diagnostics, solved (two-phase point)")
    print("=" * 70)

    d1 = core.load_idaes_helmholtz_json(core.FLUID1)
    d2 = core.load_idaes_helmholtz_json(core.FLUID2)
    x1_val = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1_val * rhoc1 + (1.0 - x1_val) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    sat300 = core.solve_pseudopure_saturation_at_t(d1, d2, x1_val, 300.0, rho_l_seed0, rho_v_seed0)
    P_test = sat300["P_Pa"]
    h_l_sat, h_v_sat = sat300["h_l_Jmol"], sat300["h_v_Jmol"]
    H_molar_test = h_l_sat + 0.5 * (h_v_sat - h_l_sat)

    m2 = build_flowsheet()
    blk2 = m2.fs.state[0]
    mw_mix = float(mw1 * x1_val + mw2 * (1.0 - x1_val))
    blk2.flow_mass.fix(1.0)
    blk2.pressure.fix(P_test)
    blk2.enth_mass.fix(H_molar_test / mw_mix)

    m2.fs.state.initialize(outlvl=0)

    dt2 = DiagnosticsToolbox(model=blk2)
    dt2.report_structural_issues()
    try:
        dt2.assert_no_structural_warnings()
        print("[structural, solved] PASS -- no structural warnings")
    except AssertionError as e:
        print(f"[structural, solved] FAIL -- structural warnings present:\n{e}")
        all_pass = False

    dt2.report_numerical_issues()
    try:
        dt2.assert_no_numerical_warnings()
        print("[numerical, solved] PASS -- no numerical warnings")
    except AssertionError as e:
        print(f"[numerical, solved] FAIL -- numerical warnings present:\n{e}")
        all_pass = False

    print(f"\nOVERALL Stage M (structural/numerical diagnostics): {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
