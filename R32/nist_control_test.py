"""
nist_control_test.py -- control experiment, 2026-09-25.

QUESTION: NIST fails to converge at 15 C (and oscillates at 10 C) in the
IDAES 2.13.0 / arm64 macOS environment, but converges 4/4 in both
IDAES 2.12.0 / arm64 and IDAES 2.13.0 / x86_64 Linux. NIST is the
reference method with published Shomate coefficients and critical
properties -- it should not fail. Is the failure caused by

    (a) the 09/24 chain (phase_1_..._0924 / vapor_compression_..._0924), or
    (b) the environment, independent of which chain is used?

METHOD: run NIST alone, at all four ambients, COLD START every time (no
warmstart seeding at all, so the monotonic-chase oscillation is removed
as a variable), through BOTH chains:

    _0910  -- the original, untouched, pre-09/19 chain
    _0924  -- the chain carrying the 09/19 delivery

NIST's parameter row is byte-identical between the two configs (verified
by diff: the files differ nowhere outside METHODS_RAW and the docstring),
so any difference in outcome here cannot come from the parameters.

READING THE RESULT:
  both chains FAIL at 15 C   -> environmental. Predates the 09/19 work
                                entirely; the _0924 chain is exonerated.
  only _0924 fails           -> something in the 09/19 chain after all;
                                report back, this would be a real bug.
  both chains PASS           -> the failure needs the warmstart chase to
                                appear, i.e. it is the chase that breaks
                                NIST, not the cold solve.

Run from the R32 directory with the environment you want to test active.

Author: Shilpa Narasimhan. Support: Claude AI.
"""

import sys
import traceback

AMBIENTS = [10, 15, 20, 25]
T_EVAP_C = -29.0
COND_APPROACH = 9.0


def run_one(CubicCycle, Mode, method, Tamb):
    """One cold-start solve. Mirrors phase6_final's run_cubic() exactly,
    minus the warmstart_points argument."""
    vc = CubicCycle("R32", compressor_efficiency=0.9999,
                    mode=Mode.IMPROVED_TPX, method=method)
    vc.specify_initial_conditions(low_side_temperature=T_EVAP_C,
                                  high_side_temperature=Tamb + COND_APPROACH)
    vc.initialize(verbose=False)
    vc.set_specifications(ambient_temperature=Tamb,
                          condenser_approach=COND_APPROACH,
                          evap_sat_temperature=T_EVAP_C,
                          superheating=0, subcooling=0,
                          max_pressure_ratio=10)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True,
                                     optimize=False)
    return cop, converged


def test_chain(label, module_name):
    print(f"\n{'='*70}")
    print(f"  {label}   ({module_name})")
    print(f"{'='*70}")
    try:
        mod = __import__(module_name)
        CubicCycle = mod.SimpleVaporCompressionCycle
        Mode = mod.Mode
    except Exception as e:
        print(f"  could not import: {type(e).__name__}: {e}")
        return {}

    results = {}
    for Tamb in AMBIENTS:
        try:
            cop, conv = run_one(CubicCycle, Mode, "NIST", Tamb)
            results[Tamb] = cop if conv else None
            status = f"COP = {cop:.4f}" if conv else "FAILED (not converged)"
        except Exception as e:
            results[Tamb] = None
            status = f"FAILED ({type(e).__name__}: {str(e)[:60]})"
        print(f"   T_amb = {Tamb:>2} C   {status}")
    n = sum(1 for v in results.values() if v is not None)
    print(f"   -> {n}/4 ambients converged")
    return results


if __name__ == "__main__":
    import idaes, pyomo, platform
    print(f"idaes {idaes.__version__} | pyomo {pyomo.version.version} | "
          f"{platform.machine()} {platform.system()}")
    print("NIST only. Cold start every solve -- no warmstart seeding.")

    old = test_chain("ORIGINAL 09/10 chain", "vapor_compression_cubic_refstate_0910")
    new = test_chain("09/24 chain (09/19 data)", "vapor_compression_cubic_refstate_0924")

    print(f"\n{'='*70}")
    print("  VERDICT")
    print(f"{'='*70}")
    print(f"   {'T_amb':>7}{'_0910':>14}{'_0924':>14}{'same?':>9}")
    for Tamb in AMBIENTS:
        a, b = old.get(Tamb), new.get(Tamb)
        fa = f"{a:.4f}" if a is not None else "FAILED"
        fb = f"{b:.4f}" if b is not None else "FAILED"
        if a is None and b is None:
            same = "yes"
        elif a is None or b is None:
            same = "NO"
        else:
            same = "yes" if abs(a - b) < 1e-6 else "NO"
        print(f"   {Tamb:>7}{fa:>14}{fb:>14}{same:>9}")

    na = sum(1 for v in old.values() if v is not None)
    nb = sum(1 for v in new.values() if v is not None)
    print()
    if na == nb == 4:
        print("   Both chains solve NIST 4/4 cold. The failure seen in")
        print("   phase6_final therefore requires the warmstart chase --")
        print("   the chase is what breaks NIST, not the cold solve.")
    elif na == nb:
        print("   Both chains fail identically. ENVIRONMENTAL -- independent")
        print("   of the 09/19 data, and predates that work entirely.")
    else:
        print("   The two chains DIFFER. This points at the 09/24 chain")
        print("   itself and needs investigating before any result is used.")
