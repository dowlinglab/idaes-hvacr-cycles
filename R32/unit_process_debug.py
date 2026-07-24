"""
unit_process_debug.py -- reusable diagnostic tool for Phase 3b initialize()
failures.

Generalizes the ad-hoc inspection scripts written during Phase 3b debugging
(evaporator, compressor, condenser, expansion valve) into one reusable tool:
build the cycle, run vc.initialize() (catching any failure instead of
crashing), then print the state (T, P, flow, composition, model-internal
saturation temperature, phase fractions, compressibility factors) for every
unit's inlet and outlet -- all in one pass, regardless of which unit is
currently the first to fail.

Usage:
    python3 unit_process_debug.py
    (edit low_side_temperature / high_side_temperature / compressor_efficiency
    below, or import inspect_state_block() / run_and_diagnose() elsewhere)

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode


def inspect_state_block(name, blk):
    """Print T, P, flow, composition, saturation temperature, phase
    fractions, and compressibility factors for one state block.

    Each property is wrapped in its own try/except -- some of these
    (e.g. temperature_bubble) are build-on-demand and may not have a value
    yet if nothing has triggered their solve, and we want one missing
    property to print an error line, not crash the whole inspection.
    """
    print(f"\n--- {name} ---")
    try:
        print(f"  T = {value(blk.temperature):.3f} K, P = {value(blk.pressure):.1f} Pa")
    except Exception as e:
        print(f"  T/P: {e}")

    try:
        print(f"  flow_mol = {value(blk.flow_mol):.4f}, "
              f"mole_frac_comp[R32] = {value(blk.mole_frac_comp['R32']):.6f}")
    except Exception as e:
        print(f"  flow/composition: {e}")

    try:
        print(f"  entr_mol = {value(blk.entr_mol):.4f} J/mol/K")
    except Exception as e:
        print(f"  entr_mol: {e}")

    try:
        tbub = value(blk.temperature_bubble["Vap", "Liq"])
        tdew = value(blk.temperature_dew["Vap", "Liq"])
        print(f"  tbub = {tbub:.3f} K, tdew = {tdew:.3f} K "
              f"(model-internal Tsat at this P; equal for a pure fluid)")
    except Exception as e:
        print(f"  tbub/tdew: {e}")

    try:
        vap_frac = value(blk.phase_frac["Vap"])
        liq_frac = value(blk.phase_frac["Liq"])
        print(f"  phase_frac: Vap = {vap_frac:.6f}, Liq = {liq_frac:.6f}, "
              f"sum = {vap_frac + liq_frac:.6f} "
              f"({'OK' if abs(vap_frac + liq_frac - 1.0) < 1e-4 else 'NOT 1.0 -- non-physical leftover'})")
    except Exception as e:
        print(f"  phase_frac: {e}")

    for ph in ["Vap", "Liq"]:
        try:
            print(f"  Z[{ph}] = {value(blk.compress_fact_phase[ph]):.4f}")
        except Exception as e:
            print(f"  Z[{ph}]: {e}")


def inspect_unit(unit_name, unit):
    """Inspect a unit's inlet and outlet control-volume state blocks."""
    inspect_state_block(f"{unit_name} INLET", unit.control_volume.properties_in[0])
    inspect_state_block(f"{unit_name} OUTLET", unit.control_volume.properties_out[0])


def run_and_diagnose(low_side_temperature=-29, high_side_temperature=29,
                      compressor_efficiency=0.9999):
    """Build the cycle, attempt vc.initialize(), then inspect every unit's
    inlet/outlet regardless of where (or whether) initialization failed.

    Returns the SimpleVaporCompressionCycle instance so the caller can keep
    poking at it afterward (e.g. vc.model.fs.<unit>...) in an interactive
    session.
    """
    vc = SimpleVaporCompressionCycle(
        "R32", compressor_efficiency=compressor_efficiency, mode=Mode.IMPROVED_TPX
    )
    vc.specify_initial_conditions(
        low_side_temperature=low_side_temperature,
        high_side_temperature=high_side_temperature,
    )

    print(f"T_init = {vc.T_init}")
    print(f"p_init = {vc.p_init}")

    try:
        vc.initialize(verbose=False)
        print("\nvc.initialize() SUCCEEDED\n")
    except Exception as e:
        print(f"\nvc.initialize() failed: {type(e).__name__}: {e}\n")

    for unit_name in ["evaporator", "compressor", "condenser", "expansion_valve"]:
        unit = getattr(vc.model.fs, unit_name)
        inspect_unit(unit_name, unit)

    return vc


if __name__ == "__main__":
    run_and_diagnose()
