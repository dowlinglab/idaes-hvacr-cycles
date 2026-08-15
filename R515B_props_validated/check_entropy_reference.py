"""
Check whether the model's own entropy at Honeywell's stated reference state
matches Honeywell's stated value.

Purpose
-------
Honeywell's Solstice N15 (R-515B) TDS p-H chart states explicitly (small print,
bottom of chart): "Reference State: h = 200 kJ/kg, s = 1.00 kJ/kg-K; sat. liq.
at 0 C". pressure_validated_model.py already uses this exact reference state
for its enthalpy chart_offset() fix: T_ref=273.15 K, rho_ref=1258.4 kg/m^3
(saturated liquid at 0 C), h_ref=200 kJ/kg.

_mix_entropy_direct() in mixture_isentrope_validation.py computes RAW entropy
straight from the Helmholtz EOS identity (s/R = tau*alpha_tau - alpha) with NO
reference-state correction applied -- unlike enthalpy, which already has a
documented offset fix. This script evaluates the model's raw entropy at that
exact same reference state and reports how far off it is from Honeywell's
stated 1.00 kJ/kg-K -- if there's a real, confirmed offset, that value is
exactly the additive correction needed to align every isentrope label with
Honeywell's convention.

Usage
-----
Run from R515B_props_validated/ (same directory as mixture_isentrope_validation.py
and its own copy of linear_model_codex.py):
    python3 check_entropy_reference.py

Does not modify anything -- read-only diagnostic.
"""

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    mix_state,
    _mix_entropy_direct,
    w1_to_x1,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911

# Honeywell's own stated reference state (matches pressure_validated_model.py's
# ChartReference: T_ref_K=273.15, h_ref_kJkg=200.0, rho_ref_kgm3=1258.4).
T_REF_K = 273.15
RHO_REF_KGM3 = 1258.4
H_REF_KJKG = 200.0
S_REF_KJKGK = 1.00  # Honeywell's stated reference entropy

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
z1 = w1_to_x1(W1, mw1, mw2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2  # kg/mol

rho_ref_molm3 = RHO_REF_KGM3 / mw_mix

print(f"Reference state: T={T_REF_K} K, rho={RHO_REF_KGM3} kg/m3 ({rho_ref_molm3:.4f} mol/m3), z1={z1:.6f}")
print(f"mw_mix = {mw_mix*1000:.4f} g/mol\n")

st = mix_state(d1, d2, T_REF_K, rho_ref_molm3, z1)
s_jmolK = _mix_entropy_direct(d1, d2, T_REF_K, rho_ref_molm3, z1)

h_kjkg_raw = st.h_jmol / mw_mix / 1000.0
s_kjkgK_raw = s_jmolK / mw_mix / 1000.0

print(f"Raw model enthalpy at reference state: {h_kjkg_raw:.4f} kJ/kg")
print(f"  Honeywell states h_ref = {H_REF_KJKG} kJ/kg")
print(f"  Implied h offset (h_ref - raw) = {H_REF_KJKG - h_kjkg_raw:+.4f} kJ/kg")
print(f"  (sanity check: pressure_validated_model.py's chart_offset() should")
print(f"   already report a consistent value close to this)\n")

print(f"Raw model entropy at reference state: {s_kjkgK_raw:.6f} kJ/(kg*K)")
print(f"  Honeywell states s_ref = {S_REF_KJKGK} kJ/(kg*K)")
print(f"  Implied s offset (s_ref - raw) = {S_REF_KJKGK - s_kjkgK_raw:+.6f} kJ/(kg*K)")
print(f"  Implied s offset in J/(mol*K) = {(S_REF_KJKGK - s_kjkgK_raw) * mw_mix * 1000.0:+.4f} J/(mol*K)")
