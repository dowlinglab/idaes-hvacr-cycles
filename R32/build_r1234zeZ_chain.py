"""
build_r1234zeZ_chain.py -- 10/08/2026.

Derives the R-1234ze(Z) IDAES chain from the R-32 chain by COPYING the _1001
files and applying a small, explicit set of in-place edits. Nothing is retyped:
every line that is not listed below is byte-identical to its R-32 original, so
a diff shows exactly what changed and nothing else.

    phase_1_cubic_eos_validation_refstate_1001.py
        -> phase_1_cubic_eos_validation_refstate_R1234zeZ_1008.py
    vapor_compression_cubic_refstate_1001.py
        -> vapor_compression_cubic_refstate_R1234zeZ_1008.py
    phase6_final_1001.py
        -> phase6_final_R1234zeZ_1008.py
    cycle_algebraic_1001.py
        -> cycle_algebraic_R1234zeZ_1008.py

Run it from the R32 directory. It refuses to run if an edit does not match the
expected number of occurrences, so a silent partial patch cannot happen.
"""

import re
import shutil
import sys

SUF_OLD = "_1001"
SUF_NEW = "_R1234zeZ_1008"

# ---- the new fluid's numbers -----------------------------------------------
MW_KG_MOL = 114.0416e-03          # R1234ZEZ.FLD
TC = 423.27                       # K      R1234ZEZ.FLD
PC_PA = 3530.6e3                  # Pa     R1234ZEZ.FLD (3530.6 kPa)
OMEGA = 0.327                     # R1234ZEZ.FLD
SHOMATE = dict(A=-25.715129, B=424.511818, C=-116.375285,
               D=-256.761855, E=0.895699)
ANTOINE = dict(A=4.33220, B=1106.78020, C=-27.20229)
R_GAS = 8.3144598                 # R1234ZEZ.FLD !Gas constant

COMP_OLD, COMP_NEW = "R32", "R1234zeZ"

errors = []


def patch(text, old, new, count, tag):
    """Replace `old` -> `new`, asserting it occurs exactly `count` times."""
    n = text.count(old)
    if n != count:
        errors.append(f"{tag}: expected {count} occurrence(s) of {old!r}, found {n}")
        return text
    return text.replace(old, new)


def patch_re(text, pattern, repl, count, tag, flags=re.M):
    out, n = re.subn(pattern, repl, text, flags=flags)
    if n != count:
        errors.append(f"{tag}: expected {count} regex match(es) for {pattern!r}, got {n}")
        return text
    return out


def add_note(text, tag, title, body):
    """Append a dated note INSIDE the module docstring.

    Learned the hard way: inserting before the opening triple-quote leaves the
    prose as bare module-level text and the file stops parsing. The note goes
    immediately before the docstring's CLOSING triple-quote instead.
    """
    m = re.match(r'\s*(?:[rubRUB]{0,2})("""|\'\'\')', text)
    if not m:
        errors.append(f"{tag}: no module docstring found, note not added")
        return text
    q = m.group(1)
    close = text.find(q, m.end())
    if close == -1:
        errors.append(f"{tag}: unterminated module docstring, note not added")
        return text
    ins = f"\n--- 10/08/2026 -- {title} ---\n{body}\n"
    return text[:close] + ins + text[close:]


# =============================================================================
# 1. phase_1  (the config module: this is the ONLY place parameters are entered)
# =============================================================================
src = "phase_1_cubic_eos_validation_refstate_1001.py"
dst = f"phase_1_cubic_eos_validation_refstate{SUF_NEW}.py"
shutil.copy(src, dst)
t = open(dst, encoding="utf-8").read()

NEW_METHODS = '''METHODS_RAW = {{
    # ---------------------------------------------------------------------
    # PR -- the NIST/REFPROP row, and the only method defined so far for
    # this fluid. Tc, Pc and omega are taken verbatim from the REFPROP fluid
    # file R1234ZEZ.FLD (NIST SRD 23 v10.0, compiled by E.W. Lemmon):
    #       423.27    !Critical temperature [K]
    #       3530.6    !Critical pressure [kPa]
    #       0.327     !Acentric factor
    #
    # NOTE Pc: Higashi et al. (2013) MEASURED 3536 kPa. NIST refitted it to
    # 3530.6 kPa when building the equation of state. The NIST number is the
    # one used here, because that is the equation everything else is scored
    # against.
    #
    # The Shomate A-E are NOT published. The NIST WebBook has no record for
    # CAS 29118-25-0 at all (checked by ID, by name and by formula), so there
    # is no Shomate set to copy. A-E here are least-squares refitted to the
    # ideal-gas Cp carried by the Akasaka & Lemmon (2019) equation itself --
    # a Planck-Einstein form with Einstein temperatures 20 K and 1335 K --
    # over 238.13-372.61 K, the experimental range. The refit reproduces that
    # Cp to 0.0047 %. verify_akasaka_lemmon_SI.py checks this against the
    # paper's own Supporting Information and must be kept passing.
    #
    # Say this out loud in any report: A-E are this work, not a citation.
    # ---------------------------------------------------------------------
    "PR": {{"Pc": {PC:.6g}, "Tc": {TC:.6g}, "omega": {OM:.6g},
           "A": {A:.6f}, "B": {B:.6f}, "C": {C:.6f},
           "D": {D:.6f}, "E": {E:.6f},
           "F": 0.0, "G": 0.0}},
}}'''.format(PC=PC_PA, TC=TC, OM=OMEGA, **SHOMATE)

# METHODS_RAW is a single top-level dict: from the line starting with
# "METHODS_RAW = {" through the first line that is exactly "}".
t = patch_re(t, r"^METHODS_RAW = \{.*?^\}", NEW_METHODS.replace("\\", "\\\\"),
             1, "phase_1/METHODS_RAW", flags=re.M | re.S)

t = patch_re(t, r"^MW = 52\.024e-03.*$",
             f"MW = {MW_KG_MOL:.6e} # kg/mol, R-1234ze(Z) -- R1234ZEZ.FLD "
             f"!Molar mass 114.0416 g/mol", 1, "phase_1/MW")
t = patch_re(t, r"^#\(NIST WebBook: 52\.024 g/mol\)\.$",
             "# (NOT from the NIST WebBook -- it has no record for CAS "
             "29118-25-0; this is REFPROP.)", 1, "phase_1/MW-comment")

t = patch_re(t, r"^R = 8\.314 .*$",
             f"R = {R_GAS} # Real gas constant, J/mol/K -- R1234ZEZ.FLD "
             f"!Gas constant (R-32 chain used 8.314, 0.0055 % low)",
             1, "phase_1/R")

t = patch_re(t, r"^Antoine_A = [\d.]+$", f"Antoine_A = {ANTOINE['A']:.5f}",
             1, "phase_1/Antoine_A")
t = patch_re(t, r"^Antoine_B = [\d.]+\s+# K$",
             f"Antoine_B = {ANTOINE['B']:.5f}   # K", 1, "phase_1/Antoine_B")
t = patch_re(t, r"^Antoine_C = -?[\d.]+\s+# K$",
             f"Antoine_C = {ANTOINE['C']:.5f}   # K", 1, "phase_1/Antoine_C")
t = patch_re(t, r"^# Antoine: log10\(P_bar\).*$",
             "# Antoine: log10(P_bar) = A - B/(T_K + C) -- saturation-pressure "
             "init guess,\n# refitted 10/08/2026 to the Akasaka-Lemmon saturation "
             "curve over 238-420 K\n# (max 1.8 % error, <0.7 % over the cycle's "
             "range). Initialisation only.", 1, "phase_1/Antoine-comment")

# component key
t = patch(t, f'"{COMP_OLD}":{{', f'"{COMP_NEW}":{{', 1, "phase_1/component")
t = patch(t, f'("{COMP_OLD}", "{COMP_OLD}")', f'("{COMP_NEW}", "{COMP_NEW}")',
          1, "phase_1/PR_kappa")
t = patch(t, f'mole_frac_comp["{COMP_OLD}"]', f'mole_frac_comp["{COMP_NEW}"]',
          1, "phase_1/mole_frac")

# state bounds: the evaporator runs at 0.169 bar, below the old 1e4 Pa floor,
# and Tc is 423 K rather than 351 K.
t = patch(t, '"temperature": (200.0, 300.0, 450.0, pyunits.K)',
          '"temperature": (220.0, 320.0, 470.0, pyunits.K)',
          1, "phase_1/T-bounds")
t = patch(t, '"pressure": (1.0e4, 1.0e5, 1.0e7, pyunits.Pa)',
          '"pressure": (5.0e3, 1.0e5, 1.0e7, pyunits.Pa)',
          1, "phase_1/P-bounds")

t = add_note(t, "phase_1", "retargeted to R-1234ze(Z)",
             "Derived from phase_1_cubic_eos_validation_refstate_1001.py by\n"
             "build_r1234zeZ_chain.py. Changed: METHODS_RAW (one PR row from\n"
             "REFPROP), MW, R, the Antoine init guess, the component key\n"
             "R32 -> R1234zeZ, and the temperature/pressure state bounds.\n"
             "Everything else is byte-identical to the R-32 original -- diff\n"
             "the two files to confirm.")

open(dst, "w", encoding="utf-8").write(t)
print(f"wrote {dst}")


# =============================================================================
# 2. vapor_compression  (cycle model: only the import and the component key)
# =============================================================================
src = "vapor_compression_cubic_refstate_1001.py"
dst = f"vapor_compression_cubic_refstate{SUF_NEW}.py"
shutil.copy(src, dst)
t = open(dst, encoding="utf-8").read()

t = patch(t, f"from phase_1_cubic_eos_validation_refstate{SUF_OLD} import",
          f"from phase_1_cubic_eos_validation_refstate{SUF_NEW} import",
          1, "vapor_compression/import")

# Every CoolProp call in this file uses self.fluid_name (set from the
# constructor), so "R32" appears here ONLY as the Pyomo component key.
n_comp = t.count(f'mole_frac_comp["{COMP_OLD}"]')
t = t.replace(f'mole_frac_comp["{COMP_OLD}"]', f'mole_frac_comp["{COMP_NEW}"]')
n_dict = t.count(f'"mole_frac_comp": {{"{COMP_OLD}": 1.0}}')
t = t.replace(f'"mole_frac_comp": {{"{COMP_OLD}": 1.0}}',
              f'"mole_frac_comp": {{"{COMP_NEW}": 1.0}}')
print(f"   component key rewritten: {n_comp} mole_frac_comp[...] + {n_dict} dict literal")
if n_comp == 0:
    errors.append("vapor_compression: no mole_frac_comp['R32'] found")

open(dst, "w", encoding="utf-8").write(t)
print(f"wrote {dst}")


# =============================================================================
# 3. phase6  (the sweep)
# =============================================================================
NEW_METHODS_LIST = '''METHODS = ["PR"]
# Only one method exists for R-1234ze(Z) so far: the NIST/REFPROP PR row.
# "Helmholtz" is NOT in this list, and must not be added -- see HELMHOLTZ
# note at the import above. When the Colon group delivers R-1234ze(Z)
# parameter sets, add their names here AND their rows to METHODS_RAW in
# phase_1_cubic_eos_validation_refstate_R1234zeZ_1008.py.'''

NEW_CONDITIONS = '''# =========================================================================
# OPERATING CONDITIONS -- the only place the duty is defined.
# =========================================================================
# R-32 ran T_evap = -29 C with T_cond = T_amb + 9 and ambients 10/15/20/25.
# That point is UNUSABLE for R-1234ze(Z). Its saturation pressure at -29 C
# is 0.169 bar -- about 0.17 atm -- and the chain failed to converge there
# at every ambient and every solver setting tried (16 attempts, 10/08/2026:
# default, tol=1e-6, tol=1e-5 + constr_viol_tol, gradient-based scaling).
#
# R-1234ze(Z) is a HIGH-TEMPERATURE fluid: NBP 282.878 K (9.7 C), Tc 423.27 K
# (150.1 C). Petr & Raabe, Energy 93 (2015) 266-274, run it between a
# near-ambient condenser and a 100-250 C heat source with a 10 K pinch; the
# window below sits inside that. The chain converges 4/4 here.
#
# AMBIENTS keeps its name for a minimal diff, but for a high-temperature heat
# pump it is the heat-SINK temperature, not the weather.
AMBIENTS = [60, 70, 80, 90]   # heat-sink temperature [C]   (was 10/15/20/25)
T_EVAP_C = 30.0               # evaporator sat. temperature [C]   (was -29)
COND_APPROACH = 10.0          # condenser approach [K]  (was 9; 10 K matches
                              #   the Petr & Raabe pinch)
SUPERHEAT = 0                 # K   (unchanged)
SUBCOOL = 0                   # K   (unchanged)
# => T_cond = 70 / 80 / 90 / 100 C'''

HELM_GUARD = '''# HELMHOLTZ IS NOT AVAILABLE FOR THIS FLUID.
# IDAES's general_helmholtz package registers 'r1234ze' -- but that is
# R-1234ze(E), the TRANS isomer, Tc = 382.513 K. Our fluid is the CIS
# isomer, R-1234ze(Z), Tc = 423.27 K. Both have molar mass 114.0416 g/mol,
# so a wrong-isomer run would look entirely plausible and be entirely wrong.
# Checked 10/08/2026: registered_components() returns
#   butane, co2, h2o, isobutane, nh3, propane, r1234ze, r125, r134a,
#   r227ea, r32
# -- R-1234ze(Z) is not among them.
#
# The Helmholtz reference for this fluid therefore comes from CoolProp
# (Akasaka & Lemmon), not from IDAES. See compare_cp_methods_R1234zeZ_1008.py
# and verify_akasaka_lemmon_SI.py. run_helm() below is left in place but is
# never reached, because "Helmholtz" is not in METHODS.
try:
    from vapor_compression import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
except Exception as _e:          # pragma: no cover
    HelmCycle = HelmMode = None
    print(f"note: IDAES Helmholtz cycle unavailable ({_e}); not needed for "
          f"R-1234ze(Z).")'''

src = "phase6_final_1001.py"
dst = f"phase6_final{SUF_NEW}.py"
shutil.copy(src, dst)
t = open(dst, encoding="utf-8").read()

t = patch(t, f"from vapor_compression_cubic_refstate{SUF_OLD} import",
          f"from vapor_compression_cubic_refstate{SUF_NEW} import",
          1, "phase6/import")

t = patch_re(t, r'^FLUID = "R32"$', 'FLUID = "R1234ze(Z)"', 1, "phase6/FLUID")
t = patch_re(t, r"^M_R32 = 0\.052024.*$",
             f"M_R32 = {MW_KG_MOL:.7f}  # kg/mol, R-1234ze(Z) (name kept for "
             f"minimal diff)", 1, "phase6/MW")

t = patch_re(t, r"^METHODS = \[.*?\]$", NEW_METHODS_LIST, 1,
             "phase6/METHODS", flags=re.M | re.S)

# ---- operating conditions: hoisted to named constants ----------------------
# In the R-32 original the duty is written as bare literals in four places
# inside run_cubic() and run_helm(). Here they become five module-level
# constants so there is exactly ONE place to read or change the duty.
t = patch_re(t, r"^AMBIENTS = \[10, 15, 20, 25\]$", NEW_CONDITIONS, 1,
             "phase6/AMBIENTS")

n = t.count("Tcond_sat = Tamb + 9")
t = t.replace("Tcond_sat = Tamb + 9", "Tcond_sat = Tamb + COND_APPROACH")
if n != 2:
    errors.append(f"phase6: expected 2 'Tcond_sat = Tamb + 9', found {n}")

n = t.count("low_side_temperature=-29,")
t = t.replace("low_side_temperature=-29,", "low_side_temperature=T_EVAP_C,")
if n != 2:
    errors.append(f"phase6: expected 2 'low_side_temperature=-29', found {n}")

n = t.count("condenser_approach=9,")
t = t.replace("condenser_approach=9,", "condenser_approach=COND_APPROACH,")
if n != 2:
    errors.append(f"phase6: expected 2 'condenser_approach=9', found {n}")

n = t.count("evap_sat_temperature=-29, superheating=0, subcooling=0,")
t = t.replace("evap_sat_temperature=-29, superheating=0, subcooling=0,",
              "evap_sat_temperature=T_EVAP_C, superheating=SUPERHEAT, "
              "subcooling=SUBCOOL,")
if n != 2:
    errors.append(f"phase6: expected 2 'evap_sat_temperature=-29 ...', found {n}")

# ---- the Helmholtz branch cannot be used for this fluid --------------------
t = patch(t,
          "from vapor_compression import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode",
          HELM_GUARD, 1, "phase6/HelmCycle-import")

# Pressure-ratio cap: at T_evap = -29 C this fluid needs up to 14.2, and the
# study is theoretical, so the cap is lifted rather than binding.
n = t.count("max_pressure_ratio=10")
t = t.replace("max_pressure_ratio=10", "max_pressure_ratio=1e6")
print(f"   max_pressure_ratio lifted in {n} place(s)")
if n != 2:
    errors.append(f"phase6: expected 2 max_pressure_ratio, found {n}")

open(dst, "w", encoding="utf-8").write(t)
print(f"wrote {dst}")


# =============================================================================
# 4. cycle_algebraic
# =============================================================================
src = "cycle_algebraic_1001.py"
dst = f"cycle_algebraic{SUF_NEW}.py"
shutil.copy(src, dst)
t = open(dst, encoding="utf-8").read()
t = patch_re(t, r"^MW = 52\.024e-3.*$",
             f"MW = {MW_KG_MOL:.6e}     # kg/mol, R-1234ze(Z)",
             1, "cycle_algebraic/MW")
t = add_note(t, "cycle_algebraic", "INCOMPLETE -- READ THIS BEFORE USING",
             "This file is NOT retargeted to R-1234ze(Z). Only the molar mass\n"
             "was changed. Every Fluid row below still carries R-32 critical\n"
             "properties and R-32 Shomate coefficients, and the __main__ block\n"
             "still sweeps the R-32 duty (T_evap = -29 C, ambients 10/15/20/25).\n"
             "Running it as-is gives R-32 numbers with an R-1234ze(Z) molar mass,\n"
             "which is meaningless.\n"
             "\n"
             "It is shipped only so the copy set is complete. The algebraic\n"
             "PR-vs-Helmholtz comparison for this fluid already exists and is\n"
             "verified -- use compare_cp_methods_R1234zeZ_1008.py instead.\n"
             "\n"
             "To finish this file: replace the Fluid rows with the REFPROP PR row\n"
             "from phase_1_cubic_eos_validation_refstate_R1234zeZ_1008.py, and set\n"
             "the __main__ sweep to T_evap = 30 C, T_cond = 70/80/90/100 C.")
open(dst, "w", encoding="utf-8").write(t)
print(f"wrote {dst}")


# =============================================================================
if errors:
    print(f"\n{len(errors)} PATCH ERROR(S):")
    for e in errors:
        print("   -", e)
    sys.exit(1)
print("\nall patches applied cleanly.")
