"""

r515b_property_package.py -- Stage L part 3 (final piece): the actual IDAES
`PhysicalParameterBlock`/`StateBlockData` classes for R-515B, wiring the
already-validated native-Pyomo EOS kernel (`r515b_pyomo_eos.py`, Stage L
parts 1-2) and the smooth PH-flash system (`r515b_pyomo_eos.
mixture_ph_flash_residuals_expr`, this same Stage) into a real,
drop-in-compatible IDAES property package.

Architecture confirmed by direct inspection of the installed IDAES 2.12.0
source before writing this file (per spec rule 2's "consult the installed
version's real API, not memory/stale docs" requirement):
- `idaes.core.base.property_base.PhysicalParameterBlock`/`StateBlock`/
  `StateBlockData` are the true generic custom-property base classes
  (distinct from the C++-external-function `general_helmholtz` framework
  `vapor_compression.py` currently uses for pure fluids).
- `idaes.models.properties.examples.saponification_thermo.py` supplied the
  modern (2.12.0) API idioms used below: `declare_process_block_class`,
  `Phase()`/`Component()` objects, `define_metadata()` classmethod,
  `fix_state_vars`/`revert_state_vars` for `initialize()`/`release_state()`.
- `idaes.models.properties.general_helmholtz.helmholtz_state.
  HelmholtzStateBlockData` (the package `vapor_compression.py` uses today)
  supplied the EXACT drop-in contract this package must match for
  `Heater`/`Compressor`/`PressureChanger` (property-package-agnostic
  generic IDAES unit models) to work unchanged:
  - `StateVars.PH` + `AmountBasis.MASS` -> exactly 3 state Vars:
    `flow_mass`, `pressure`, `enth_mass`. `temperature`/`vapor_frac` are
    DERIVED (not independently fixable in PH mode) -- confirmed both from
    that source and from `vapor_compression.py`'s own Mode.PH usage
    (temperature bounds go through separate flowsheet-level Tmin/Tmax
    Params+Constraints, never `.temperature[0].setlb()` directly).
  - `PhaseType.MIX` -> a single `Phase()` named "Mix" exposed to the
    framework (`self.Mix = Phase()`), matching how `vapor_compression.py`'s
    existing pure-fluid package presents itself (confirmed: it does not
    pass `phase_presentation`, so it uses the MIX default).
  - `get_material_flow_terms(p,j)` = `flow_mass` (Mix phase, mass basis);
    `get_enthalpy_flow_terms(p)` = `flow_mass * enth_mass` (Mix phase).
  - `define_port_members()` overridden to expose `temperature`/
    `vapor_frac` on the Port in addition to the 3 state vars (both are
    proper Pyomo components -- an Expression and a Var respectively --
    which `StateBlock.build_port`'s `Reference(...)` mechanism accepts
    for either type).

R-515B specifics: single pseudo-pure `Component()` ("r515b") at the fixed
composition x1=x1_from_w1(w1=0.911) (per the 2026-08-17 confirmed
architectural decision), two-phase logic via the x1=y1=z1 near-azeotropic
simplification (per the 2026-08-18 rule-22 resolution) using the smooth
8-equation PH-flash system built and validated in `r515b_pyomo_eos.
mixture_ph_flash_residuals_expr` (Stage L part 3, this same 2026-08-18
session).

Traceability / what NOT to expect here: this file wires together already-
validated pieces -- it does not introduce new thermodynamic math. Every
equation inside `mixture_ph_flash_residuals_expr` (and everything it calls)
was independently validated against `r515b_helmholtz_core.py`'s SciPy
reference (itself validated exact vs. the oracle, Stages D-K) before this
file was written. Validated here (`validate_state_block_construction.py`):
that the StateBlockData actually CONSTRUCTS, has exactly 3 unfixed DOF
before state vars are fixed and 0 after, and reproduces
`flash_ph_pseudopure`'s answers once solved -- i.e. this file adds
plumbing/structure, and that plumbing is what gets tested.
"""

import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from pyomo.environ import Var, Param, Constraint, Expression, units as pyunits, value  # noqa: E402
from pyomo.common.config import ConfigValue  # noqa: E402

from idaes.core import (  # noqa: E402
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MaterialFlowBasis,
    Component,
)
from idaes.core.base.phases import Phase  # noqa: E402
from idaes.core.util.model_statistics import degrees_of_freedom  # noqa: E402
from idaes.core.util.initialization import fix_state_vars, revert_state_vars  # noqa: E402
import idaes.logger as idaeslog  # noqa: E402

import r515b_helmholtz_core as core  # noqa: E402
from r515b_pyomo_eos import mixture_ph_flash_residuals_expr  # noqa: E402

_log = idaeslog.getLogger(__name__)

# Molar-enthalpy scale sanity bounds (J/mol), used only to set generous Var
# bounds for numerical conditioning -- NOT a physical restriction. Derived
# from the observed 25-50 kJ/mol range around 250-380K seen throughout this
# project's validation work (Sections 5-13); widened generously on both
# sides since flowsheet conditions may run outside that narrow window.
H_MOLAR_MIN_JMOL = -50_000.0
H_MOLAR_MAX_JMOL = 150_000.0
T_MIN_K = 180.0
T_MAX_K = 420.0
RHO_MIN_MOLM3 = 1.0e-6
RHO_MAX_MOLM3 = 2.0e4
P_MIN_PA = 1.0e3
P_MAX_PA = 1.0e7


@declare_process_block_class("R515BParameterBlock")
class R515BParameterBlockData(PhysicalParameterBlock):
    """
    Parameter block for the R-515B pseudo-pure (x1=y1=z1 fixed) property
    package. Single "Mix" phase (PhaseType.MIX convention, matching the
    existing pure-fluid `HelmholtzParameterBlock` usage in
    `vapor_compression.py`), single pseudo-pure Component ("r515b").
    """

    def build(self):
        super().build()
        self._state_block_class = R515BStateBlock

        self.Mix = Phase()
        self.r515b = Component()

        d1 = core.load_idaes_helmholtz_json(core.FLUID1)
        d2 = core.load_idaes_helmholtz_json(core.FLUID2)
        x1_val = core.r515b_x1()
        mw1 = core.mw_from_json(d1)
        mw2 = core.mw_from_json(d2)
        mw_mix_val = x1_val * mw1 + (1.0 - x1_val) * mw2
        tc1_val = float(d1["basic"]["Tc"])
        tc2_val = float(d2["basic"]["Tc"])
        rhoc1_val = float(d1["basic"]["rhoc"]) / mw1
        rhoc2_val = float(d2["basic"]["rhoc"]) / mw2

        # Store the loaded JSON dicts directly (not Params -- they are
        # plain-Python nested dicts of floats/strings used to build Pyomo
        # EXPRESSIONS in r515b_pyomo_eos.py, not Pyomo components
        # themselves; there is nothing to make mutable/scalable about
        # them, they are fixed pure-component reference data).
        self._d1 = d1
        self._d2 = d2

        self.x1 = Param(
            initialize=x1_val, mutable=False,
            doc="Fixed pseudo-pure liquid=vapor mole fraction of R-1234ze(E) "
                "(x1=y1=z1 near-azeotropic simplification, spec rule-22 "
                "resolution 2026-08-18)",
            units=pyunits.dimensionless,
        )
        self.mw = Param(
            initialize=mw_mix_val, mutable=False,
            doc="Fixed pseudo-pure mixture molecular weight",
            units=pyunits.kg / pyunits.mol,
        )
        self.temperature_crit_1 = Param(initialize=tc1_val, mutable=False, units=pyunits.K)
        self.temperature_crit_2 = Param(initialize=tc2_val, mutable=False, units=pyunits.K)
        self.dens_mol_crit_1 = Param(initialize=rhoc1_val, mutable=False, units=pyunits.mol / pyunits.m**3)
        self.dens_mol_crit_2 = Param(initialize=rhoc2_val, mutable=False, units=pyunits.mol / pyunits.m**3)
        self.vol_mol_crit_1 = Param(initialize=1.0 / rhoc1_val, mutable=False, units=pyunits.m**3 / pyunits.mol)
        self.vol_mol_crit_2 = Param(initialize=1.0 / rhoc2_val, mutable=False, units=pyunits.m**3 / pyunits.mol)

    @classmethod
    def define_metadata(cls, obj):
        # `vapor_frac` is not among IDAES's StandardPropertySet names (the
        # closest standard name, `phase_frac`, is indexed by phase and does
        # not fit our single-Mix-phase/lever-rule-on-branches convention),
        # so it is declared via define_custom_properties() rather than
        # add_properties() -- using add_properties() for a non-standard name
        # still works but emits a deprecation warning (removal targeted for
        # 3.0.0) recommending exactly this split. All other properties below
        # are recognized standard names.
        obj.add_properties(
            {
                "flow_mass": {"method": None, "units": "kg/s"},
                "flow_mol": {"method": None, "units": "mol/s"},
                "pressure": {"method": None, "units": "Pa"},
                "enth_mass": {"method": None, "units": "J/kg"},
                "temperature": {"method": None, "units": "K"},
                "temperature_sat": {"method": None, "units": "K"},
                "dens_mass": {"method": None, "units": "kg/m^3"},
                "dens_mol": {"method": None, "units": "mol/m^3"},
                "entr_mass": {"method": None, "units": "J/kg/K"},
                "entr_mol": {"method": None, "units": "J/mol/K"},
                "enth_mol": {"method": None, "units": "J/mol"},
            }
        )
        obj.define_custom_properties(
            {
                "vapor_frac": {"method": None, "units": None},
            }
        )
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _R515BStateBlock(StateBlock):
    """
    Methods that apply to the R515B StateBlock as a whole (all time/space
    indices at once) rather than to a single R515BStateBlockData -- mirrors
    `saponification_thermo.py`'s `_StateBlock` pattern.
    """

    def fix_initialization_states(self):
        fix_state_vars(self)

    def initialize(
        blk,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine. Per the hard seeding requirement recorded
        in helmholtz_prop_validation.md Section 22 (2026-08-18 entry on
        the PH-flash system's seed-basin non-uniqueness), this method MUST
        seed every auxiliary Var (T_sat, rho_l_sat, rho_v_sat, vapor_frac,
        T_liq, rho_liq, T_vap, rho_vap) from the validated SciPy reference
        `r515b_helmholtz_core.flash_ph_pseudopure` -- NOT from generic
        defaults -- before any Pyomo/IPOPT solve is attempted. This is the
        rule-44 carve-out in action: the SciPy solve here is legitimate
        because it happens at initialize()-time, not inside the active
        NLP's Constraint bodies.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)
        else:
            flags = None
            for k in blk.values():
                if degrees_of_freedom(k) != 0:
                    raise RuntimeError(
                        "State vars fixed but degrees of freedom for state "
                        "block is not zero during initialization."
                    )

        d1 = blk[blk.index_set().first()].params._d1
        d2 = blk[blk.index_set().first()].params._d2
        for k in blk.values():
            x1_val = value(k.params.x1)
            mw_val = value(k.params.mw)
            p_val = value(k.pressure)
            h_molar_val = value(k.enth_mass) * mw_val

            rhoc1 = value(k.params.dens_mol_crit_1)
            rhoc2 = value(k.params.dens_mol_crit_2)
            rho_l_seed = 0.8 * (x1_val * rhoc1 + (1.0 - x1_val) * rhoc2)
            rho_v_seed = 0.01 * rho_l_seed
            t_seed = 0.5 * (T_MIN_K + T_MAX_K)

            ref = core.flash_ph_pseudopure(d1, d2, x1_val, p_val, h_molar_val, t_seed, rho_l_seed, rho_v_seed)
            if ref["status"] != "CONVERGED":
                init_log.warning(
                    f"{k.name}: pseudo-pure reference flash did not converge during "
                    f"initialization seeding; falling back to generic seed values."
                )
                t_sat_seed, rho_l_sat_seed, rho_v_sat_seed = t_seed, rho_l_seed, rho_v_seed
                vf_seed = 0.0
                t_liq_seed, rho_liq_seed = t_seed, rho_l_seed
                t_vap_seed, rho_vap_seed = t_seed, rho_v_seed
            else:
                t_sat_seed = ref["T_sat_K"]
                rho_l_sat_seed = ref["rho_l_sat_molm3"]
                rho_v_sat_seed = ref["rho_v_sat_molm3"]
                vf_seed = ref["vapor_frac"]
                if ref["region"] == "subcooled_liquid":
                    t_liq_seed, rho_liq_seed = ref["T_K"], ref["rho_molm3"]
                    t_vap_seed, rho_vap_seed = t_sat_seed, rho_v_sat_seed
                elif ref["region"] == "superheated_vapor":
                    t_liq_seed, rho_liq_seed = t_sat_seed, rho_l_sat_seed
                    t_vap_seed, rho_vap_seed = ref["T_K"], ref["rho_molm3"]
                else:
                    t_liq_seed, rho_liq_seed = t_sat_seed, rho_l_sat_seed
                    t_vap_seed, rho_vap_seed = t_sat_seed, rho_v_sat_seed

            k.T_sat.set_value(t_sat_seed)
            k.rho_l_sat.set_value(rho_l_sat_seed)
            k.rho_v_sat.set_value(rho_v_sat_seed)
            k.vapor_frac.set_value(vf_seed)
            k.T_liq.set_value(t_liq_seed)
            k.rho_liq.set_value(rho_liq_seed)
            k.T_vap.set_value(t_vap_seed)
            k.rho_vap.set_value(rho_vap_seed)

        if solver is None:
            from idaes.core.solvers import get_solver
            opt = get_solver(solver, optarg)
        else:
            from idaes.core.solvers import get_solver
            opt = get_solver(solver, optarg)

        for k in blk.values():
            res = opt.solve(k, tee=False)
            if str(res.solver.termination_condition) != "optimal":
                init_log.warning(f"{k.name}: initialization solve did not report optimal termination.")

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)
        init_log.info("Initialization Complete.")

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        if flags is None:
            return
        revert_state_vars(blk, flags)
        init_log.info("State Released.")


@declare_process_block_class("R515BStateBlock", block_class=_R515BStateBlock)
class R515BStateBlockData(StateBlockData):
    """
    StateBlockData for R-515B, PH state variables (flow_mass, pressure,
    enth_mass), MASS amount basis -- matching the exact contract
    `vapor_compression.py`'s generic Heater/Compressor/PressureChanger unit
    models expect from `HelmholtzParameterBlock`'s own PH+MASS mode (see
    module docstring for the confirmed interface details).

    Internally wires the 8-equation smooth PH-flash system from
    `r515b_pyomo_eos.mixture_ph_flash_residuals_expr` (Stage L part 3,
    validated 2026-08-18) as 8 always-active Constraints on 8 auxiliary
    Vars (T_sat/rho_l_sat/rho_v_sat, vapor_frac, T_liq/rho_liq, T_vap/
    rho_vap), then derives temperature/density/etc. from those via the
    same (1-vapor_frac)/vapor_frac blend proven correct in all three
    regions (subcooled liquid / two-phase / superheated vapor).
    """

    def build(self):
        super().build()
        params = self.params

        # ---- State variables (PH + MASS, matching the confirmed contract) ----
        self.flow_mass = Var(
            initialize=1.0, doc="Total mass flow rate", units=pyunits.kg / pyunits.s,
        )
        self.pressure = Var(
            initialize=5.0e5, bounds=(P_MIN_PA, P_MAX_PA), doc="State pressure",
            units=pyunits.Pa,
        )
        self.enth_mass = Var(
            initialize=3.0e5,
            bounds=(H_MOLAR_MIN_JMOL / value(params.mw), H_MOLAR_MAX_JMOL / value(params.mw)),
            doc="Total specific enthalpy (mass basis)", units=pyunits.J / pyunits.kg,
        )

        # ---- Auxiliary Vars for the smooth PH-flash system ----
        # DELIBERATELY UNITLESS (no `units=` kwarg -> Pyomo dimensionless):
        # `mixture_ph_flash_residuals_expr` and everything it calls
        # (`r515b_pyomo_eos.py`, ~30KB of Bell-2023 EOS/departure-function
        # algebra transcribed from the oracle) was built and validated
        # end-to-end (Stage L parts 1-3, `validate_pyomo_*_vs_core.py`)
        # entirely in raw-SI-value Pyomo Vars with NO `pyunits` anywhere in
        # that module -- confirmed by grep (zero `pyunits`/`units=` hits).
        # Attaching real units to these 8 Vars and feeding them straight
        # into that already-validated, units-naive expression tree is
        # exactly what triggered Stage M's first `DiagnosticsToolbox.
        # report_structural_issues()` finding ("Units problem with
        # expression..." on the saturation-pressure residual, since the
        # gas constant and every EOS coefficient inside are bare Python
        # floats with no unit annotation). Rather than retrofitting units
        # through the entire validated kernel (high risk of introducing a
        # NEW bug into already-machine-precision-matched math, for zero
        # physical benefit -- these are purely internal/auxiliary
        # quantities never exposed on the StateBlock's Port), the fix is
        # architectural: keep the internal 8-Var/8-Constraint system fully
        # unitless (by SI-value convention: T in K, rho in mol/m^3, exactly
        # as validated), and convert at the boundary where it talks to the
        # two unit-bearing PUBLIC state vars it needs (`pressure`, and the
        # molar enthalpy derived from `enth_mass`) by dividing by the
        # relevant unit below -- a standard, unit-system-safe way to
        # interface unit-naive legacy expressions with a unit-aware model
        # (dividing a units-bearing quantity by its own unit yields a
        # correctly-scaled dimensionless Pyomo expression, not a plain
        # Python float, so this remains exact and symbolic, not a
        # `value()` evaluation). Units are restored on every derived
        # Expression actually exposed on the Port (`temperature`,
        # `dens_mol`, `dens_mass`) below.
        self.T_sat = Var(initialize=300.0, bounds=(T_MIN_K, T_MAX_K), doc="Pseudo-pure saturation temperature at `pressure` [K, unitless Var by convention -- see note above]")
        self.rho_l_sat = Var(initialize=8000.0, bounds=(RHO_MIN_MOLM3, RHO_MAX_MOLM3), doc="Saturated-liquid molar density at `pressure` [mol/m^3, unitless Var by convention]")
        self.rho_v_sat = Var(initialize=500.0, bounds=(RHO_MIN_MOLM3, RHO_MAX_MOLM3), doc="Saturated-vapor molar density at `pressure` [mol/m^3, unitless Var by convention]")
        self.vapor_frac = Var(initialize=0.5, bounds=(-0.01, 1.01), doc="Vapor quality (mol/mol, == mass quality since composition is fixed)", units=pyunits.dimensionless)
        self.T_liq = Var(initialize=300.0, bounds=(T_MIN_K, T_MAX_K), doc="Liquid-branch temperature [K, unitless Var by convention] (real subcooled T when applicable; pinned to T_sat otherwise)")
        self.rho_liq = Var(initialize=8000.0, bounds=(RHO_MIN_MOLM3, RHO_MAX_MOLM3), doc="Liquid-branch molar density [mol/m^3, unitless Var by convention]")
        self.T_vap = Var(initialize=300.0, bounds=(T_MIN_K, T_MAX_K), doc="Vapor-branch temperature [K, unitless Var by convention] (real superheated T when applicable; pinned to T_sat otherwise)")
        self.rho_vap = Var(initialize=500.0, bounds=(RHO_MIN_MOLM3, RHO_MAX_MOLM3), doc="Vapor-branch molar density [mol/m^3, unitless Var by convention]")

        # ---- Derived molar quantities used to interface with the EOS kernel ----
        self.flow_mol = Expression(expr=self.flow_mass / params.mw, doc="Total mole flow rate")
        self._h_molar_jmol = Expression(
            expr=self.enth_mass * params.mw,
            doc="Total molar enthalpy [J/mol] -- internal, for the EOS kernel only "
                "(the public state var is enth_mass, mass basis, per the confirmed contract)",
        )
        # Unitless (SI-value) views of the two unit-bearing public
        # quantities the flash system needs, per the boundary-conversion
        # note above.
        self._p_pa_dimensionless = Expression(expr=self.pressure / pyunits.Pa)
        self._h_molar_jmol_dimensionless = Expression(
            expr=self._h_molar_jmol / (pyunits.J / pyunits.mol)
        )

        # ---- The 8-equation smooth PH-flash system (fully unitless, as validated) ----
        res = mixture_ph_flash_residuals_expr(
            params._d1, params._d2, params.x1,
            self._p_pa_dimensionless, self._h_molar_jmol_dimensionless,
            self.T_sat, self.rho_l_sat, self.rho_v_sat, self.vapor_frac,
            self.T_liq, self.rho_liq, self.T_vap, self.rho_vap,
            value(params.temperature_crit_1), value(params.temperature_crit_2),
            value(params.vol_mol_crit_1), value(params.vol_mol_crit_2),
        )
        self.eq_sat_p_liq = Constraint(expr=res["sat_p_liq"] == 0)
        self.eq_sat_p_vap = Constraint(expr=res["sat_p_vap"] == 0)
        self.eq_sat_gibbs = Constraint(expr=res["sat_gibbs"] == 0)
        self.eq_complementarity = Constraint(expr=res["complementarity"] == 0)
        self.eq_liq_p = Constraint(expr=res["liq_p"] == 0)
        self.eq_liq_h = Constraint(expr=res["liq_h"] == 0)
        self.eq_vap_p = Constraint(expr=res["vap_p"] == 0)
        self.eq_vap_h = Constraint(expr=res["vap_h"] == 0)

        # ---- Derived properties (Expressions, always valid in all 3 regions) ----
        # Units restored here -- these are the quantities actually exposed
        # on the StateBlock's Port (see define_port_members()/
        # define_display_vars() below), so they must be unit-bearing to
        # match vapor_compression.py's units-aware flowsheet arithmetic
        # (Tmin/Tmax Params/Constraints in K, P_low/P_high Vars in Pa,
        # etc. -- confirmed by direct inspection of vapor_compression.py).
        self.temperature = Expression(expr=res["T_actual"] * pyunits.K, doc="Actual state temperature")
        self._v_mol = Expression(
            expr=(1.0 - self.vapor_frac) / self.rho_liq + self.vapor_frac / self.rho_vap,
            doc="Actual state molar volume (lever rule on the liquid/vapor branches) [m^3/mol, unitless Expression by convention -- internal only]",
        )
        self.dens_mol = Expression(expr=(1.0 / self._v_mol) * (pyunits.mol / pyunits.m**3), doc="Actual state molar density")
        self.dens_mass = Expression(expr=self.dens_mol * params.mw, doc="Actual state mass density")

        # `temperature_sat` / `entr_mass`: added 2026-08-18 for Stage N
        # (vapor_compression.py integration) parity -- `vapor_compression.
        # py`'s evaporator/condenser superheat/subcool/approach constraints
        # and compressor vapor_constraint all reference `properties_out[0].
        # temperature_sat` (the general_helmholtz name), and `report_
        # solution()` plots entropy via `.entr_mass()`. Both are free
        # extensions of already-validated fields -- `temperature_sat` is
        # just `T_sat` (already solved by the 8-equation system above) with
        # units restored the same way as `temperature`; `entr_mass` uses
        # the entropy fields `mixture_ph_flash_residuals_expr` now also
        # returns (added same session, additive-only, no change to the 8
        # existing residuals -- see r515b_pyomo_eos.py). No new EOS math;
        # both are pure unit-restoring/mass-basis-converting wrappers
        # around already-validated symbolic-AD quantities (Section 24b).
        self.temperature_sat = Expression(expr=self.T_sat * pyunits.K, doc="Pseudo-pure saturation temperature at `pressure`")
        self._s_molar_jmolk = Expression(expr=res["S_actual_jmolk"] * (pyunits.J / pyunits.mol / pyunits.K))
        self.entr_mass = Expression(expr=self._s_molar_jmolk / params.mw, doc="Actual state specific entropy (mass basis)")
        # `entr_mol` (molar entropy): required by generic IDAES unit models
        # UNCONDITIONALLY, regardless of `amount_basis` -- e.g.
        # `PressureChanger.add_isentropic()` (used by `Compressor`) writes
        # its isentropic-assumption Constraint directly against `entr_mol`
        # on both the inlet StateBlock and an internal `properties_
        # isentropic` StateBlock (found the hard way: `Compressor(...)`
        # construction raised `PropertyNotSupportedError: ... entr_mol is
        # not supported`, Stage N's first construction attempt). Trivial
        # to support since it's the exact same already-validated molar
        # entropy Expression `entr_mass` is built from, just without the
        # `/params.mw` mass-basis conversion.
        self.entr_mol = Expression(expr=self._s_molar_jmolk, doc="Actual state molar entropy")
        # `enth_mol` (molar enthalpy): same rationale as `entr_mol` above --
        # `PressureChanger.init_isentropic()` (the Compressor's own
        # `initialize()` heuristic) accesses `properties_isentropic[t].
        # enth_mol` directly, unconditionally of `amount_basis` (found the
        # same way, one construction/initialize step later). Already have
        # the exact molar enthalpy this needs as `self._h_molar_jmol`
        # (used internally to drive the flash residuals) -- just exposing
        # it under the standard property name.
        self.enth_mol = Expression(expr=self._h_molar_jmol, doc="Actual state molar enthalpy")

    # ---- Interface methods required by generic IDAES unit models ----
    def get_material_flow_terms(self, p, j):
        return self.flow_mass

    def get_enthalpy_flow_terms(self, p):
        return self.flow_mass * self.enth_mass

    def get_material_density_terms(self, p, j):
        return self.dens_mass

    def get_energy_density_terms(self, p):
        return self.dens_mass * self.enth_mass

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "pressure": self.pressure,
            "enth_mass": self.enth_mass,
        }

    # define_port_members(): DELIBERATELY left at the StateBlockData base
    # class default (which returns define_state_vars() -- just flow_mass/
    # pressure/enth_mass), NOT overridden to also expose temperature/
    # vapor_frac. An earlier draft (this session, before Stage N exercised
    # it) DID override this to add temperature/vapor_frac, reasoning from
    # `StateBlock.build_port()`'s `Reference()` mechanism accepting both
    # Var and Expression port members at BUILD time. That is true, but
    # irrelevant to a different consumer: `idaes.core.util.initialization.
    # propagate_state()`, used by every unit model's own `initialize()`
    # method to copy values across an Arc, requires every port member to
    # be settable (a Var) -- it raises TypeError on an Expression member
    # (`temperature`). This was only caught here in Stage N (the first
    # point in this whole project where an actual Arc-connected flowsheet
    # calls `initialize()`/`propagate_state()`; Stage L/M's validation
    # never exercised it). Direct inspection of the REFERENCE package this
    # project must drop into confirms the fix: `HelmholtzStateBlockData.
    # define_port_members()` is NOT overridden at all -- it uses the base
    # class default, i.e. the port carries ONLY the 3 state vars, exactly
    # like ours now. `temperature`/`vapor_frac`/`entr_mass`/`entr_mol`/
    # `dens_mass`/`dens_mol`/`temperature_sat` remain fully accessible the
    # SAME way `vapor_compression.py` already accesses them for the
    # reference package -- directly off `control_volume.properties_out[0]`
    # /`properties_in[0]`, bypassing the Port entirely (confirmed: neither
    # `vapor_compression.py` nor `vapor_compression_r515b_integration.py`
    # ever reference `.outlet.temperature[0]`/`.inlet.vapor_frac[0]`
    # through a Port in PH mode).

    def define_display_vars(self):
        return {
            "Mass Flow": self.flow_mass,
            "Pressure": self.pressure,
            "Mass Enthalpy": self.enth_mass,
            "Temperature": self.temperature,
            "Saturation Temperature": self.temperature_sat,
            "Vapor Fraction": self.vapor_frac,
            "Mass Density": self.dens_mass,
            "Mass Entropy": self.entr_mass,
        }
