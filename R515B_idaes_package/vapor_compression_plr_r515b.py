"""
vapor_compression_plr_r515b.py -- PLR-enabled R-515B cycle, built the same
way the project's own R1234yf PLR benchmark (`R1234yf/vapor_compression_
plr_r1234yf.py`, found on the user's live `~/idaes-hvacr-cycles` device
folder, not the stale cloud snapshot) was built from `vapor_compression_
plr.py`: this is a copy of the already-validated Stage N `vapor_
compression_r515b_integration.py`, with the SAME small, additive PLR/CD
diff (`_compute_plf`, `fs.plr`/`fs.cd` Params, `get_full_load_cop()`/
`get_part_load_cop()`, plus the superheat/subcool floor+cap fix already
present in `vapor_compression_plr_fixed.py` and `vapor_compression_plr_
r1234yf.py`) grafted on, confirmed line-for-line against a diff of those
two real files rather than reimplemented from a description.

This exists to run the project's actual "PLR" COP-vs-ambient benchmark
(see `run_plr_only_cop_vs_ambient_with_r1234yf.py`'s setpoints: ambient
10-45C step 5, PLR=0.75, CD=0.13, evap_sat fixed at -29C, condenser
approach 9C, superheat/subcool 3C, compressor_efficiency 0.75) for
R-515B alongside R134a/R1234ze(E)/R1234yf, per the user's explicit
correction that the earlier ad hoc COP sweep (`cop_sweep_r515b_vs_
r134a.py`) used the wrong, self-chosen conditions instead of this
established benchmark.

One new addition beyond the straight PLR-diff graft: `_sat_at_t_robust`,
a small-step continuation-seeded wrapper around the already-validated
`core.solve_pseudopure_saturation_at_t` (Option 1 -- same equations, no
new methodology), because -29C = 244.15K sits below Stage L's originally
validated 255-375K saturation grid. Investigated directly before use
(see PROJECT_CONTEXT.md, 2026-08-19 entry): the raw solver is fragile
across ~240-255K (its default multiplicative retry seeds have a narrow
basin and occasionally land on a spurious collapsed root, correctly
rejected by the existing separation gate), but the TRUE physical branch
converges reliably and consistently when seeded from the nearest
previously-converged neighbor via careful small-step continuation from
the validated 255K anchor. Confirmed exact value at -29C this way:
P_sat=66,937.17 Pa, rho_liq=11,186.59 mol/m^3, rho_vap=34.06 mol/m^3.

--- Original Stage N docstring below, still accurate for everything not
called out above ---

vapor_compression_r515b_integration.py -- Stage N: a NEW integration copy
of the vapor-compression-cycle flowsheet (`vapor_compression.py`, which
remains read-only and UNMODIFIED per the master task's hard constraint),
substituting the validated `R515BParameterBlock` (Stage L, this project)
for `HelmholtzParameterBlock`.

Why this is a new, self-contained file rather than a thin subclass/import
of `vapor_compression.py`: the master task requires `vapor_compression.py`
to remain strictly read-only, and "a new integration copy" was the
explicitly scoped Stage N deliverable -- so this file duplicates and
adapts the relevant logic rather than importing internals from that file.

Scope: `vapor_compression.py`'s `Mode.PH` branch is the ONLY mode
`R515BParameterBlock` can drop into -- it implements PH+MASS state vars
exclusively (no TPX state-var support was built; this was the confirmed
contract researched earlier this session by direct inspection of
`vapor_compression.py`'s own Mode.PH usage and `general_helmholtz.
helmholtz_state`'s PH+MASS `_state_vars()`). This file therefore hardcodes
PH-mode behavior throughout and omits the ORIGINAL_TPX/IMPROVED_TPX
branches entirely (dead code that could never execute against this
property package would be worse than no code).

Two further adaptations, both because R-515B is a fixed-composition BLEND
(not a literal CoolProp pure fluid), and both explicitly Option-1/Option-2
compliant per the 2026-08-18 standing commitment (helmholtz_prop_
validation.md Section 22: no fitted R-515B ancillary correlations without
pausing for approval):

1. `vapor_compression.py`'s `specify_initial_conditions()` calls
   `CP.PropsSI(...)` against `self.fluid_name` for initial-guess
   pressures/enthalpies. Even if CoolProp happens to carry an R-515B
   model, it would be a SEPARATE, unvalidated implementation from this
   project's oracle-traced one -- using it here would silently
   reintroduce exactly the "different model" risk this whole project
   exists to avoid. Replaced below with calls into this project's own
   validated `r515b_helmholtz_core.solve_pseudopure_saturation_at_t`
   (Option 1: the real VLE solve, never a fitted correlation).
2. `draw_thermodynamic_diagrams()`/`report_solution()` call
   `self.model.fs.properties.hp_diagram()` etc. -- convenience plotting
   methods specific to `HelmholtzParameterBlock`'s compiled ancillary-
   equation background, not part of the generic `PhysicalParameterBlock`
   interface at all (any custom property package integrated this way
   would need the same graceful substitution). Replaced with an
   equivalent envelope built from repeated REAL saturation solves, swept
   over the SAME 255-375K/5K grid already validated in `validate_pyomo_
   saturation_vs_core.py` (no new, unvalidated T-range explored) -- again
   Option 1, not a fitted correlation.

Every other method (`_define_flowsheet`, `initialize`, `set_specifications`,
`optimize_COP`, `report_solution`) is a faithful line-for-line port of
`vapor_compression.py`'s PH-mode code path, since that is exactly the
already-confirmed drop-in contract: generic property-package-agnostic
`Heater`/`Compressor`/`PressureChanger` unit models via `control_volume`,
PH+MASS state vars, and `temperature`/`temperature_sat`/`vapor_frac`/
`entr_mass` all present on `R515BParameterBlock` (Stage L, plus this
session's Stage-N-prep additions to `r515b_property_package.py`).
"""
import sys
import logging
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from idaes.core import FlowsheetBlock  # noqa: E402
from idaes.models.unit_models import Heater, Compressor, PressureChanger  # noqa: E402
from idaes.core.solvers import get_solver  # noqa: E402
from idaes.core.util.scaling import calculate_scaling_factors  # noqa: E402
from idaes.core.util.initialization import propagate_state  # noqa: E402
from idaes.core.util import DiagnosticsToolbox  # noqa: E402
from pyomo.environ import ConcreteModel, value, Param, Var, TransformationFactory, maximize  # noqa: E402
from pyomo.network import Arc  # noqa: E402
from pyomo.environ import units as pyunits  # noqa: E402

import r515b_helmholtz_core as core  # noqa: E402
from r515b_property_package import R515BParameterBlock  # noqa: E402

C_to_K = 273.15

# Same continuation grid already validated end-to-end in
# `validate_pyomo_saturation_vs_core.py` (25/25 PASS at machine epsilon) --
# reused here for the envelope sweep so no NEW, unvalidated temperature
# range is explored just to draw a picture.
_ENVELOPE_T_GRID_K = [255.0 + 5.0 * i for i in range(25)]  # 255..375 K


class R515BVaporCompressionCyclePLR:
    """PLR-enabled variant of `R515BVaporCompressionCycle` (Stage N),
    mirroring `SimpleVaporCompressionCyclePLRR1234yf`'s relationship to
    `SimpleVaporCompressionCyclePLR`. Simple vapor-compression cycle
    (evaporator -> compressor -> condenser -> expansion valve), PH state
    vars, R-515B pseudo-pure property package."""

    def __init__(self, compressor_efficiency=0.75, PLR=0.5, CD=0.25):
        """
        Parameters:
            compressor_efficiency : float
                Isentropic efficiency of the compressor.
            PLR : float
                Part-load ratio in [0, 1] (metadata only; does not alter
                the thermodynamic equations -- same as the original PLR
                files' own design).
            CD : float
                Degradation coefficient in [0, 1] used by PLF = 1 - CD*(1-PLR).
        """
        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        self.model.fs.properties = R515BParameterBlock()

        assert 0 < compressor_efficiency < 1, "Compressor efficiency must be between 0 and 1"
        self.compressor_efficiency = compressor_efficiency
        assert 0.0 <= PLR <= 1.0, "PLR must be in [0,1]"
        self.plr = PLR
        assert 0.0 <= CD <= 1.0, "CD must be in [0, 1]"
        self.cd = CD

        self.optimization_converged = None

        # Reference data (mole fraction, molecular weight, saturation
        # envelope) used by `specify_initial_conditions`/diagram helpers --
        # all sourced from the already-validated SciPy core, never fitted.
        self._d1 = core.load_idaes_helmholtz_json(core.FLUID1)
        self._d2 = core.load_idaes_helmholtz_json(core.FLUID2)
        self._x1 = core.r515b_x1()
        self._mw1 = core.mw_from_json(self._d1)
        self._mw2 = core.mw_from_json(self._d2)
        self._mw_mix = self._x1 * self._mw1 + (1.0 - self._x1) * self._mw2
        rhoc1 = float(self._d1["basic"]["rhoc"]) / self._mw1
        rhoc2 = float(self._d2["basic"]["rhoc"]) / self._mw2
        self._rho_l_seed0 = 0.8 * (self._x1 * rhoc1 + (1.0 - self._x1) * rhoc2)
        self._rho_v_seed0 = 0.01 * self._rho_l_seed0
        self._envelope = None  # lazily built by _compute_envelope()
        self._cold_seed_cache = {}  # t_k (rounded) -> (rho_l, rho_v), see _sat_at_t_robust

        self._define_flowsheet()

    @staticmethod
    def _compute_plf(plr: float, cd: float) -> float:
        """Compute part-load factor using PLF = 1 - CD * (1 - PLR).
        Identical formula/signature to `SimpleVaporCompressionCyclePLR.
        _compute_plf` in `vapor_compression_plr.py` (confirmed via diff)."""
        if not (0.0 <= plr <= 1.0):
            raise ValueError("PLR must be in [0, 1]")
        if not (0.0 <= cd <= 1.0):
            raise ValueError("CD must be in [0, 1]")
        plf = 1.0 - cd * (1.0 - plr)
        return max(0.0, min(1.0, plf))

    # ------------------------------------------------------------------
    # Own validated-solver helpers (Option 1: real VLE solves, never a
    # fitted correlation -- see module docstring)
    # ------------------------------------------------------------------
    def _sat_at_t(self, t_k, rho_l_seed=None, rho_v_seed=None):
        """Real saturation solve at fixed T via the validated SciPy core."""
        rl = self._rho_l_seed0 if rho_l_seed is None else rho_l_seed
        rv = self._rho_v_seed0 if rho_v_seed is None else rho_v_seed
        return core.solve_pseudopure_saturation_at_t(self._d1, self._d2, self._x1, t_k, rl, rv)

    def _sat_at_t_robust(self, t_k, anchor_t_k=255.0, step_k=1.0):
        """Cold-continuation wrapper for T below the Stage L validated
        255-375K grid (see module docstring for the investigation behind
        this). Walks down from the validated `anchor_t_k` in `step_k`
        steps, always calling the SAME already-validated `core.solve_
        pseudopure_saturation_at_t` (no new equations, no manual seed
        grid), carrying forward only successfully-converged (rho_l,
        rho_v) as the seed for the next step and simply skipping
        (not fitting/faking) isolated divergent steps -- then solves
        exactly at the requested `t_k` using the last good seed. For
        t_k >= 255K this is equivalent to (and slightly more expensive
        than) calling `_sat_at_t` directly, so it's only used for the
        below-grid case. Results are cached per rounded target T since
        the PLR benchmark holds evap_sat fixed across its whole ambient
        sweep."""
        key = round(float(t_k), 3)
        if key in self._cold_seed_cache:
            rl, rv = self._cold_seed_cache[key]
            return core.solve_pseudopure_saturation_at_t(self._d1, self._d2, self._x1, t_k, rl, rv)

        if t_k >= anchor_t_k:
            return self._sat_at_t(t_k)

        rl, rv = self._rho_l_seed0, self._rho_v_seed0
        row_anchor = self._sat_at_t(anchor_t_k, rl, rv)
        if row_anchor["status"] != "CONVERGED":
            raise RuntimeError(f"Cold-continuation anchor at {anchor_t_k}K did not converge")
        rl, rv = row_anchor["rho_l_molm3"], row_anchor["rho_v_molm3"]

        n_steps = int(round((anchor_t_k - t_k) / step_k))
        for i in range(1, n_steps + 1):
            t_step = anchor_t_k - i * step_k
            row = core.solve_pseudopure_saturation_at_t(self._d1, self._d2, self._x1, t_step, rl, rv)
            if row["status"] == "CONVERGED":
                rl, rv = row["rho_l_molm3"], row["rho_v_molm3"]
            # else: isolated divergent step, skip -- keep the last good seed

        final = core.solve_pseudopure_saturation_at_t(self._d1, self._d2, self._x1, t_k, rl, rv)
        if final["status"] == "CONVERGED":
            self._cold_seed_cache[key] = (final["rho_l_molm3"], final["rho_v_molm3"])
        return final

    def _compute_envelope(self):
        """Sweeps the same validated 255-375K/5K grid used in
        `validate_pyomo_saturation_vs_core.py`, with continuation seeding
        (each point seeds the next), to build a real saturation-envelope
        dataset for plotting -- P_Pa, h_l/h_v [J/kg], s_l/s_v [J/kg/K],
        T [K] arrays. Cached after first call."""
        if self._envelope is not None:
            return self._envelope

        rl, rv = self._rho_l_seed0, self._rho_v_seed0
        T, P, Hl, Hv, Sl, Sv = [], [], [], [], [], []
        for t_k in _ENVELOPE_T_GRID_K:
            row = self._sat_at_t(t_k, rl, rv)
            if row["status"] != "CONVERGED":
                continue
            rl, rv = row["rho_l_molm3"], row["rho_v_molm3"]
            s_l = core.mix_entropy_direct(self._d1, self._d2, t_k, rl, self._x1)
            s_v = core.mix_entropy_direct(self._d1, self._d2, t_k, rv, self._x1)
            T.append(t_k)
            P.append(row["P_Pa"])
            Hl.append(row["h_l_Jmol"] / self._mw_mix)
            Hv.append(row["h_v_Jmol"] / self._mw_mix)
            Sl.append(s_l / self._mw_mix)
            Sv.append(s_v / self._mw_mix)

        self._envelope = {
            "T_K": np.array(T), "P_Pa": np.array(P),
            "h_l_Jkg": np.array(Hl), "h_v_Jkg": np.array(Hv),
            "s_l_JkgK": np.array(Sl), "s_v_JkgK": np.array(Sv),
        }
        return self._envelope

    def draw_thermodynamic_diagrams(self):
        """Real-envelope P-h / T-P / T-s diagrams (no unit model state
        required) -- built from repeated validated VLE solves rather than
        `HelmholtzParameterBlock`'s compiled ancillary-correlation
        `hp_diagram`/`pt_diagram`/`ts_diagram` (see module docstring)."""
        env = self._compute_envelope()

        plt.figure()
        plt.plot(env["h_l_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "b-", label="bubble (liquid)")
        plt.plot(env["h_v_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "r-", label="dew (vapor)")
        plt.xlabel("h [kJ/kg]"); plt.ylabel("P [kPa]"); plt.yscale("log")
        plt.title("R-515B P-h envelope (real VLE solves)"); plt.legend()
        plt.show()

        plt.figure()
        plt.plot(env["T_K"], env["P_Pa"] / 1000.0, "k-")
        plt.xlabel("T [K]"); plt.ylabel("P [kPa]"); plt.yscale("log")
        plt.title("R-515B saturation P-T curve"); plt.show()

        plt.figure()
        plt.plot(env["s_l_JkgK"] / 1000.0, env["T_K"], "b-", label="bubble (liquid)")
        plt.plot(env["s_v_JkgK"] / 1000.0, env["T_K"], "r-", label="dew (vapor)")
        plt.xlabel("s [kJ/kg/K]"); plt.ylabel("T [K]")
        plt.title("R-515B T-s envelope (real VLE solves)"); plt.legend()
        plt.show()

    def _define_flowsheet(self):
        """Faithful port of `vapor_compression.py`'s `_define_flowsheet`
        (Mode.PH code path only)."""
        logging.basicConfig(level=logging.WARNING)
        self.logger = logging.getLogger(__name__)

        self.model.fs.evaporator = Heater(property_package=self.model.fs.properties)
        self.model.fs.compressor = Compressor(property_package=self.model.fs.properties)
        self.model.fs.condenser = Heater(property_package=self.model.fs.properties)
        self.model.fs.expansion_valve = PressureChanger(
            property_package=self.model.fs.properties,
            thermodynamic_assumption="adiabatic",
            compressor=False,
        )

        self.model.fs.evaporator_to_compressor = Arc(
            source=self.model.fs.evaporator.outlet, destination=self.model.fs.compressor.inlet)
        self.model.fs.compressor_to_condenser = Arc(
            source=self.model.fs.compressor.outlet, destination=self.model.fs.condenser.inlet)
        self.model.fs.condenser_to_expansion_valve = Arc(
            source=self.model.fs.condenser.outlet, destination=self.model.fs.expansion_valve.inlet)
        self.model.fs.expansion_valve_to_evaporator = Arc(
            source=self.model.fs.expansion_valve.outlet, destination=self.model.fs.evaporator.inlet)

        TransformationFactory('network.expand_arcs').apply_to(self.model)

        self.model.fs.evaporator_to_compressor_expanded.flow_mass_equality.deactivate()

        self.model.fs.cop = Var(initialize=1, units=pyunits.dimensionless, bounds=(0.1, 100))

        @self.model.fs.Constraint(doc="COP constraint")
        def compute_cop(b):
            return b.cop * b.compressor.work_mechanical[0] == b.evaporator.heat_duty[0]

        @self.model.fs.Objective(doc="Maximize COP", sense=maximize)
        def obj(b):
            return b.cop

        self.model.fs.compute_cop.deactivate()
        self.model.fs.obj.deactivate()

        # PLR metadata only; does not alter thermodynamic equations.
        self.model.fs.plr = Param(initialize=self.plr, units=pyunits.dimensionless, mutable=True)
        self.model.fs.cd = Param(initialize=self.cd, units=pyunits.dimensionless, mutable=True)

        self.model.fs.evaporator.superheating = Param(initialize=0, units=pyunits.K, mutable=True)
        self.model.fs.evaporator.superheating_max = Param(initialize=8, units=pyunits.K, mutable=True)
        self.model.fs.evaporator.T_sat_set = Param(initialize=C_to_K, units=pyunits.K, mutable=True)

        @self.model.fs.evaporator.Constraint(doc="Superheat evaporator outlet")
        def superheating_constraint(b):
            return b.control_volume.properties_out[0].temperature >= \
                b.control_volume.properties_out[0].temperature_sat + b.superheating

        self.model.fs.evaporator.superheating_constraint.deactivate()

        @self.model.fs.evaporator.Constraint(doc="Superheat evaporator outlet upper cap")
        def superheating_upper_constraint(b):
            return b.control_volume.properties_out[0].temperature <= \
                b.control_volume.properties_out[0].temperature_sat + b.superheating_max

        self.model.fs.evaporator.superheating_upper_constraint.deactivate()

        @self.model.fs.evaporator.Constraint(doc="Evaporator saturation temperature setpoint")
        def evap_sat_constraint(b):
            return b.control_volume.properties_out[0].temperature_sat == b.T_sat_set

        self.model.fs.evaporator.evap_sat_constraint.deactivate()

        self.model.fs.condenser.subcooling = Param(initialize=0, units=pyunits.K, mutable=True)
        self.model.fs.condenser.subcooling_max = Param(initialize=8, units=pyunits.K, mutable=True)
        self.model.fs.condenser.ambient_T = Param(initialize=C_to_K, units=pyunits.K, mutable=True)
        self.model.fs.condenser.approach_T = Param(initialize=0, units=pyunits.K, mutable=True)

        @self.model.fs.condenser.Constraint(doc="Subcool condenser outlet")
        def subcooling_constraint(b):
            return b.control_volume.properties_out[0].temperature <= \
                b.control_volume.properties_out[0].temperature_sat - b.subcooling

        self.model.fs.condenser.subcooling_constraint.deactivate()

        @self.model.fs.condenser.Constraint(doc="Subcool condenser outlet lower cap")
        def subcooling_lower_constraint(b):
            return b.control_volume.properties_out[0].temperature >= \
                b.control_volume.properties_out[0].temperature_sat - b.subcooling_max

        self.model.fs.condenser.subcooling_lower_constraint.deactivate()

        @self.model.fs.condenser.Constraint(doc="Condenser saturation temperature setpoint (sat = ambient + approach)")
        def approach_constraint(b):
            return b.control_volume.properties_out[0].temperature_sat == b.ambient_T + b.approach_T

        self.model.fs.condenser.approach_constraint.deactivate()

        @self.model.fs.compressor.Constraint(doc="Must be a vapor")
        def vapor_constraint(b):
            return b.control_volume.properties_out[0].temperature >= \
                b.control_volume.properties_out[0].temperature_sat

        self.model.fs.compressor.vapor_constraint.deactivate()

        # Explicit high/low pressure levels (used when arc pressure equalities are disabled)
        self.model.fs.P_low = Var(initialize=200e3, units=pyunits.Pa, bounds=(50e3, 2000e3))
        self.model.fs.P_high = Var(initialize=1200e3, units=pyunits.Pa, bounds=(100e3, 5000e3))
        self.model.fs.P_low_target = Param(initialize=200e3, units=pyunits.Pa, mutable=True)
        self.model.fs.P_high_target = Param(initialize=1200e3, units=pyunits.Pa, mutable=True)

        @self.model.fs.Constraint(doc="Low-side pressure target")
        def P_low_target_constraint(b):
            return b.P_low == b.P_low_target

        @self.model.fs.Constraint(doc="High-side pressure target")
        def P_high_target_constraint(b):
            return b.P_high == b.P_high_target

        @self.model.fs.Constraint(doc="Evaporator inlet pressure equals P_low")
        def P_low_evap_in(b):
            return b.evaporator.inlet.pressure[0] == b.P_low

        @self.model.fs.Constraint(doc="Evaporator outlet pressure equals P_low")
        def P_low_evap_out(b):
            return b.evaporator.outlet.pressure[0] == b.P_low

        @self.model.fs.Constraint(doc="Compressor inlet pressure equals P_low")
        def P_low_comp_in(b):
            return b.compressor.inlet.pressure[0] == b.P_low

        @self.model.fs.Constraint(doc="Valve outlet pressure equals P_low")
        def P_low_valve_out(b):
            return b.expansion_valve.outlet.pressure[0] == b.P_low

        @self.model.fs.Constraint(doc="Compressor outlet pressure equals P_high")
        def P_high_comp_out(b):
            return b.compressor.outlet.pressure[0] == b.P_high

        @self.model.fs.Constraint(doc="Condenser inlet pressure equals P_high")
        def P_high_cond_in(b):
            return b.condenser.inlet.pressure[0] == b.P_high

        @self.model.fs.Constraint(doc="Condenser outlet pressure equals P_high")
        def P_high_cond_out(b):
            return b.condenser.outlet.pressure[0] == b.P_high

        @self.model.fs.Constraint(doc="Valve inlet pressure equals P_high")
        def P_high_valve_in(b):
            return b.expansion_valve.inlet.pressure[0] == b.P_high

        for c in [self.model.fs.P_low_target_constraint, self.model.fs.P_high_target_constraint,
                  self.model.fs.P_low_evap_in, self.model.fs.P_low_evap_out,
                  self.model.fs.P_low_comp_in, self.model.fs.P_low_valve_out,
                  self.model.fs.P_high_comp_out, self.model.fs.P_high_cond_in,
                  self.model.fs.P_high_cond_out, self.model.fs.P_high_valve_in]:
            c.deactivate()

        units = [self.model.fs.evaporator, self.model.fs.compressor,
                  self.model.fs.condenser, self.model.fs.expansion_valve]

        for u in units:
            u.Tmin = Param(initialize=C_to_K, mutable=True)
            u.Tmax = Param(initialize=C_to_K + 10, mutable=True)

            @u.Constraint(doc="Temperature lower bound")
            def T_lower_bound(b):
                return b.control_volume.properties_out[0].temperature >= u.Tmin

            u.T_lower_bound.deactivate()

            @u.Constraint(doc="Temperature upper bound")
            def T_upper_bound(b):
                return b.control_volume.properties_out[0].temperature <= u.Tmax

            u.T_upper_bound.deactivate()

    def specify_initial_conditions(self, low_side_temperature=-20, high_side_temperature=30):
        """Port of `vapor_compression.py`'s `specify_initial_conditions`,
        with `CP.PropsSI` replaced by this project's own validated
        saturation solver (see module docstring, adaptation 1)."""
        low_side_temperature += C_to_K
        high_side_temperature += C_to_K
        superheat = 3
        subcool = 3

        # Use the cold-robust wrapper throughout -- identical to _sat_at_t
        # for T>=255K (Stage L's validated grid), and correctly handles
        # low_side_temperature down to the PLR benchmark's -29C target
        # via careful continuation (see module docstring/_sat_at_t_robust).
        sat_low = self._sat_at_t_robust(low_side_temperature)
        sat_high = self._sat_at_t_robust(high_side_temperature)
        sat_low_sh = self._sat_at_t_robust(low_side_temperature + superheat)
        sat_high_sh = self._sat_at_t_robust(high_side_temperature + superheat)
        sat_high_sc = self._sat_at_t_robust(high_side_temperature - subcool)
        for row, name in [(sat_low, "low_side"), (sat_high, "high_side"),
                           (sat_low_sh, "low_side+superheat"),
                           (sat_high_sh, "high_side+superheat"),
                           (sat_high_sc, "high_side-subcool")]:
            if row["status"] != "CONVERGED":
                raise RuntimeError(f"Saturation solve for initial conditions did not converge at {name}")

        low_side_pressure = sat_low["P_Pa"]
        high_side_pressure = sat_high["P_Pa"]

        # Expansion valve outlet (assume slightly vaporized, Q=0.2 lever rule)
        low_side_liquid_H = (sat_low["h_l_Jmol"] + 0.2 * (sat_low["h_v_Jmol"] - sat_low["h_l_Jmol"])) / self._mw_mix

        # Evaporator outlet (saturated vapor at T+superheat)
        low_side_vapor_H = sat_low_sh["h_v_Jmol"] / self._mw_mix

        # Compressor outlet seed (saturated vapor at high_side_T - subcool)
        high_side_vapor_H = sat_high_sc["h_v_Jmol"] / self._mw_mix

        # Condenser outlet seed (saturated liquid at high_side_T + superheat)
        high_side_liquid_H = sat_high_sh["h_l_Jmol"] / self._mw_mix

        self.h_init = np.array([low_side_vapor_H, high_side_vapor_H, high_side_liquid_H, low_side_liquid_H])
        self.p_init = np.array([low_side_pressure, high_side_pressure, high_side_pressure, low_side_pressure])
        self.T_init = np.array([
            low_side_temperature + superheat,
            high_side_temperature + superheat,
            high_side_temperature - subcool,
            low_side_temperature,
        ])

        env = self._compute_envelope()
        plt.figure()
        plt.plot(env["h_l_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "b-")
        plt.plot(env["h_v_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "r-")
        plt.plot(self.h_init / 1000.0, self.p_init / 1000.0, "ko")
        plt.xlabel("h [kJ/kg]"); plt.ylabel("P [kPa]"); plt.yscale("log")
        plt.title("Initial conditions on R-515B P-h envelope"); plt.show()

    def initialize(self, verbose=False):
        """Faithful port of `vapor_compression.py`'s `initialize` (Mode.PH branch)."""
        p_scale = 1

        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)
        self.model.fs.evaporator.inlet.pressure[0].fix(self.p_init[-1] * p_scale)
        self.model.fs.evaporator.inlet.enth_mass[0].fix(self.h_init[-1])
        self.model.fs.evaporator.outlet.enth_mass[0].fix(self.h_init[0])

        self.logger.info("Initializing evaporator...")
        self.model.fs.evaporator.initialize(outlvl=logging.WARNING)
        if verbose:
            self.model.fs.evaporator.report()
        propagate_state(self.model.fs.evaporator_to_compressor)

        self.model.fs.compressor.inlet.pressure[0].fix(self.p_init[0] * p_scale)
        self.model.fs.compressor.inlet.enth_mass[0].fix(self.h_init[0])
        self.model.fs.compressor.outlet.pressure[0].fix(self.p_init[1] * p_scale)
        self.model.fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)

        self.logger.info("Initializing compressor...")
        self.model.fs.compressor.initialize(outlvl=logging.WARNING)
        if verbose:
            self.model.fs.compressor.report()
        propagate_state(self.model.fs.compressor_to_condenser)

        self.model.fs.condenser.inlet.pressure[0].fix(self.p_init[1] * p_scale)
        self.model.fs.condenser.inlet.enth_mass[0].fix(self.h_init[1])
        self.model.fs.condenser.outlet.enth_mass[0].fix(self.h_init[2])

        self.logger.info("Initializing condenser...")
        self.model.fs.condenser.initialize(outlvl=logging.WARNING)
        if verbose:
            self.model.fs.condenser.report()
        propagate_state(self.model.fs.condenser_to_expansion_valve)

        self.model.fs.expansion_valve.inlet.pressure[0].fix(self.p_init[2] * p_scale)
        self.model.fs.expansion_valve.outlet.pressure[0].fix(self.p_init[3] * p_scale)

        self.logger.info("Initializing expansion valve...")
        self.model.fs.expansion_valve.initialize(outlvl=logging.WARNING)
        if verbose:
            self.model.fs.expansion_valve.report()
        propagate_state(self.model.fs.expansion_valve_to_evaporator)

        if verbose:
            print("\nFinished initialization. Stream summary:")
            self.model.fs.report()

    def set_specifications(self,
                            low_side_pressure=(200, 500),
                            high_side_pressure=(1000, 3000),
                            evaporator_temperature=(-20, 0),
                            compressor_temperature=None,
                            condenser_temperature=(30, 50),
                            expansion_valve_temperature=None,
                            subcooling=3,
                            superheating=3,
                            max_pressure_ratio=4,
                            ambient_temperature=None,
                            condenser_approach=None,
                            evap_sat_temperature=None,
                            debug_disable_arc_pressure_eq=False,
                            plr=None,
                            cd=None,
                            superheating_max=8,
                            subcooling_max=8):
        """Faithful port of `vapor_compression.py`'s `set_specifications`
        (Mode.PH branch only -- TPX-only vapor_frac bound/fix branches
        removed, since R515BParameterBlock has no TPX state vars).

        ONE deliberate, independent deviation from the original (found via
        an end-to-end test optimization that ran away to an unphysical
        COP~34.5 by superheating the evaporator outlet to 412 K): the
        original file guards each temperature-bound element with a bare
        `if evaporator_temperature[0]:`-style truthiness check. In Python,
        `0` is falsy, so a caller who legitimately wants an upper (or
        lower) bound of exactly 0 degC -- a perfectly ordinary refrigerant-
        cycle setpoint -- gets that bound SILENTLY skipped, leaving the
        corresponding T_lower_bound/T_upper_bound Constraint permanently
        deactivated with no error or warning. This is a pre-existing
        quirk in the read-only `vapor_compression.py` (never modified
        there, per the master task's hard constraint), reproduced
        faithfully at first in this file too -- but since this integration
        file is new/independent code, not a constraint on the original,
        fixing it here is a legitimate, low-risk, clearly-documented
        improvement: every `if X[i]:` below is `if X[i] is not None:`
        instead, so an explicit 0 degC bound is honored like any other
        value."""
        assert superheating >= 0, "Superheating must be greater than or equal to 0"
        assert subcooling >= 0, "Subcooling must be greater than or equal to 0"
        assert max_pressure_ratio > 1.2, "Maximum pressure ratio must be greater than 1.2"
        if plr is not None:
            assert 0.0 <= plr <= 1.0, "PLR must be in [0, 1]"
            self.plr = plr
        if cd is not None:
            assert 0.0 <= cd <= 1.0, "CD must be in [0, 1]"
            self.cd = cd
        self.model.fs.plr.set_value(self.plr)
        self.model.fs.cd.set_value(self.cd)

        self.unit_operations = [self.model.fs.evaporator, self.model.fs.compressor,
                                 self.model.fs.condenser, self.model.fs.expansion_valve]

        for unit in self.unit_operations:
            unit.inlet.flow_mass[0].unfix()
            unit.outlet.flow_mass[0].unfix()
            unit.inlet.enth_mass[0].unfix()
            unit.outlet.enth_mass[0].unfix()
            unit.inlet.pressure[0].unfix()
            unit.outlet.pressure[0].unfix()
            unit.T_lower_bound.deactivate()
            unit.T_upper_bound.deactivate()

        def check_input(bounds):
            return bounds is not None and len(bounds) == 2

        evap_sat_active = evap_sat_temperature is not None
        cond_sat_active = (ambient_temperature is not None) and (condenser_approach is not None)

        if evap_sat_active or cond_sat_active:
            low_side_pressure_min, low_side_pressure_max = 50 * 1000, 2000 * 1000
            high_side_pressure_min, high_side_pressure_max = 100 * 1000, 5000 * 1000
        else:
            if check_input(low_side_pressure):
                low_side_pressure_min, low_side_pressure_max = low_side_pressure[0] * 1000, low_side_pressure[1] * 1000
            else:
                low_side_pressure_min = low_side_pressure_max = None
            if check_input(high_side_pressure):
                high_side_pressure_min, high_side_pressure_max = high_side_pressure[0] * 1000, high_side_pressure[1] * 1000
            else:
                high_side_pressure_min = high_side_pressure_max = None

        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)

        if low_side_pressure_min:
            self.model.fs.evaporator.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.evaporator.inlet.pressure[0].setub(low_side_pressure_max)

        if check_input(evaporator_temperature):
            if evaporator_temperature[0] is not None:
                self.model.fs.evaporator.Tmin.set_value(evaporator_temperature[0] + C_to_K)
                self.model.fs.evaporator.T_lower_bound.activate()
            if evaporator_temperature[1] is not None:
                self.model.fs.evaporator.Tmax.set_value(evaporator_temperature[1] + C_to_K)
                self.model.fs.evaporator.T_upper_bound.activate()

        if superheating > 0.1:
            self.model.fs.evaporator.superheating.set_value(superheating)
            self.model.fs.evaporator.superheating_constraint.activate()
            self.model.fs.evaporator.superheating_max.set_value(superheating_max)
            self.model.fs.evaporator.superheating_upper_constraint.activate()
        else:
            self.model.fs.evaporator.superheating_constraint.deactivate()
            self.model.fs.evaporator.superheating_upper_constraint.deactivate()

        if evap_sat_temperature is not None:
            self.model.fs.evaporator.T_sat_set.set_value(evap_sat_temperature + C_to_K)
            self.model.fs.evaporator.evap_sat_constraint.activate()
        else:
            self.model.fs.evaporator.evap_sat_constraint.deactivate()

        if low_side_pressure_min:
            self.model.fs.compressor.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.compressor.inlet.pressure[0].setub(low_side_pressure_max)
        if high_side_pressure_min:
            self.model.fs.compressor.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.compressor.outlet.pressure[0].setub(high_side_pressure_max)

        if check_input(compressor_temperature):
            if compressor_temperature[0] is not None:
                self.model.fs.compressor.Tmin.set_value(compressor_temperature[0] + C_to_K)
                self.model.fs.compressor.T_lower_bound.activate()
            if compressor_temperature[1] is not None:
                self.model.fs.compressor.Tmax.set_value(compressor_temperature[1] + C_to_K)
                self.model.fs.compressor.T_upper_bound.activate()

        self.model.fs.compressor.ratioP.setub(max_pressure_ratio)
        self.model.fs.compressor.ratioP.setlb(1.1)

        if high_side_pressure_min:
            self.model.fs.condenser.inlet.pressure[0].setlb(high_side_pressure_min)
            self.model.fs.condenser.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.condenser.inlet.pressure[0].setub(high_side_pressure_max)
            self.model.fs.condenser.outlet.pressure[0].setub(high_side_pressure_max)

        use_condenser_temperature_bounds = not ((ambient_temperature is not None) and (condenser_approach is not None))
        if use_condenser_temperature_bounds and check_input(condenser_temperature):
            if condenser_temperature[0] is not None:
                self.model.fs.condenser.Tmin.set_value(condenser_temperature[0] + C_to_K)
                self.model.fs.condenser.T_lower_bound.activate()
            if condenser_temperature[1] is not None:
                self.model.fs.condenser.Tmax.set_value(condenser_temperature[1] + C_to_K)
                self.model.fs.condenser.T_upper_bound.activate()
        else:
            self.model.fs.condenser.T_lower_bound.deactivate()
            self.model.fs.condenser.T_upper_bound.deactivate()

        if subcooling > 0.1:
            self.model.fs.condenser.subcooling.set_value(subcooling)
            self.model.fs.condenser.subcooling_constraint.activate()
            self.model.fs.condenser.subcooling_max.set_value(subcooling_max)
            self.model.fs.condenser.subcooling_lower_constraint.activate()
        else:
            self.model.fs.condenser.subcooling_constraint.deactivate()
            self.model.fs.condenser.subcooling_lower_constraint.deactivate()

        if (ambient_temperature is not None) and (condenser_approach is not None):
            self.model.fs.condenser.ambient_T.set_value(ambient_temperature + C_to_K)
            self.model.fs.condenser.approach_T.set_value(condenser_approach)
            self.model.fs.condenser.approach_constraint.activate()
        else:
            self.model.fs.condenser.approach_constraint.deactivate()

        if debug_disable_arc_pressure_eq:
            if evap_sat_temperature is not None:
                sat_evap = self._sat_at_t_robust(evap_sat_temperature + C_to_K)
                if sat_evap["status"] == "CONVERGED":
                    self.model.fs.P_low_target.set_value(sat_evap["P_Pa"])
                    self.model.fs.P_low.set_value(sat_evap["P_Pa"])
            if (ambient_temperature is not None) and (condenser_approach is not None):
                t_cond_sat = ambient_temperature + condenser_approach + C_to_K
                sat_cond = self._sat_at_t_robust(t_cond_sat)
                if sat_cond["status"] == "CONVERGED":
                    self.model.fs.P_high_target.set_value(sat_cond["P_Pa"])
                    self.model.fs.P_high.set_value(sat_cond["P_Pa"])

            for arc in [self.model.fs.evaporator_to_compressor_expanded,
                        self.model.fs.compressor_to_condenser_expanded,
                        self.model.fs.condenser_to_expansion_valve_expanded,
                        self.model.fs.expansion_valve_to_evaporator_expanded]:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.deactivate()

            self.model.fs.P_low_target_constraint.activate()
            self.model.fs.P_high_target_constraint.activate()
            self.model.fs.P_low_evap_in.activate()
            self.model.fs.P_low_evap_out.deactivate()
            self.model.fs.P_low_comp_in.activate()
            self.model.fs.P_low_valve_out.activate()
            self.model.fs.P_high_comp_out.activate()
            self.model.fs.P_high_cond_in.activate()
            self.model.fs.P_high_cond_out.deactivate()
            self.model.fs.P_high_valve_in.activate()
        else:
            for arc in [self.model.fs.evaporator_to_compressor_expanded,
                        self.model.fs.compressor_to_condenser_expanded,
                        self.model.fs.condenser_to_expansion_valve_expanded,
                        self.model.fs.expansion_valve_to_evaporator_expanded]:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.activate()
            for c in [self.model.fs.P_low_target_constraint, self.model.fs.P_high_target_constraint,
                      self.model.fs.P_low_evap_in, self.model.fs.P_low_evap_out,
                      self.model.fs.P_low_comp_in, self.model.fs.P_low_valve_out,
                      self.model.fs.P_high_comp_out, self.model.fs.P_high_cond_in,
                      self.model.fs.P_high_cond_out, self.model.fs.P_high_valve_in]:
                c.deactivate()

        if low_side_pressure_min:
            self.model.fs.expansion_valve.outlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.expansion_valve.outlet.pressure[0].setub(low_side_pressure_max)
        if high_side_pressure_min:
            self.model.fs.expansion_valve.inlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.expansion_valve.inlet.pressure[0].setub(high_side_pressure_max)

        if check_input(expansion_valve_temperature):
            if expansion_valve_temperature[0] is not None:
                self.model.fs.expansion_valve.Tmin.set_value(expansion_valve_temperature[0] + C_to_K)
                self.model.fs.expansion_valve.T_lower_bound.activate()
            if expansion_valve_temperature[1] is not None:
                self.model.fs.expansion_valve.Tmax.set_value(expansion_valve_temperature[1] + C_to_K)
                self.model.fs.expansion_valve.T_upper_bound.activate()

        calculate_scaling_factors(self.model)

    def optimize_COP(self, verbose, initialize=True, optimize=True):
        """Faithful port of `vapor_compression.py`'s `optimize_COP`."""
        solver = get_solver()
        solver.options = {'max_iter': 1000, 'tol': 1e-6}

        if initialize:
            self.logger.info("Initializing the flowsheet by solving with no objective...")
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()
            results = solver.solve(self.model, tee=verbose)
            if results.solver.termination_condition == "optimal":
                self.logger.info("Initialization successful")
            else:
                self.logger.error("Initialization failed")
            if verbose:
                self.model.fs.report()
            self.model.fs.cop.set_value(
                self.model.fs.evaporator.heat_duty[0].value / self.model.fs.compressor.work_mechanical[0].value
            )

        self.logger.info("Setting up the optimization problem...")

        if optimize:
            self.model.fs.compute_cop.activate()
            self.model.fs.obj.activate()
        else:
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()

        results = solver.solve(self.model, tee=verbose)

        if not optimize:
            if (self.model.fs.compressor.work_mechanical[0].value is not None
                    and self.model.fs.compressor.work_mechanical[0].value != 0):
                self.model.fs.cop.set_value(
                    self.model.fs.evaporator.heat_duty[0].value / self.model.fs.compressor.work_mechanical[0].value
                )

        if results.solver.termination_condition != "optimal":
            results = solver.solve(self.model, tee=verbose)
        if results.solver.termination_condition != "optimal":
            results = solver.solve(self.model, tee=verbose)

        if results.solver.termination_condition == "optimal":
            self.logger.info("Optimization successful")
            self.logger.info("COP: {:.2f}".format(value(self.model.fs.cop)))
            optimization_converged = True
        else:
            self.logger.error("Optimization failed")
            diag = DiagnosticsToolbox(self.model, constraint_residual_tolerance=1e-6)
            diag.display_constraints_with_large_residuals()
            optimization_converged = False

        if verbose:
            self.model.fs.report()

        self.optimization_converged = optimization_converged

        cop_full = value(self.model.fs.cop)
        plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        self._last_cop_full = cop_full
        self._last_plf = plf
        self._last_cop_part = plf * cop_full
        return value(self.model.fs.cop), optimization_converged

    def get_full_load_cop(self):
        """Return solved full-load COP. Identical to `SimpleVaporCompressionCyclePLR.get_full_load_cop`."""
        return value(self.model.fs.cop)

    def get_part_load_cop(self):
        """Return PLR-corrected COP. Identical to `SimpleVaporCompressionCyclePLR.get_part_load_cop`."""
        plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        return plf * value(self.model.fs.cop)

    def report_solution(self):
        """Port of `vapor_compression.py`'s `report_solution`, using
        `entr_mass`/the real envelope in place of `HelmholtzParameterBlock`'s
        diagram helpers (see module docstring, adaptation 2)."""
        print("Optimized COP:", round(value(self.model.fs.cop), 3))

        n = len(self.h_init)
        h_sol = np.zeros(n); p_sol = np.zeros(n); T_sol = np.zeros(n); S_sol = np.zeros(n)

        for i, unit in enumerate(self.unit_operations):
            h_sol[i] = unit.control_volume.properties_out[0].enth_mass()
            p_sol[i] = unit.outlet.pressure[0].value
            T_sol[i] = unit.control_volume.properties_out[0].temperature()
            S_sol[i] = unit.control_volume.properties_out[0].entr_mass()

        def add_warning():
            if self.optimization_converged is False:
                xlim = plt.gca().get_xlim(); ylim = plt.gca().get_ylim()
                x = xlim[1] - (xlim[1] - xlim[0]) * 0.1
                y = ylim[0] + (ylim[1] - ylim[0]) * 0.1
                plt.text(x, y, "Warning: did not converge", color="red", fontsize=12,
                         bbox=dict(facecolor='white', alpha=0.8), va='bottom', ha='right')

        env = self._compute_envelope()

        plt.figure()
        plt.plot(env["h_l_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "b-")
        plt.plot(env["h_v_Jkg"] / 1000.0, env["P_Pa"] / 1000.0, "r-")
        plt.plot(h_sol / 1000.0, p_sol / 1000.0, 'ko')
        add_warning()
        plt.xlabel("h [kJ/kg]"); plt.ylabel("P [kPa]"); plt.yscale("log")
        plt.title("Solution on R-515B P-h envelope"); plt.show()

        plt.figure()
        plt.plot(env["T_K"], env["P_Pa"] / 1000.0, "k-")
        plt.plot(T_sol, p_sol / 1000.0, 'ko')
        add_warning()
        plt.xlabel("T [K]"); plt.ylabel("P [kPa]"); plt.yscale("log")
        plt.title("Solution on R-515B T-P curve"); plt.show()

        plt.figure()
        plt.plot(env["s_l_JkgK"] / 1000.0, env["T_K"], "b-")
        plt.plot(env["s_v_JkgK"] / 1000.0, env["T_K"], "r-")
        plt.plot(S_sol / 1000.0, T_sol, 'ko')
        add_warning()
        plt.xlabel("s [kJ/kg/K]"); plt.ylabel("T [K]")
        plt.title("Solution on R-515B T-s envelope"); plt.show()

        for unit in self.unit_operations:
            unit.report()
