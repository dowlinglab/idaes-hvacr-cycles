# COPY for reference-state isolation testing -- see BREADCRUMB_07-20.md,
# 2026-08-10 update. Identical to vapor_compression_cubic.py except it
# imports make_config/METHODS from phase_1_cubic_eos_validation_refstate
# (which carries the F/G reference-state calibration) instead of the
# original phase_1_cubic_eos_validation. Original file is untouched.
#
# 08/26/2026: retargeted to be the _GCN working copy -- the only change
# from vapor_compression_cubic_refstate.py (production) is the import
# below, now pointing at phase_1_cubic_eos_validation_refstate_GCN
# instead of phase_1_cubic_eos_validation_refstate, so method="first_
# principle"/"gcn" resolve to real METHODS entries (the Colon group's
# two new Shomate fits, replacing SPGP) instead of raising
# `assert method in METHODS`. Used by phase6_final_GCN.py.

# Import required IDAES-PSE modules
from idaes.core import FlowsheetBlock
# Phase 3b (3g): Helmholtz property-package imports removed -- the cycle now uses
# the generic cubic PR package via make_config (imported below); general_helmholtz
# is no longer a dependency of this file.
from idaes.models.unit_models import (Heater, Turbine, Compressor, 
                                      Mixer, Separator, PressureChanger,
                                      Valve)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.exceptions import InitializationError
from pyomo.environ import ConcreteModel, value, Objective, SolverFactory, maximize, minimize, TransformationFactory, Param, Var, Constraint
from idaes.core.util.initialization import propagate_state
from pyomo.network import Arc
from pyomo.environ import units as pyunits
import numpy as np
import matplotlib.pyplot as plt

import logging
import CoolProp.CoolProp as CP

from idaes.core.util import DiagnosticsToolbox
from enum import Enum

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
import os, sys; sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from phase_1_cubic_eos_validation_refstate_GCN import make_config, METHODS

class Mode(Enum):
    ORIGINAL_TPX = "original_TPx"
    IMPROVED_TPX = "improved_TPx"
    PH = "PH"

# Conversion factor
C_to_K = 273.15

class SimpleVaporCompressionCycle:



    def __init__(self,fluid_name, compressor_efficiency=0.75, mode=Mode.IMPROVED_TPX, method="NIST"):
        ''' Simple Vapor Compression Cycle

        Parameters:
            fluid_name : str
                Name of the fluid
            compressor_efficiency : float
                Isentropic efficiency of the compressor
            method : str
                Which R-32 critical-property/Shomate parameter set to use
                (Colon-group collaboration): "NIST", "GCGP", or "SPGP".
                See phase_1_cubic_eos_validation.py's METHODS dict. Defaults
                to "NIST" (the method Phases 1-2 validated as best matching
                the Linde reference dome). Phase 4 comparison: run this
                class once per method at the same spec to see how much the
                property-set choice itself moves the cycle's COP.

        '''

        self.fluid_name = fluid_name

        # Create the ConcreteModel and Flowsheet
        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        # Save the mode
        self.mode = mode

        # Phase 3b (3f): FTPx has no pressure-enthalpy state, so PH mode is invalid
        # with the generic package. Fail fast rather than relying on dead branches.
        assert mode != Mode.PH, "Mode.PH is not supported with the generic PR (FTPx) package"

        # Phase 3b: cubic PR generic property package (Phases 0-2) instead of
        # Helmholtz. FTPx (molar) state; no mass-basis / PH state option.
        assert method in METHODS, f"method must be one of {list(METHODS.keys())}, got {method!r}"
        self.method = method
        self.model.fs.properties = GenericParameterBlock(**make_config(METHODS[method]))

        # Save the compressor efficiency
        assert 0 < compressor_efficiency < 1, "Compressor efficiency must be between 0 and 1"
        self.compressor_efficiency = compressor_efficiency

        self.optimization_converged = None
        
        self._define_flowsheet()


    def _define_flowsheet(self):

        # Set up logging
        logging.basicConfig(level=logging.WARNING)
        self.logger = logging.getLogger(__name__)

        # Add unit models to the flowsheet
        self.model.fs.evaporator = Heater(property_package=self.model.fs.properties)
        self.model.fs.compressor = Compressor(property_package=self.model.fs.properties)
        self.model.fs.condenser = Heater(property_package=self.model.fs.properties)

        self.model.fs.expansion_valve = PressureChanger(property_package=self.model.fs.properties, 
                                thermodynamic_assumption="adiabatic",
                                compressor=False)

        # Let's use a valve instead of a pressure changer
        # self.model.fs.expansion_valve = Valve(property_package=self.model.fs.properties)
                            
        # Connect components with arcs
        self.model.fs.evaporator_to_compressor = Arc(source=self.model.fs.evaporator.outlet, 
                                 destination=self.model.fs.compressor.inlet)
        self.model.fs.compressor_to_condenser = Arc(source=self.model.fs.compressor.outlet, 
                                destination=self.model.fs.condenser.inlet)
        self.model.fs.condenser_to_expansion_valve = Arc(source=self.model.fs.condenser.outlet, 
                                 destination=self.model.fs.expansion_valve.inlet)
        self.model.fs.expansion_valve_to_evaporator = Arc(source=self.model.fs.expansion_valve.outlet, 
                                  destination=self.model.fs.evaporator.inlet)

        # Expand arcs to build the connectivity
        TransformationFactory('network.expand_arcs').apply_to(self.model)

        # Deactivate the flowrate constraint on one of the arcs
        # Out flowsheet is a closed, circular loop
        self.model.fs.evaporator_to_compressor_expanded.flow_mol_equality.deactivate()

        # Let's see if this helps with convergence
        # self.model.fs.evaporator_to_compressor_expanded.pressure_equality.deactivate()
        # self.model.fs.expansion_valve_to_evaporator_expanded.pressure_equality.deactivate()
        # self.model.fs.expansion_valve_to_evaporator_expanded.temperature_equality.deactivate()

        # Phase 3b: same closed-loop redundancy as flow_mol_equality above, but for
        # COMPOSITION instead of flow. R32 is pure, so mole_frac_comp["R32"] == 1
        # everywhere, always -- it never changes across any unit. Each unit's own
        # OUTLET state block carries a "sum_mole_frac_out" constraint enforcing that
        # same trivial fact locally, on top of the 4 arc-level mole_frac_comp_equality
        # constraints already tying composition together all the way around the
        # closed loop. Going all the way around, that's one redundant equation per
        # unit -- confirmed via DiagnosticsToolbox: overall model DOF was -4 after
        # set_specifications(), and deactivating these 4 constraints (one per unit
        # outlet) brings it to exactly 0. Unlike flow (only needs one arc broken),
        # composition redundancy shows up locally at every unit's own outlet block,
        # so all 4 need deactivating, not just one arc.
        #
        # CAVEAT: this deactivation does NOT reliably stick past initialize().
        # Confirmed via diagnostic: each unit's own initialize() call reactivates
        # its own properties_out sum_mole_frac_out as part of its internal
        # bootstrapping, regardless of what we set here. The deactivation is
        # re-asserted at the end of set_specifications() (right before
        # calculate_scaling_factors), which is the one spot guaranteed to run
        # last, before the real solve -- that's the one that actually matters.
        # Kept here too so the model's DOF is also correct immediately after
        # construction, before any initialize() call.
        for unit_name in ["evaporator", "compressor", "condenser", "expansion_valve"]:
            unit = getattr(self.model.fs, unit_name)
            unit.control_volume.properties_out[0.0].sum_mole_frac_out.deactivate()

        # Set up the objective function
        '''
        self.model.fs.COP = Objective(expr=(self.model.fs.evaporator.heat_duty[0]) /
                                (self.model.fs.compressor.work_mechanical[0]), sense=maximize)
        '''
        self.model.fs.cop = Var(initialize=1, units=pyunits.dimensionless, bounds=(0.1, 100))

        @self.model.fs.Constraint(doc="COP constraint")
        def compute_cop(b):
            return b.cop * b.compressor.work_mechanical[0] == b.evaporator.heat_duty[0] 
        
        @self.model.fs.Objective(doc="Maximize COP", sense=maximize)
        def obj(b):
            return b.cop
                                
        self.model.fs.compute_cop.deactivate()
        self.model.fs.obj.deactivate()
        
        self.model.fs.evaporator.superheating = Param(initialize=0, units=pyunits.K, mutable=True)
        self.model.fs.evaporator.T_sat_set = Param(initialize=C_to_K, units=pyunits.K, mutable=True)

        @self.model.fs.evaporator.Constraint(doc="Superheat evaporator outlet")
        def superheating_constraint(b):
            return b.control_volume.properties_out[0].temperature >= b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"] + b.superheating
        
        self.model.fs.evaporator.superheating_constraint.deactivate()

        @self.model.fs.evaporator.Constraint(doc="Evaporator saturation temperature setpoint")
        def evap_sat_constraint(b):
            return b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"] == b.T_sat_set

        self.model.fs.evaporator.evap_sat_constraint.deactivate()

        self.model.fs.condenser.subcooling = Param(initialize=0, units=pyunits.K, mutable=True)
        self.model.fs.condenser.ambient_T = Param(initialize=C_to_K, units=pyunits.K, mutable=True)
        self.model.fs.condenser.approach_T = Param(initialize=0, units=pyunits.K, mutable=True)

        @self.model.fs.condenser.Constraint(doc="Subcool condenser outlet") 
        def subcooling_constraint(b):
            return b.control_volume.properties_out[0].temperature <= b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"] - b.subcooling
        
        self.model.fs.condenser.subcooling_constraint.deactivate()

        @self.model.fs.condenser.Constraint(doc="Condenser saturation temperature setpoint (sat = ambient + approach)")
        def approach_constraint(b):
            return b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"] == b.ambient_T + b.approach_T

        self.model.fs.condenser.approach_constraint.deactivate()

        @self.model.fs.compressor.Constraint(doc="Must be a vapor")
        def vapor_constraint(b):
            return b.control_volume.properties_out[0].temperature >= b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"]
            # return b.outlet.pressure[0]/1e3 >= b.control_volume.properties_out[0].pressure_sat/1e3

        self.model.fs.compressor.vapor_constraint.deactivate()

        # FIX (2026-07-27, pass 5): `vapor_constraint` above (T_out >= Tsat)
        # is too weak on its own -- confirmed via compressor_fix_
        # regression_check.py for BOTH NIST and GCGP that the real outlet
        # can land close to (NIST) or right against (GCGP) a purely
        # numeric Tsat-based floor while reporting essentially the SAME
        # enthalpy as the correctly-superheated isentropic state -- a
        # physically impossible T/h combination for a clean single-phase
        # vapor, and the exact "wrong root near the two-phase boundary"
        # pattern behind nearly every bug this session. A fixed "Tsat+3K"
        # margin (pass 3) isn't reliably far enough from the dome for every
        # method/ambient combination -- GCGP's real outlet landed only
        # 1.46 K above that floor at T_amb=20, while NIST's didn't. Add a
        # STRONGER, physically-exact constraint instead of tuning a numeric
        # margin: the real (inefficient) compression's outlet temperature
        # can never be below the ideal (isentropic) outlet temperature --
        # T_out >= T_isen is a genuine thermodynamic fact, and since both
        # are live model variables (not captured numbers), this tracks
        # correctly regardless of how far the isentropic temperature ends
        # up from Tsat for any given method/ambient.
        @self.model.fs.compressor.Constraint(doc="Real outlet T >= isentropic T")
        def superheat_vs_isentropic_constraint(b):
            return (b.control_volume.properties_out[0].temperature
                    >= b.properties_isentropic[0].temperature)

        self.model.fs.compressor.superheat_vs_isentropic_constraint.deactivate()

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

        # Deactivate by default; activated in set_specifications when requested
        self.model.fs.P_low_target_constraint.deactivate()
        self.model.fs.P_high_target_constraint.deactivate()
        self.model.fs.P_low_evap_in.deactivate()
        self.model.fs.P_low_evap_out.deactivate()
        self.model.fs.P_low_comp_in.deactivate()
        self.model.fs.P_low_valve_out.deactivate()
        self.model.fs.P_high_comp_out.deactivate()
        self.model.fs.P_high_cond_in.deactivate()
        self.model.fs.P_high_cond_out.deactivate()
        self.model.fs.P_high_valve_in.deactivate()

        '''
        # This constraint is redundant with another constraint
        @self.model.fs.expansion_valve.Constraint(doc="Must be two-phase")
        def two_phase_constraint(b):
            # return b.control_volume.properties_out[0].temperature[0] == b.control_volume.properties_out[0].temperature_bubble["Vap", "Liq"]
            return b.outlet.pressure[0]/1e3 == b.control_volume.properties_out[0].pressure_sat/1e3

        self.model.fs.expansion_valve.two_phase_constraint.deactivate()
        '''

        ## Temperature Bounds (needed for PH state variables)

        units = [self.model.fs.evaporator, 
                 self.model.fs.compressor,
                 self.model.fs.condenser,
                 self.model.fs.expansion_valve
        ]

        for u in units:
            # These default values will get reset in the setup method
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

    def draw_thermodynamic_diagrams(self):
        # Phase 3b: hp/pt/ts_diagram are Helmholtz-only, not on the generic PR
        # package. Disabled.
        pass

    def specify_initial_conditions(self,
                                   low_side_temperature = -20, # degC
                                   high_side_temperature = 30 # degC
    ):
    
        ''' Specify initial conditions for the

        Arguments:
            low_side_temperature : float
                Low side temperature in degC
            high_side_temperature : float
                High side temperature in degC

        '''

        # Convert temperatures to Kelvin
        low_side_temperature += C_to_K
        high_side_temperature += C_to_K

        superheat = 3
        subcool = 3

        # Compute saturation pressures
        # Q = 0 for saturated liquid
        # Q = 1 for saturated vapor
        # Does not matter here for a pure component
        # Output is in P
        low_side_pressure = CP.PropsSI('P', 'T', low_side_temperature, 'Q', 0, self.fluid_name)
        high_side_pressure = CP.PropsSI('P', 'T', high_side_temperature, 'Q', 0, self.fluid_name)

        # Expansion valve outlet (assume slightly vaporized)
        low_side_liquid_H = CP.PropsSI('H', 'T', low_side_temperature, 'Q', 0.2, self.fluid_name)

        # Evaporator outlet
        low_side_vapor_H = CP.PropsSI('H', 'T', low_side_temperature + superheat, 'Q', 1, self.fluid_name)
        
        # Compressor outlet
        high_side_liquid_H = CP.PropsSI('H', 'T', high_side_temperature + superheat, 'Q', 0, self.fluid_name)

        # Condenser outlet
        high_side_vapor_H = CP.PropsSI('H', 'T', high_side_temperature - subcool, 'Q', 1, self.fluid_name)

        '''
        self.T_init = np.zeros(len(self.h_init))
        for i in range(len(self.h_init)):
            # Inputs are in J/kg and Pa, hence 1000 is needed for unit conversion
            self.T_init[i] = CP.PropsSI('T', 'H', self.h_init[i]*1000, 'P', self.p_init[i]*1000, self.fluid_name)
        '''

        '''
        Unit Operation (outlets):
            0: Evaporator
            1: Compressor
            2: Condenser
            3: Expansion Valve
        '''

        self.h_init = np.array([low_side_vapor_H, # Evaporator
                       high_side_vapor_H, # Compressor
                       high_side_liquid_H, # Condenser 
                       low_side_liquid_H# Expansion Valve
                       ]) 
        
        self.p_init = np.array([low_side_pressure,
                       high_side_pressure,
                       high_side_pressure,
                       low_side_pressure
                       ])
        
        self.T_init = np.array([low_side_temperature + superheat, # Evaporator
                       high_side_temperature + superheat, # Compressor
                       high_side_temperature - subcool, # Condenser
                       low_side_temperature # Expansion valve
                       ])

        # Phase 3b: Helmholtz-only p-h / p-T / T-s diagrams and their overlays
        # removed (not available on the generic PR package). The h_init/p_init/
        # T_init guess arrays above are still used for initialization.

    def _initialize_compressor_with_retry(self, compressor_state_args, verbose=False):
        """Initialize the compressor by bypassing IDAES's built-in
        ``compressor.initialize()`` entirely, solving the isentropic
        reference state and the real outlet directly and in isolation.

        FIX (2026-07-24): the previous version of this method called IDAES's
        built-in ``init_isentropic()`` (via ``compressor.initialize()``),
        with a single retry using an improved temperature guess if the first
        attempt raised ``InitializationError``. That was NOT robust across
        the Phase 4 ambient sweep: NIST failed at T_amb=15/25 and GCGP failed
        at T_amb=10/25 (both methods' T_amb=20 case, and NIST's T_amb=10/
        GCGP's T_amb=15, happened to work -- a fragile, guess-dependent
        pattern, not a real fix).

        Root cause (confirmed by reading IDAES's own
        ``idaes/models/unit_models/pressure_changer.py`` source directly):
        ``init_isentropic()`` does a brittle 4-step sequence -- (1) init
        inlet/outlet with a shared guess (same T_init[1] guess leaking into
        both, the same class of problem diagnosed for the expansion valve
        earlier this session), (2) init the isentropic pseudo-state off
        THAT SAME guess, (3) TEMPORARILY FIX the isentropic temperature at
        the (possibly bad) guess, deactivate the entropy-matching
        constraint, and solve everything else, (4) unfix/reactivate and
        solve the whole unit for real. If step 3 lands in a bad basin (e.g.
        outlet phase_frac split non-physically, confirmed once at 1.13, a
        real "leftover" sum), step 4 inherits a bad starting point and can
        fail outright -- and simply changing the temperature GUESS (the
        previous retry's only lever) doesn't reliably fix that, since the
        other guessed quantities (pressure ratio behavior, phase split) are
        just as capable of steering step 3 into a bad basin.

        This method instead applies the SAME "bound, not fix; solve blocks
        directly and in isolation" recipe already proven for the expansion
        valve outlet:
          1. Solve the inlet block directly (cheap safety net; it should
             already be consistent from upstream propagation).
          2. Solve the isentropic reference state ALONE against an explicit
             entropy-target constraint (entr_mol == the inlet's CAPTURED
             entropy value -- a captured number, not a live cross-block
             reference, exactly mirroring the valve's isenthalpic-target
             technique). This block is known to be robust regardless of
             guess quality (confirmed via diagnostic earlier this session:
             entropy-matching residual = 0.000000 exactly) BECAUSE it now
             gets a real, unshared, well-bounded warm start instead of
             inheriting the outlet's guess.
          3. Compute the REAL outlet's target enthalpy from the compressor's
             own efficiency relation (matches IDAES's ``actual_work``
             constraint: work_isentropic == work_mechanical *
             efficiency_isentropic, i.e. h_out = h_in + (h_isen - h_in) /
             efficiency), then solve the real outlet ALONE against that
             enthalpy target, with temperature bound-not-fixed within a
             generous but still protective window.
        Composition on both blocks is fixed only TEMPORARILY (to make each
        isolated solve well-posed) and explicitly unfixed again before this
        method returns -- leaving a second, permanently-fixed composition
        variable in the closed loop would reproduce the exact fixed-vs-fixed
        contradiction bug already found and fixed for the expansion valve.

        Args:
            compressor_state_args: dict of initial guesses (flow_mol,
                temperature, pressure, mole_frac_comp). Only "temperature"
                is used here, as a warm-start seed for the isentropic block;
                kept as a parameter for call-site compatibility.
            verbose: if True, print the compressor's report after success.
        """
        comp = self.model.fs.compressor
        comp_in = comp.control_volume.properties_in[0]
        comp_out = comp.control_volume.properties_out[0]
        comp_isen = comp.properties_isentropic[0]

        # Step 1: make sure the inlet is genuinely solved (cheap safety net;
        # it should already be consistent from upstream propagation). Only
        # T/P are fixed by the caller -- flow_mol/mole_frac_comp arrive via
        # propagate_state's value-only copy from the evaporator outlet, so
        # comp_in isn't DOF=0 on its own yet. IDAES's own init_isentropic()
        # handles this via hold_state=True/release_state(); since we're
        # bypassing init_isentropic() entirely, replicate that here:
        # temporarily fix whichever of these aren't already fixed, solve,
        # then release only what we ourselves fixed.
        flow_was_fixed = comp_in.flow_mol.fixed
        comp_was_fixed = comp_in.mole_frac_comp["R32"].fixed
        if not flow_was_fixed:
            comp_in.flow_mol.fix()
        if not comp_was_fixed:
            comp_in.mole_frac_comp["R32"].fix()
        get_solver().solve(comp_in)

        flow_target = value(comp_in.flow_mol)
        comp_target = value(comp_in.mole_frac_comp["R32"])
        entr_in = value(comp_in.entr_mol)
        h_in = value(comp_in.enth_mol)
        T_in = value(comp_in.temperature)
        P_out = value(comp.outlet.pressure[0])
        efficiency = value(comp.efficiency_isentropic[0])

        if not flow_was_fixed:
            comp_in.flow_mol.unfix()
        if not comp_was_fixed:
            comp_in.mole_frac_comp["R32"].unfix()

        relaxed_solver = get_solver(solver_options={
            "tol": 1e-4, "constr_viol_tol": 1e-4, "acceptable_tol": 1e-3
        })

        # Step 2: isentropic reference state, solved alone against a
        # captured-value entropy target (mirrors the valve's isenthalpic
        # seed-constraint technique).
        # FIX (2026-07-24, second pass): the FIRST fix (tightening the real
        # outlet's bound in Step 3) wasn't enough on its own -- re-testing
        # showed the SAME wrong-root-near-the-dome failure had simply moved
        # UP a level, into THIS block. T_isen came out at 302.18 K (~Tsat at
        # P_out) instead of the correct ~310-311 K, with an entropy-match
        # residual just as good (-1.5e-07) as the correct root -- confirming
        # the entropy-matching constraint genuinely has (at least) two
        # nearby solutions here, and a loose [T_in, T_in+200] bound doesn't
        # reliably exclude the spurious one. Use an independent, physically-
        # grounded floor instead of just a generous margin: the isentropic
        # compression of a vapor is ALWAYS superheated relative to the
        # saturation temperature at the discharge pressure, so bound T_isen
        # (and, transitively, T_out) to sit clearly above Tsat at P_out --
        # computed via CoolProp (already used elsewhere in this file for
        # independent real-fluid checks), not the model's own EoS, so this
        # bound can't be fooled by the same degenerate root it's meant to
        # exclude.
        T_sat_high = CP.PropsSI('T', 'P', P_out, 'Q', 1, self.fluid_name)
        comp_isen.flow_mol.fix(flow_target)
        comp_isen.mole_frac_comp["R32"].fix(comp_target)
        comp_isen.pressure.fix(P_out)
        comp_isen.temperature.unfix()
        comp_isen.temperature.setlb(max(T_in, T_sat_high + 3.0))
        comp_isen.temperature.setub(T_in + 200.0)
        T_guess = compressor_state_args.get("temperature", T_sat_high + 20.0)
        comp_isen.temperature.set_value(max(T_sat_high + 3.0, min(T_in + 200.0, T_guess)))

        # TRIED (2026-07-27, pass 8) AND REVERTED: fixed
        # `comp_isen.phase_frac["Vap"]=1.0` HERE (only for this isolated,
        # temporary solve, unfixed again right after) on the theory that
        # even this well-bounded warm-start solve was landing on a wrong
        # (T, quality) pair since phase_frac was left free. Measured
        # effect on GCGP was negligible (T_isen barely moved), which was
        # read as "harmless" and left in place while iterating on pass 7/
        # 9 -- WRONG. After reverting passes 7 and 9, NIST T_amb=20 was
        # STILL badly broken (COP=484980, nonsensical) with only pass 8
        # still active -- proving pass 8 alone, even in this fully
        # isolated/temporary context, corrupts something that propagates
        # forward through the rest of initialize() for NIST specifically.
        # REVERTED. Current shipped state is back to pass-1/2/2b/3 ONLY --
        # phase_frac is not touched anywhere on the isentropic block,
        # isolated or coupled. This is the fully-original state confirmed
        # safe for NIST across the whole grid multiple times earlier this
        # session, before any of passes 4-9 were attempted.
        comp_isen.isentropic_seed_con = Constraint(expr=comp_isen.entr_mol == entr_in)
        res_isen = relaxed_solver.solve(comp_isen)
        self.logger.info(
            f"Compressor isentropic-state direct solve: "
            f"{res_isen.solver.termination_condition}"
        )
        comp_isen.del_component(comp_isen.isentropic_seed_con)
        T_isen = value(comp_isen.temperature)
        h_isen = value(comp_isen.enth_mol)
        # Release everything we temporarily fixed/bounded on this block --
        # properties_isentropic isn't touched by set_specifications()'s
        # general unit unfix loop (that only covers inlet/outlet), so if we
        # don't release these ourselves here, they'd stay fixed/bounded
        # PERMANENTLY through the real coupled solve. The unit's own
        # isentropic_pressure/isentropic (entropy) constraints already tie
        # this block to the real outlet/inlet during that solve; nothing
        # here needs to remain fixed once we have our captured values.
        comp_isen.flow_mol.unfix()
        comp_isen.mole_frac_comp["R32"].unfix()
        comp_isen.pressure.unfix()
        comp_isen.temperature.setlb(None)
        comp_isen.temperature.setub(None)

        # Step 3: real outlet, solved alone against the efficiency-derived
        # enthalpy target, temperature bound-not-fixed around a physically
        # motivated guess.
        h_out_target = h_in + (h_isen - h_in) / efficiency
        # Estimate the real temperature rise by applying the SAME
        # efficiency-derived scaling to the temperature rise as to the
        # enthalpy rise (exact for a constant-Cp ideal gas; a good estimate
        # otherwise -- only used as a warm start/bound center, not a hard
        # constraint).
        T_guess_real = T_in + (T_isen - T_in) / efficiency

        comp_out.flow_mol.fix(flow_target)
        comp_out.mole_frac_comp["R32"].fix(comp_target)
        comp_out.pressure.fix(P_out)
        comp_out.temperature.unfix()
        # FIX (2026-07-24): a wide bound here (previously [T_in,
        # T_isen+100]) let Newton settle on a WRONG root -- the real outlet
        # landing almost exactly at the condensing saturation temperature
        # (a two-phase dome point) despite having essentially the SAME
        # enthalpy as the isentropic state, which is only possible if it's
        # sitting in the degenerate two-phase region rather than
        # superheated vapor (confirmed via
        # compressor_fix_regression_check.py: T_isen=311.29 K, matching the
        # earlier validated ~310 K outlet, but the real outlet came out at
        # 302.18 K -- essentially Tsat -- with h differing by only 0.3
        # J/mol from the isentropic state). Same "wrong root near the
        # two-phase boundary" failure mode already fixed for the valve/
        # evaporator/condenser outlets earlier this session; fixed the same
        # way, with a bound tight enough to physically EXCLUDE the dome
        # branch rather than just warm-starting near the right answer:
        # for a compressor, actual work is always >= ideal (isentropic)
        # work, so the real outlet temperature can never be BELOW T_isen --
        # use that as a hard, physically-justified lower bound (not just a
        # margin), which safely excludes the dome as long as T_isen itself
        # is superheated (guaranteed, since isentropic compression of vapor
        # is always superheated).
        margin_hi = max(30.0, 0.25 * abs(T_guess_real - T_in))
        comp_out.temperature.setlb(T_isen)
        comp_out.temperature.setub(T_guess_real + margin_hi)
        comp_out.temperature.set_value(T_guess_real)
        comp_out.actual_outlet_seed_con = Constraint(expr=comp_out.enth_mol == h_out_target)
        res_out = relaxed_solver.solve(comp_out)
        self.logger.info(
            f"Compressor real-outlet direct solve: {res_out.solver.termination_condition}"
        )
        comp_out.del_component(comp_out.actual_outlet_seed_con)
        comp_out.mole_frac_comp["R32"].unfix()

        # Protect the solved temperature until set_specifications()'s
        # general unfix loop releases it (same pattern as every other unit).
        comp_out.temperature.setlb(None)
        comp_out.temperature.setub(None)
        comp_out.temperature.fix(value(comp_out.temperature))

        if verbose:
            comp.report()
        if verbose:
            self.model.fs.compressor.report()

    def initialize(self, verbose=False):
        ''' Initialize the flowsheet '''


        # p_init has units of Pa, hence the scale factor is 1 to get Pa
        p_scale = 1

        ## Evaporator
        self.model.fs.evaporator.inlet.flow_mol[0].fix(1)   # Example value

        self.model.fs.evaporator.inlet.pressure[0].fix(self.p_init[-1]*p_scale)

        if self.mode == Mode.PH:
            # Initialize H
            self.model.fs.evaporator.inlet.enth_mol[0].fix(self.h_init[-1])
            self.model.fs.evaporator.outlet.enth_mol[0].fix(self.h_init[0])
        else:
            # The evaporator INLET is the post-expansion-valve state: it sits ON
            # the saturation dome (T = T_init[-1], P = Psat, already fixed above)
            # and is genuinely TWO-PHASE. For a pure fluid with 2 coexisting
            # phases, Gibbs phase rule gives F = C - P + 2 = 1 - 2 + 2 = 1: only
            # ONE independent variable is free on the dome (T and P are not
            # independent there). Fixing BOTH T and P together (as the old code
            # did) over-specifies this state and the flash goes "locally
            # infeasible" -- confirmed root cause, see
            # phase_3b_flash_seed_test.py Attempts 1-2.
            #
            # phase_frac is not exposed on the generic FTPx port, so we operate
            # on the control-volume property block directly (not evaporator.inlet).
            evap_in = self.model.fs.evaporator.control_volume.properties_in[0]

            T_target = self.T_init[-1]   # expected saturation temperature (K), our guess

            # Leave T FREE instead of fixed -- P is already fixed, so on the
            # dome T is determined by the EoS's own equal-fugacity condition,
            # not by us.
            evap_in.temperature.unfix()

            # Anchor T near the true value with a TIGHT temporary bound
            # (+/- 5 K). Without this, freeing T lets the solver wander off to
            # a completely different, unphysical single-phase point far from
            # the dome (confirmed: Attempt 4 without a bound converged
            # "optimal" but at T = 173 C, with both phases collapsing to the
            # same value -- the "trivial solution").
            evap_in.temperature.setlb(T_target - 5.0)
            evap_in.temperature.setub(T_target + 5.0)

            # Warm-start AT the expected value so Newton's first step is
            # already in the correct basin (the real two-phase solution), not
            # the trivial one.
            evap_in.temperature.set_value(T_target)

            # Fix the vapor fraction as a NUMERICAL SEED ONLY -- not a
            # physical constraint. This is what actually prevents the solver
            # from collapsing both phases to the same value (the trivial
            # root, which mathematically satisfies equal-fugacity but is
            # unphysical). It gets UNFIXED right after evaporator.initialize()
            # below; the real quality is whatever the upstream expansion
            # valve's isenthalpic balance determines once the cycle is
            # coupled -- not this 0.2 guess.
            evap_in.phase_frac["Vap"].fix(0.2)

            # Evaporator outlet is single-phase (superheated vapor) -- fixing
            # both T and P here is fine, no degeneracy, unchanged from the
            # original code.
            self.model.fs.evaporator.outlet.temperature[0].fix(self.T_init[0])
            evap_in.mole_frac_comp["R32"].fix(1.0)    

            assert degrees_of_freedom(evap_in) ==0,"evaporator inlet DOF !=0 before solve"
            res = get_solver().solve(evap_in)
            self.logger.info(f"Evaporator inlet 2-phase solve:{res.solver.termination_condition}")

            # Revert the evaporator inlet to the STANDARD (T, P)-given spec now
            # that a real, converged two-phase point has been found. The
            # quality/bound scaffolding above was only needed to survive this
            # one initialization sub-solve -- discard it so the rest of the
            # cycle sees a normally specified state.
            evap_in = self.model.fs.evaporator.control_volume.properties_in[0]

            # Quality is no longer fixed -- it becomes flowsheet-determined
            # once the evaporator is coupled to the upstream expansion valve.
            evap_in.phase_frac["Vap"].unfix()

            # Remove the temporary +/- 5 K box so it can't clip the real
            # operating range once the full cycle is solved.
            evap_in.temperature.setlb(None)
            evap_in.temperature.setub(None)

            # Fix T at whatever value it actually converged to (should be
            # close to T_target, e.g. within ~1 K per Attempt 5) -- restores
            # the same (T, P)-fixed pattern used for every other, single-phase
            # state.
            evap_in.temperature.fix(value(evap_in.temperature))

        self.logger.info("Initializing evaporator...")
        self.model.fs.evaporator.initialize(outlvl=logging.WARNING)

        

        if verbose:
            self.model.fs.evaporator.report()

        propagate_state(self.model.fs.evaporator_to_compressor)

        ## Compressor
        self.model.fs.compressor.inlet.pressure[0].fix(self.p_init[0]*p_scale)

        if self.mode == Mode.PH:
            self.model.fs.compressor.inlet.enth_mol[0].fix(self.h_init[0])
        else:
            self.model.fs.compressor.inlet.temperature[0].fix(self.T_init[0])

        self.model.fs.compressor.outlet.pressure[0].fix(self.p_init[1]*p_scale) # Set to target
        # self.model.fs.compressor.control_volume.properties_out[0].phase_frac["Vap"].fix(1.0)  # Ensure vapor phase

        self.model.fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)

        self.logger.info("Initializing compressor...")
        compressor_state_args = {
            "flow_mol": 1.0,
            "temperature": self.T_init[1],
            "pressure": value(self.model.fs.compressor.inlet.pressure[0]),
            "mole_frac_comp": {"R32": 1.0},
        }
        self._initialize_compressor_with_retry(compressor_state_args, verbose=verbose)

        propagate_state(self.model.fs.compressor_to_condenser)

        ## Condenser
        
        self.model.fs.condenser.inlet.pressure[0].fix(self.p_init[1]*p_scale)

        if self.mode == Mode.PH:
            self.model.fs.condenser.inlet.enth_mol[0].fix(self.h_init[1])
            self.model.fs.condenser.outlet.enth_mol[0].fix(self.h_init[2])
        else:
            self.model.fs.condenser.inlet.temperature[0].fix()
            self.model.fs.condenser.outlet.temperature[0].fix(self.T_init[2]-15.0)

            # self.model.fs.condenser.control_volume.properties_out[0].phase_frac["Vap"].fix(0.0)  # Ensure liquid phase (is this needed?)

        self.logger.info("Initializing condenser...")
        # self.model.fs.condenser.initialize(outlvl=logging.WARNING)
        # Modified to rectify numerical tolerance issue with residuals for IPOPT.
        # Diagnosed via a DEBUG-level trace: the "locally infeasible" failure here
        # had a tiny residual (~2.5e-6) -- consistent with the same pure-fluid
        # log_mole_frac_tbub/tdew floating-point dust seen as warnings throughout
        # this file, not a real structural problem (compare to the compressor bug,
        # where the residual was ~1.0, a genuinely broken state). Confirmed in
        # isolated testing: relaxing tol/constr_viol_tol/acceptable_tol for just
        # this init-only solve lets it converge to the physically correct,
        # near-pure-liquid outlet state (phase_frac[Vap] ~ 0.00005). Does not
        # affect the precision of the real, final coupled solve later.

        self.model.fs.condenser.initialize(
            outlvl=logging.WARNING,
            optarg={"tol": 1e-4, "constr_viol_tol": 1e-4, "acceptable_tol": 1e-3},
        )
        if verbose:
            self.model.fs.condenser.report()

        propagate_state(self.model.fs.condenser_to_expansion_valve)

        ## Expansion Valve
        self.model.fs.expansion_valve.inlet.pressure[0].fix(self.p_init[2]*p_scale)
        self.model.fs.expansion_valve.outlet.pressure[0].fix(self.p_init[3]*p_scale)

        # Phase 3b: the valve's inlet is the SAME physical state as the
        # condenser's outlet (subcooled liquid, T=284.15K here, comfortably
        # below the dome at this pressure) -- but despite inheriting the
        # exact same T, P, composition via propagate_state(), its own
        # unconstrained flash was wandering to a wrong-branch answer
        # (phase_frac[Vap] -> ~1.0, entropy off by ~50 J/mol/K from the
        # condenser outlet feeding it) rather than the correct near-pure-
        # liquid answer. Confirmed via diagnostic (valve_branch_debug.py):
        # the INITIAL GUESS was correct (phase_frac[Vap] ~ 1e-5, matching
        # FTPx's own tbub/tdew-based rule) -- it's the numerical SOLVE that
        # drifts away from it, not a bad starting point. This is the same
        # general failure mode as the evaporator's original "trivial root"
        # bug: a mathematically valid but unphysical solution exists nearby,
        # and a good guess alone doesn't stop Newton's method from wandering
        # to it.
        #
        # Fix: bound (not fix) phase_frac["Vap"] with a tight ceiling during
        # this one initialization call, fencing the solver away from the
        # wrong branch. A bound, unlike a fix, doesn't consume a degree of
        # freedom -- fixing it outright was tried first and immediately hit
        # a DOF=-1 BurntToast error, since T/P/flow/composition are already
        # fully fixed here (unlike the evaporator's genuinely two-phase
        # inlet, this state isn't actually on the dome, so there's no real
        # nonzero quality to fix to in the first place -- a bound is the
        # thermodynamically appropriate tool, not a fix).
        valve_in = self.model.fs.expansion_valve.control_volume.properties_in[0]

        # Confirmed by reading IDAES's own init_adiabatic() (pressure_changer.py):
        # it passes the SAME state_args dict we give expansion_valve.initialize()
        # to BOTH properties_in.initialize() and properties_out.initialize()
        # (deriving state_args_out mostly by adjusting pressure). fix_state_vars()
        # only skips a variable if it's ALREADY fixed -- otherwise it force-fixes
        # it using state_args's value. propagate_state() (confirmed via its
        # source) only ever sets .value, never .fix() -- so flow_mol/
        # mole_frac_comp/temperature on this inlet were never actually fixed,
        # just value-copied from the condenser outlet. That's harmless as long
        # as expansion_valve.initialize() is called with no state_args (the
        # inlet's unfixed vars just keep their already-correct propagated
        # value) -- but now that we pass an explicit state_args below (needed
        # for the OUTLET's warm start), it gets applied to this inlet too,
        # overwriting its correct propagated temperature with the outlet's
        # target. Fix: explicitly fix these at their current (already correct)
        # values, same as pressure already is, so fix_state_vars() skips them.
        valve_in.flow_mol.fix()
        valve_in.mole_frac_comp["R32"].fix()
        valve_in.temperature.fix()

        valve_in.phase_frac["Vap"].setub(0.01)

        # FINAL FIX (2026-07-23, night): expansion_valve.initialize()'s
        # built-in joint (whole-unit) solve proved fundamentally unreliable
        # for this outlet, no matter how T/phase_frac were bounded.
        # Diagnosed directly with tee=True + DiagnosticsToolbox
        # (valve_joint_solve_debug.py): that built-in solve re-solves the
        # ENTIRE unit -- inlet and outlet together. Even though the inlet's 4
        # canonical state vars are fixed, its DERIVED flash variables
        # (phase-specific flows, log mole fractions) are still free, and a
        # badly-guessed outlet warm start dragged them into a bad
        # restoration spiral (50 Ipopt iterations) ending in "locally
        # infeasible" -- confirmed by a huge residual specifically on the
        # INLET's own component_flow_balances constraint (~0.989), nothing
        # to do with the outlet's bounds (two rounds of bound-tweaking were
        # tried and ruled out first).
        #
        # Fix: stop relying on expansion_valve.initialize() for this unit at
        # all. Solve the inlet's own flash directly first (0 DOF, given
        # T/P/flow/comp already fixed above), then solve the OUTLET
        # directly and in ISOLATION -- never touching the inlet's internals
        # -- against the one equation that actually governs an adiabatic
        # throttle: inlet enthalpy == outlet enthalpy, using the valve's OWN
        # actual inlet enthalpy (not a value borrowed from the evaporator,
        # which was the earlier, wrong approach). Confirmed working in
        # isolation via valve_outlet_test5.py ("optimal", isenthalpic gap
        # ~0.007 J/mol).
        valve_out = self.model.fs.expansion_valve.control_volume.properties_out[0]
        evap_in_solved = self.model.fs.evaporator.control_volume.properties_in[0]
        T_target = value(evap_in_solved.temperature)  # bound center/warm-start only

        res_in = get_solver().solve(valve_in)
        self.logger.info(f"Expansion valve inlet direct solve: {res_in.solver.termination_condition}")

        h_target = value(valve_in.enth_mol)
        flow_target = value(valve_in.flow_mol)
        comp_target = value(valve_in.mole_frac_comp["R32"])

        # Same composition-loop redundancy as Task #24 -- must stay
        # deactivated, not reactivated (confirmed: reactivating it creates a
        # genuine contradiction, 1.0 exactly vs the inlet's actual 0.999988,
        # tiny numerical dust visible in every diagnostic run this session).
        valve_out.sum_mole_frac_out.deactivate()

        valve_out.flow_mol.unfix()
        valve_out.mole_frac_comp["R32"].unfix()
        valve_out.temperature.unfix()
        valve_out.phase_frac["Vap"].unfix()

        valve_out.flow_mol.fix(flow_target)
        valve_out.mole_frac_comp["R32"].fix(comp_target)

        # T bound is just a fence to keep Newton in the right basin -- the
        # ACTUAL value is determined by the isenthalpic constraint below,
        # not by this guess (unlike the earlier, wrong approach that fixed T
        # outright at a value borrowed from the evaporator).
        valve_out.temperature.setlb(T_target - 10.0)
        valve_out.temperature.setub(T_target + 10.0)
        valve_out.temperature.set_value(T_target)

        # The actual physical requirement for an adiabatic throttle, added
        # as an explicit, temporary constraint.
        valve_out.isenthalpic_seed_con = Constraint(expr=valve_out.enth_mol == h_target)

        assert degrees_of_freedom(valve_out) == 0, "valve outlet DOF != 0 before isenthalpic solve"
        # Relaxed tolerance -- same fix already confirmed for the condenser:
        # the near-converged point here trips Ipopt's strict default
        # tolerance on tiny residuals (~1e-5), not a real infeasibility.
        outlet_solver = get_solver(solver_options={"tol": 1e-4, "constr_viol_tol": 1e-4,
                                                    "acceptable_tol": 1e-3})
        res_out = outlet_solver.solve(valve_out)
        self.logger.info(f"Expansion valve outlet isenthalpic solve: {res_out.solver.termination_condition}")

        # Clean up: remove the temporary constraint, revert T to the
        # standard fixed pattern (matches every other unit's state before
        # the real coupled solve; set_specifications() unfixes as needed).
        valve_out.del_component(valve_out.isenthalpic_seed_con)
        valve_out.temperature.setlb(None)
        valve_out.temperature.setub(None)
        valve_out.temperature.fix(value(valve_out.temperature))

        # Remove the temporary bound on the inlet now that a real converged
        # point has been found for both ends.
        valve_in.phase_frac["Vap"].setub(None)

        if verbose:
            self.model.fs.expansion_valve.report()

        propagate_state(self.model.fs.expansion_valve_to_evaporator)

        if verbose:
            print("\nFinished initialization. Stream summary:")
            self.model.fs.report()

    def set_specifications(self,
                           low_side_pressure = (200, 500), # kPa
                           high_side_pressure = (1000, 3000), # kPa
                           evaporator_temperature = (-20, 0), # degC
                           compressor_temperature = None, # degC
                           condenser_temperature = (30, 50), # degC
                           expansion_valve_temperature = None, # degC
                           subcooling = 3, # degC
                           superheating = 3, # degC
                           max_pressure_ratio = 4,
                           ambient_temperature = None, # degC
                           condenser_approach = None, # degC
                           evap_sat_temperature = None, # degC
                           debug_disable_arc_pressure_eq = False
                           ):
        # print(superheating)
        assert superheating >= 0, "Superheating must be greater than or equal to 0"
        assert subcooling >= 0, "Subcooling must be greater than or equal to 0"
        assert max_pressure_ratio > 1.2, "Maximum pressure ratio must be greater than 1.2"
        # assert not (Tsat_constraints and bound_vapor_frac), "Cannot have both Tsat constraints and vapor fraction bounds"

        ## Unfix the conditions from initialization
        self.unit_operations = [self.model.fs.evaporator, self.model.fs.compressor, self.model.fs.condenser, self.model.fs.expansion_valve]

        for unit in self.unit_operations:
            # Unfix all variables
            unit.inlet.flow_mol[0].unfix()
            unit.outlet.flow_mol[0].unfix()

            # (2026-07-23, night: tried blanket-unfixing mole_frac_comp for
            # every unit's inlet/outlet here -- REVERTED, made things worse.
            # With sum_mole_frac_out deactivated everywhere (Task #24) AND
            # no fixed composition anywhere, nothing in the closed loop says
            # "composition = 1.0" anymore -- the 4 arc equalities only
            # enforce "all equal to each other," under-constraining the
            # composition subsystem even though the global DOF count still
            # showed 0. Confirmed empirically: residuals got dramatically
            # WORSE across the whole model (one pressure_balance residual
            # hit 15.6). Targeted fix instead, right after the valve outlet
            # section below: unfix ONLY the specific pair that was actually
            # contradictory (valve_out vs evap_in), leaving evap_in as the
            # loop's sole composition anchor.)

            if self.mode == Mode.PH:
                unit.inlet.enth_mol[0].unfix()
                unit.outlet.enth_mol[0].unfix()
            else:
                unit.inlet.temperature[0].unfix()
                unit.control_volume.properties_in[0].phase_frac["Vap"].unfix()
                unit.outlet.temperature[0].unfix()
                unit.control_volume.properties_out[0].phase_frac["Vap"].unfix()
            
            unit.inlet.pressure[0].unfix()
            unit.outlet.pressure[0].unfix()

            unit.T_lower_bound.deactivate()
            unit.T_upper_bound.deactivate()

            # Set bounds for the vapor fraction to ensure it is within [0,1]
            if self.mode == Mode.ORIGINAL_TPX or self.mode == Mode.IMPROVED_TPX:
                unit.control_volume.properties_in[0].phase_frac["Vap"].setlb(0)
                unit.control_volume.properties_in[0].phase_frac["Vap"].setub(1)

        if self.mode == Mode.IMPROVED_TPX:
            self.model.fs.compressor.vapor_constraint.activate()
            # self.model.fs.expansion_valve.two_phase_constraint.activate()
        elif self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.compressor.vapor_constraint.deactivate()
            # self.model.fs.expansion_valve.two_phase_constraint.deactivate()
        # Need to decide what to do for PH here

        def check_input(bounds):
            if bounds is not None and len(bounds) == 2:
                return True
            else:
                return False

        # Determine whether saturation constraints are active
        evap_sat_active = evap_sat_temperature is not None
        cond_sat_active = (ambient_temperature is not None) and (condenser_approach is not None)

        # Convert pressures from kPa to Pa (use wide safety bounds if sat constraints are active)
        if evap_sat_active or cond_sat_active:
            low_side_pressure_min = 50 * 1000
            low_side_pressure_max = 2000 * 1000
            high_side_pressure_min = 100 * 1000
            high_side_pressure_max = 5000 * 1000
        else:
            if check_input(low_side_pressure):
                low_side_pressure_min = low_side_pressure[0]*1000
                low_side_pressure_max = low_side_pressure[1]*1000
            else:
                low_side_pressure_min = None
                low_side_pressure_max = None

            if check_input(high_side_pressure):
                high_side_pressure_min = high_side_pressure[0]*1000
                high_side_pressure_max = high_side_pressure[1]*1000
            else:
                high_side_pressure_min = None
                high_side_pressure_max = None

        # Set mass flowrate to 1 kg/s because we only care about thermodynamic efficiency
        self.model.fs.evaporator.inlet.flow_mol[0].fix(1)

        ## Evaporator

        # Evaporator pressure bounds
        if low_side_pressure_min:
            self.model.fs.evaporator.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.evaporator.inlet.pressure[0].setub(low_side_pressure_max)

        # Evaporator temperature bounds (outlet)
        ## TODO: Set these bounds on the control variable if PH mode
        # FIX (2026-07-23, night): this legacy bound-setting block has a
        # non-None default tuple (-20, 0), so it fired UNCONDITIONALLY even
        # when the newer evap_sat_temperature spec was also active -- the
        # two mechanisms directly contradicted each other (evap_sat_
        # constraint wants T~244K, this legacy bound forbade T<253.15K).
        # Confirmed via coupled_solve_debug.py: evaporator outlet's lb was
        # 253.15 while its actual (correct) value was 247.15K, BEFORE the
        # coupled solve even started -- Ipopt was fighting a contradiction
        # from iteration 0, which is why it ran away to 328.96K (+84.8K off
        # the dome) instead of converging. Guard this block so it only
        # applies when evap_sat_temperature is NOT being used.
        if check_input(evaporator_temperature) and not evap_sat_active:
            if evaporator_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.evaporator.Tmin.set_value(evaporator_temperature[0] + C_to_K)
                    self.model.fs.evaporator.T_lower_bound.activate()

                else:
                    self.model.fs.evaporator.outlet.temperature[0].setlb(evaporator_temperature[0] + C_to_K)
            
            if evaporator_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.evaporator.Tmax.set_value(evaporator_temperature[1] + C_to_K)
                    self.model.fs.evaporator.T_upper_bound.activate()
                else:
                    self.model.fs.evaporator.outlet.temperature[0].setub(evaporator_temperature[1] + C_to_K)

        # Evaporator outlet must be a vapor
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.evaporator.control_volume.properties_out[0].phase_frac["Vap"].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.evaporator.control_volume.properties_out[0].phase_frac["Vap"].fix(1.0)

            # FIX (2026-07-23, night): phase_frac=1.0 EXACTLY is still a
            # point ON the saturation dome (zero liquid, but still a
            # two-phase boundary) -- same Gibbs-phase-rule degeneracy
            # (T,P not independent) behind nearly every bug this session.
            # With T left completely free/unfenced, the coupled solve ran
            # away to 322.9K (+78.8K off the dome) instead of landing near
            # Tsat. Same fix as the valve outlet: a tight bound + warm
            # start around the real saturation target, not a fix (T is
            # still determined by the real energy/mass balance within the
            # fence).
            if evap_sat_active:
                T_evap_sat_K = evap_sat_temperature + C_to_K
                evap_out_state = self.model.fs.evaporator.control_volume.properties_out[0]
                evap_out_state.temperature.setlb(T_evap_sat_K - 10.0)
                evap_out_state.temperature.setub(T_evap_sat_K + 10.0)
                evap_out_state.temperature.set_value(T_evap_sat_K)

        # Phase 3b: eq_complementarity is Helmholtz-only; generic SmoothVLE has no
        # complementarity var -- nothing to deactivate here.

        # Activate superheating constraint
        if superheating > 0.1:
            self.model.fs.evaporator.superheating.set_value(superheating)
            self.model.fs.evaporator.superheating_constraint.activate()
        else:
            self.model.fs.evaporator.superheating_constraint.deactivate()

        # Activate evaporator saturation temperature constraint
        if evap_sat_temperature is not None:
            self.model.fs.evaporator.T_sat_set.set_value(evap_sat_temperature + C_to_K)
            self.model.fs.evaporator.evap_sat_constraint.activate()
        else:
            self.model.fs.evaporator.evap_sat_constraint.deactivate()

        ## Compressor

        # Compressor pressure bounds
        if low_side_pressure_min:
            self.model.fs.compressor.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.compressor.inlet.pressure[0].setub(low_side_pressure_max)

        if high_side_pressure_min:
            self.model.fs.compressor.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.compressor.outlet.pressure[0].setub(high_side_pressure_max)

        # Compressor temperature bounds (outlet)
        ## TODO: Set temperature bounds on control volume in PH mode
        if check_input(compressor_temperature):
            if compressor_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.compressor.Tmin.set_value(compressor_temperature[0] + C_to_K)
                    self.model.fs.compressor.T_lower_bound.activate()
                else:
                    self.model.fs.compressor.outlet.temperature[0].setlb(compressor_temperature[0] + C_to_K)

            if compressor_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.compressor.Tmax.set_value(compressor_temperature[1] + C_to_K)
                    self.model.fs.compressor.T_upper_bound.activate()
                else:
                    self.model.fs.compressor.outlet.temperature[0].setub(compressor_temperature[1] + C_to_K)
        
        # Compressor outlet must be a vapor
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.compressor.control_volume.properties_out[0].phase_frac["Vap"].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.compressor.control_volume.properties_out[0].phase_frac["Vap"].fix(1.0)

            # Phase 3b: eq_complementarity removed (Helmholtz-only)

            # Add inequality constraint to ensure the compressor outlet is only vapor
            self.model.fs.compressor.vapor_constraint.activate()

            # FIX (2026-07-24, third pass): phase_frac["Vap"].fix(1.0) is
            # STILL a point ON the dome (same Gibbs-phase-rule degeneracy as
            # the evaporator/condenser outlets, fixed earlier this session
            # with a +-10K bound around the true target) -- but the
            # compressor never got that same protection. The existing
            # `vapor_constraint` (T_out >= Tsat) is too WEAK: it allows T to
            # sit exactly AT Tsat, which is precisely the degenerate branch
            # confirmed via compressor_fix_regression_check.py (T_isen and
            # T_out both landing within ~0.03 K of Tsat at the discharge
            # pressure, with an entropy-match residual just as clean as the
            # correct ~310-311 K root). Bounding T only during
            # initialize()'s isolated bypass solve (the two earlier fixes
            # this session) wasn't enough -- set_specifications() unfixes
            # everything and gives NO bound back, so the real coupled solve
            # is completely free to (and reliably does) drift back to the
            # degenerate root. Fix the same way as evap/cond: an
            # independent, physically-grounded floor (CoolProp Tsat at the
            # discharge pressure, not the model's own EoS -- can't be fooled
            # by the same degenerate root it's meant to exclude), applied to
            # BOTH the isentropic reference block and the real outlet, and
            # left in place for the whole coupled solve (not reset to None
            # afterward like the transient initialize()-time bounds).
            P_high_now = value(self.model.fs.compressor.outlet.pressure[0])
            T_sat_high_now = CP.PropsSI('T', 'P', P_high_now, 'Q', 1, self.fluid_name)
            comp_out_state = self.model.fs.compressor.control_volume.properties_out[0]
            comp_isen_state = self.model.fs.compressor.properties_isentropic[0]

            # TRIED (2026-07-27, pass 6) AND REVERTED: replaced the flat
            # `Tsat + 3.0 K` floor with a per-case floor anchored to
            # CoolProp's REAL R32 isentropic temperature at this
            # compressor's actual inlet state/discharge pressure
            # (T_in, P_in -> real entropy -> real T at P_high), on the
            # assumption it would sit close to whatever each method's own
            # cubic-PR model predicts, off by only a small EoS-deviation
            # margin (5K). WRONG assumption for GCGP specifically: GCGP
            # deliberately fits DIFFERENT critical constants (Tc/Pc/omega)
            # than real R32 (that's the whole point of comparing multiple
            # parameter sets), so CoolProp's real-fluid estimate is not a
            # valid stand-in for what GCGP's own EoS should predict. The
            # gap turned out to be ~72K (T_isen_real=397.8K vs the
            # model's own ~325.6K), not the assumed ~5K -- so the "safety
            # margin" forced GCGP's solution up to an artificial branch
            # far from its actual self-consistent answer. Single-point
            # check (GCGP, T_amb=20) already showed this clearly, no full
            # grid needed to confirm: COP/Carnot dropped to 0.6847 (worse
            # than pass 5's 0.7165, well below the 0.75-0.78 band),
            # superheat ballooned to +2.089K (target 0), and subcool
            # flipped sign to -1.681K (condenser outlet no longer even
            # subcooled -- a clear spec violation, worse than passes 3-5).
            # REVERTED back to the flat pass-3 floor.
            #
            # Lesson for any future attempt: an external real-fluid
            # anchor (CoolProp) is only valid for methods whose fitted
            # critical constants are close to the real fluid (NIST,
            # basically by construction) -- NOT for methods that
            # deliberately differ (GCGP, SPGP). A next attempt should
            # stay entirely inside each method's OWN model -- e.g. scale
            # the flat margin by that method's own computed pressure
            # ratio (`self.model.fs.compressor.ratioP[0]`), with no
            # external real-fluid reference at all.
            # FIX (2026-07-27, pass 7): the real anomaly (see the MAJOR
            # REFRAME breadcrumb entry) is that the ISENTROPIC block's own
            # entropy-matching solve can land on a spurious root, because
            # its `phase_frac["Vap"]` is left FREE (unlike the real
            # outlet, which gets it fixed to 1.0) -- a two-phase mixture's
            # entropy depends on both T AND quality, so a wide-open T
            # range gives Newton room to match the target entropy at a
            # wrong (T, quality) pair far from the correct single-phase
            # answer. Passes 4-6 all tried to fix this by touching physics
            # (fixing phase_frac, adding a cross-block constraint, or
            # anchoring to external real-fluid data) and each broke other
            # cells. Pass 7 instead does the SAME kind of thing pass 3
            # already proved safe for the real outlet -- narrow the
            # TEMPERATURE bound -- but anchors it to the isentropic
            # block's OWN already-good warm-started value (captured
            # HERE, before this bound is applied) rather than a flat
            # Tsat-derived number. `_initialize_compressor_with_retry()`
            # already does a dedicated, well-posed, CoolProp-Tsat-anchored
            # solve for this exact block during initialize() -- by the
            # time set_specifications() runs, comp_isen_state.temperature
            # should already hold a good value; the wide flat bound
            # (Tsat+3 to Tsat+150, ~147K range) pass 3 gave it was likely
            # what let the SUBSEQUENT coupled solve wander away from that
            # good starting point to a distant spurious root. Tightening
            # the bound around the warm start directly removes that room
            # to wander, without touching phase_frac and without any
            # external real-fluid reference. The real outlet's own bound
            # is left exactly as pass 3 had it (already proven fine on
            # its own -- its degenerate-root problem was always inherited
            # FROM the isentropic block, not independent).
            # TRIED (2026-07-27, pass 7) AND REVERTED: tightened this
            # block's bound to `[max(Tsat+3, T_isen_warm-15), T_isen_warm
            # +30]` (anchored to whatever _initialize_compressor_with_
            # retry() already had as a warm start) instead of the flat
            # Tsat+3/Tsat+150 range. Result for GCGP: T_out-T_isen shrank
            # but stayed large and wrong at 3 of 4 ambients (-14 to -16K),
            # confirming the warm start itself was already on a bad
            # branch (see pass 8). TRIED (2026-07-27, pass 9) AND
            # REVERTED: ALSO fixing `comp_isen_state.phase_frac["Vap"]
            # =1.0` here (combined with pass 7's tight bound, testing
            # whether pass 4's original wide-bound combination was what
            # caused ITS breakage, not the phase_frac fix itself) DID
            # make T_out track T_isen correctly for GCGP (tiny positive
            # gaps at 3 of 4 ambients) -- but shifted the WHOLE compressor
            # calculation onto a different, much colder, still-wrong
            # branch (GCGP T_isen dropped to ~301-333K vs the ~358-376K
            # smooth trend established via phase4_property_state_
            # comparison.py), causing a severe GCGP T15 outlier
            # (COP=2.3331, -34% vs Helmholtz) AND a catastrophic NIST
            # T_amb=20 regression (COP=-12980.08, nonsensical -- T_out-
            # T_isen blew out to +58.5K, the OPPOSITE direction from every
            # prior bug this session). NIST T_amb=20 is the original
            # reference point this entire compressor-fix effort has been
            # validated against -- breaking it this badly is decisive.
            # REVERTED both passes 7 and 9 back to the flat, pass-3-only
            # bound (confirmed safe for NIST across the whole grid
            # multiple times this session). Pass 8 (phase_frac fix inside
            # the ISOLATED warm-start solve in _initialize_compressor_
            # with_retry()) was INITIALLY assessed as harmless/inert based
            # on a single comparison -- that assessment was WRONG. After
            # reverting passes 7 and 9, NIST T_amb=20 was STILL badly
            # broken (COP=484980, nonsensical) with only pass 8 still
            # active, proving it alone corrupts something that propagates
            # forward through initialize() for NIST specifically. Pass 8
            # has ALSO been fully reverted (see the comment at its
            # original location, ~line 551) -- phase_frac is not touched
            # anywhere on the isentropic block, isolated or coupled. This
            # is the fully-original state confirmed safe for NIST across
            # the whole grid multiple times this session.
            comp_isen_state.temperature.setlb(T_sat_high_now + 3.0)
            comp_isen_state.temperature.setub(T_sat_high_now + 150.0)

            comp_out_state.temperature.setlb(T_sat_high_now + 3.0)
            comp_out_state.temperature.setub(T_sat_high_now + 150.0)

            # TRIED (2026-07-27, pass 5) AND REVERTED: activating
            # `superheat_vs_isentropic_constraint` (T_out >= T_isen -- see
            # its definition/rationale in _define_flowsheet()) did fix the
            # GCGP T_amb=20 sign flip in isolation (constraint confirmed
            # genuinely active and binding: T_out-T_isen -> 0.017, both
            # landing at ~374K instead of the old 306.6/325.6K split) --
            # but the full grid came back WORSE overall, not better:
            #   NIST  -9.68/-3.76/-1.83/-3.76%  -> -12.76/-7.99/-3.6/-2.8%
            #   GCGP -11.57/-8.09/+7.65/-1.21%  -> -12.15/-21.7/-5.45/-14.9%
            # NIST degraded at 3 of 4 ambients (T20 roughly doubled, from
            # -1.83% to -3.6%); GCGP's T20 sign flip is gone but overshot
            # past zero to -5.45%, while T15 and T25 both blew out much
            # worse (-21.7%, -14.9%). Root cause of why this generalizes
            # badly: NIST's OWN already-accepted T_amb=20 solution
            # (COP=3.1316, the reference this whole compressor saga was
            # validated against) has T_isen=330.01K > T_out=319.80K --
            # i.e. it ALREADY violates T_out>=T_isen, so turning this
            # constraint on globally forces essentially every cell onto a
            # different branch than the ones already validated, not just
            # the one GCGP cell it was aimed at. Same lesson as pass 4:
            # a fix confirmed at one diagnosed point is not safe to ship
            # without checking the whole grid. REVERTED -- constraint
            # left defined+deactivated in _define_flowsheet() (not
            # activated here) for future revisiting.
            #
            # If revisiting: the underlying T/h inconsistency at GCGP's
            # T_amb=20 (confirmed via compressor_fix_regression_check.py
            # GCGP) is real and unresolved -- a next attempt might scale
            # the Tsat+3K margin (pass 3) by pressure ratio instead of
            # using a flat 3K for every method, or investigate why NIST's
            # own accepted solution tolerates T_isen > T_out while still
            # reporting a low-error COP (is that inconsistency actually
            # harmless there, and if so why does forcing it away hurt
            # NIST's accuracy?).

            # REVERTED (2026-07-24): tried fixing comp_isen_state.phase_frac
            # ["Vap"] to 1.0 here (mirroring the real outlet's treatment),
            # since the T/h inconsistency it was meant to fix DID resolve in
            # the single-point regression check -- but the user reported the
            # full ambient-grid sweep got WORSE with this change in place,
            # not better. Reverted. The Tsat-anchored T bounds above (pass 3)
            # stay -- only this specific phase_frac fix (pass 4) is undone.
            # If revisiting: the T/h inconsistency this was targeting is
            # real (confirmed via compressor_fix_regression_check.py at
            # T_amb=20), but forcing phase_frac=1.0 on the isentropic block
            # apparently over-constrains or shifts the branch selection
            # badly at OTHER ambients -- needs a different fix, not a blanket
            # revert-and-ignore.

        # Compressor only allows input work
        # self.model.fs.compressor.work_mechanical.setlb(0)

        # Set the pressure ratio bounds
        self.model.fs.compressor.ratioP.setub(max_pressure_ratio)
        self.model.fs.compressor.ratioP.setlb(1.1)

        ## Condenser

        # Condenser pressure bounds
        if high_side_pressure_min:
            self.model.fs.condenser.inlet.pressure[0].setlb(high_side_pressure_min)
            self.model.fs.condenser.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.condenser.inlet.pressure[0].setub(high_side_pressure_max)
            self.model.fs.condenser.outlet.pressure[0].setub(high_side_pressure_max)

        # Condenser temperature bounds (outlet)
        # TODO: Write constraint with control volume in PH mode
        # If approach-to-ambient is enforced, skip outlet temperature bounds to avoid over-constraining.
        use_condenser_temperature_bounds = not (
            (ambient_temperature is not None)
            and (condenser_approach is not None)
        )

        if use_condenser_temperature_bounds and check_input(condenser_temperature):
            if condenser_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.condenser.Tmin.set_value(condenser_temperature[0] + C_to_K)
                    self.model.fs.condenser.T_lower_bound.activate()
                else:
                    self.model.fs.condenser.outlet.temperature[0].setlb(condenser_temperature[0]+ C_to_K)
            if condenser_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.condenser.Tmax.set_value(condenser_temperature[1] + C_to_K)
                    self.model.fs.condenser.T_upper_bound.activate()
                else:
                    self.model.fs.condenser.outlet.temperature[0].setub(condenser_temperature[1]+ C_to_K)
        else:
            # Deactivate any previously-activated bounds in PH mode
            if self.mode == Mode.PH:
                self.model.fs.condenser.T_lower_bound.deactivate()
                self.model.fs.condenser.T_upper_bound.deactivate()

        # Condenser outlet must be a liquid
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.condenser.control_volume.properties_out[0].phase_frac["Vap"].setub(0.01)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.condenser.control_volume.properties_out[0].phase_frac["Vap"].fix(0.0)

            # FIX (2026-07-23, night): same reasoning as the evaporator
            # outlet above -- phase_frac=0.0 EXACTLY is still a dome
            # boundary point, and T was left completely unfenced (only the
            # generic [200,450] package default), which is exactly why the
            # condenser outlet was landing 6.9K above its own Tsat instead
            # of at it. Tight bound + warm start around the real
            # ambient+approach saturation target, not a fix.
            if cond_sat_active:
                T_cond_sat_K = ambient_temperature + condenser_approach + C_to_K
                cond_out_state = self.model.fs.condenser.control_volume.properties_out[0]
                cond_out_state.temperature.setlb(T_cond_sat_K - 10.0)
                cond_out_state.temperature.setub(T_cond_sat_K + 10.0)
                cond_out_state.temperature.set_value(T_cond_sat_K)

            # Phase 3b: eq_complementarity removed (Helmholtz-only)

        # Activate subcooling constraint
        if subcooling > 0.1:
            self.model.fs.condenser.subcooling.set_value(subcooling)
            self.model.fs.condenser.subcooling_constraint.activate()
        else:
            self.model.fs.condenser.subcooling_constraint.deactivate()

        # Activate condenser saturation temperature constraint (sat = ambient + approach)
        if (ambient_temperature is not None) and (condenser_approach is not None):
            self.model.fs.condenser.ambient_T.set_value(ambient_temperature + C_to_K)
            self.model.fs.condenser.approach_T.set_value(condenser_approach)
            self.model.fs.condenser.approach_constraint.activate()
        else:
            self.model.fs.condenser.approach_constraint.deactivate()


        ## Expansion Valve

        # FIX (2026-07-23, night): valve_out.mole_frac_comp["R32"] is left
        # FIXED by initialize()'s isenthalpic-solve recipe (0.99998754,
        # snapshotted from the valve's own inlet at that moment) and never
        # unfixed anywhere -- while evap_in.mole_frac_comp["R32"] is
        # SEPARATELY fixed at exactly 1.0 (the original evaporator recipe).
        # Since these two are directly linked by the expansion_valve_to_
        # evaporator arc's mole_frac_comp_equality constraint, having BOTH
        # permanently fixed at different values created a residual that
        # could never shrink no matter how long the solve ran (confirmed
        # bit-for-bit identical across multiple otherwise-very-different
        # solve attempts via coupled_solve_debug.py). Unfix ONLY this one
        # side -- evap_in stays fixed at 1.0 as the closed loop's sole
        # composition anchor (unlike blanket-unfixing composition
        # everywhere, which was tried and reverted above: with
        # sum_mole_frac_out deactivated everywhere too, that left nothing
        # in the whole loop anchoring composition to 1.0 at all).
        self.model.fs.expansion_valve.control_volume.properties_out[0].mole_frac_comp["R32"].unfix()

        # FIX (2026-07-23, later night): the valve outlet is ALSO a
        # dome-boundary two-phase state (isenthalpic throttle end point --
        # same Gibbs-phase-rule T/P degeneracy as the evaporator/condenser
        # outlets above). During initialize()'s isenthalpic-solve bypass it
        # was protected by a tight +-10K bound + warm start, but that bound
        # gets stripped right after (temperature.fix()'d instead), and the
        # general unfix loop above (for unit in self.unit_operations: ...
        # unit.outlet.temperature[0].unfix()) unfixes it again with NO bound
        # reapplied -- leaving it free across only the generic property
        # bounds (~[200,450]). Confirmed via optimize_COP(): after relaxing
        # the solver tol/linear_solver to match coupled_solve_debug.py, the
        # solve reached "Solved To Acceptable Level" in only 10 iterations
        # but landed on a WRONG branch -- expansion_valve_to_evaporator
        # stream at T=291.90K where Tsat~=244K (evaporator inlet should sit
        # at essentially the same Tsat as the evaporator itself, since the
        # valve outlet feeds directly into it). Apply the identical
        # bound-not-fix recipe here, centered on the same T_evap_sat_K
        # already computed above (still in scope -- no block scoping in
        # Python) for the evaporator outlet fix.
        if evap_sat_active:
            valve_out_state = self.model.fs.expansion_valve.control_volume.properties_out[0]
            valve_out_state.temperature.setlb(T_evap_sat_K - 10.0)
            valve_out_state.temperature.setub(T_evap_sat_K + 10.0)
            valve_out_state.temperature.set_value(T_evap_sat_K)

        # Debug: optionally disable arc pressure equalities in the loop
        if debug_disable_arc_pressure_eq:
            # Set pressure targets from saturation setpoints (Pa)
            if evap_sat_temperature is not None:
                P_low = CP.PropsSI('P', 'T', evap_sat_temperature + C_to_K, 'Q', 1, self.fluid_name)
                self.model.fs.P_low_target.set_value(P_low)
                self.model.fs.P_low.set_value(P_low)
            if (ambient_temperature is not None) and (condenser_approach is not None):
                T_cond_sat = ambient_temperature + condenser_approach + C_to_K
                P_high = CP.PropsSI('P', 'T', T_cond_sat, 'Q', 1, self.fluid_name)
                self.model.fs.P_high_target.set_value(P_high)
                self.model.fs.P_high.set_value(P_high)

            for arc in [
                self.model.fs.evaporator_to_compressor_expanded,
                self.model.fs.compressor_to_condenser_expanded,
                self.model.fs.condenser_to_expansion_valve_expanded,
                self.model.fs.expansion_valve_to_evaporator_expanded,
            ]:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.deactivate()

            # Activate explicit pressure-level constraints
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
            # Keep arc pressure equalities; deactivate explicit pressure-level constraints
            for arc in [
                self.model.fs.evaporator_to_compressor_expanded,
                self.model.fs.compressor_to_condenser_expanded,
                self.model.fs.condenser_to_expansion_valve_expanded,
                self.model.fs.expansion_valve_to_evaporator_expanded,
            ]:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.activate()
            self.model.fs.P_low_target_constraint.deactivate()
            self.model.fs.P_high_target_constraint.deactivate()
            self.model.fs.P_low_evap_in.deactivate()
            self.model.fs.P_low_evap_out.deactivate()
            self.model.fs.P_low_comp_in.deactivate()
            self.model.fs.P_low_valve_out.deactivate()
            self.model.fs.P_high_comp_out.deactivate()
            self.model.fs.P_high_cond_in.deactivate()
            self.model.fs.P_high_cond_out.deactivate()
            self.model.fs.P_high_valve_in.deactivate()

        # Expansion valve pressure bounds
        if low_side_pressure_min:
            self.model.fs.expansion_valve.outlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.expansion_valve.outlet.pressure[0].setub(low_side_pressure_max)
        if high_side_pressure_min:
            self.model.fs.expansion_valve.inlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.expansion_valve.inlet.pressure[0].setub(high_side_pressure_max)

        # Expansion valve temperature bounds (outlet)
        # TODO: Write constraint with control volume if PH mode
        if check_input(expansion_valve_temperature):
            if expansion_valve_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.expansion_valve.Tmin.set_value(expansion_valve_temperature[0] + C_to_K)
                    self.model.fs.expansion_valve.T_lower_bound.activate()
                else:
                    self.model.fs.expansion_valve.outlet.temperature[0].setlb(expansion_valve_temperature[0]+ C_to_K)
            if expansion_valve_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.expansion_valve.Tmax.set_value(expansion_valve_temperature[1] + C_to_K)
                    self.model.fs.expansion_valve.T_upper_bound.activate()
                else:
                    self.model.fs.expansion_valve.outlet.temperature[0].setub(expansion_valve_temperature[1]+ C_to_K)

        # Expansion valve outlet must be two-phase
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.expansion_valve.control_volume.properties_out[0].phase_frac["Vap"].setlb(0.01)
            self.model.fs.expansion_valve.control_volume.properties_out[0].phase_frac["Vap"].setub(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            # Phase 3b: eq_complementarity / eq_sat are Helmholtz-only; the generic
            # SmoothVLE handles the phase transition -- nothing to toggle here.
            pass

        # Phase 3b: re-assert the sum_mole_frac_out deactivation from __init__.
        # Confirmed via diagnostic: each unit's OWN initialize() call reactivates
        # its own properties_out[0.0].sum_mole_frac_out as part of its internal
        # bootstrapping (evaporator and compressor came back active=True after
        # vc.initialize(), even though __init__ deactivated all 4 -- condenser and
        # expansion_valve only stayed deactivated because initialize() failed
        # before reaching them). set_specifications() runs after initialize() and
        # is the last thing before the real solve, so re-deactivating here, right
        # before calculate_scaling_factors, is the one place guaranteed to stick.
        for unit_name in ["evaporator", "compressor", "condenser", "expansion_valve"]:
            unit = getattr(self.model.fs, unit_name)
            unit.control_volume.properties_out[0.0].sum_mole_frac_out.deactivate()

        # Calculate scaling factors
        calculate_scaling_factors(self.model)

    def optimize_COP(self, verbose, initialize=True, optimize=True, solver_options=None):
        """
        Args:
            solver_options: optional dict to OVERRIDE/EXTEND the default
                Ipopt options below (added 2026-07-24 so per-method scripts,
                e.g. the SPGP sweep, can try a more permissive tolerance
                without duplicating this whole method). If None, behavior
                is unchanged from before this parameter existed.
        """

        solver = get_solver()
        # FIX (2026-07-23, night): added constr_viol_tol/acceptable_tol.
        # Confirmed via coupled_solve_debug.py that the real coupled solve,
        # at the actual target spec (SH=SC=0), reaches a point with ZERO
        # large constraint residuals but still gets reported as "infeasible"
        # with the strict default tol=1e-6 alone -- many log_mole_frac_tbub/
        # tdew variables sit EXACTLY at 0.0 against a (None,0) bound, which
        # is correct for a pure fluid (log(1.0)=0) but numerically
        # degenerate right at the solution, tripping Ipopt's dual-
        # feasibility check even with a perfect primal solution. Same fix
        # already confirmed for the condenser and valve outlet earlier.
        default_options = {'max_iter': 1000, 'tol': 1e-4,
                            'constr_viol_tol': 1e-4, 'acceptable_tol': 1e-3}
        if solver_options:
            default_options.update(solver_options)
        solver.options = default_options
        
        if initialize:
            self.logger.info("Initializing the flowsheet by solving with no objective...")

            # Disable the objective function
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()

            results = solver.solve(self.model, tee=verbose)

            if results.solver.termination_condition == "optimal":
                self.logger.info("Initialization successful")
            else:
                self.logger.error("Initialization failed")

            if verbose:
                self.model.fs.report()

            # Compute the COP
            self.model.fs.cop.set_value(
                self.model.fs.evaporator.heat_duty[0].value
                / self.model.fs.compressor.work_mechanical[0].value
            )

        self.logger.info("Setting up the optimization problem...")

        if optimize:
            # Activate COP constraint and objective
            self.model.fs.compute_cop.activate()
            self.model.fs.obj.activate()
        else:
            # Feasibility solve only; compute COP after
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()
            
        # Solve the problem
        results = solver.solve(self.model, tee=verbose)

        if not optimize:
            if (self.model.fs.compressor.work_mechanical[0].value is not None
                and self.model.fs.compressor.work_mechanical[0].value != 0):
                self.model.fs.cop.set_value(
                    self.model.fs.evaporator.heat_duty[0].value
                    / self.model.fs.compressor.work_mechanical[0].value
                )

        # Resolve if the optimizer got stuck
        if results.solver.termination_condition != "optimal":
            results = solver.solve(self.model, tee=verbose)

        # Resolve if the optimizer got stuck
        if results.solver.termination_condition != "optimal":
            results = solver.solve(self.model, tee=verbose)

        # Check the solver status
        if results.solver.termination_condition == "optimal":    
            self.logger.info("Optimization successful")
            self.logger.info("COP: {:.2f}".format(value(self.model.fs.cop)))
            optimization_converged = True
        else:
            self.logger.error("Optimization failed")

            # Create a diagnostics toolbox instance
            diag = DiagnosticsToolbox(self.model, constraint_residual_tolerance=1e-6)

            diag.display_constraints_with_large_residuals()
            optimization_converged = False

        if verbose:
            self.model.fs.report()

        # Save the status
        self.optimization_converged = optimization_converged

        return value(self.model.fs.cop), optimization_converged

    def report_solution(self):

        print("Optimized COP:", round(value(self.model.fs.cop),3))

        n = len(self.h_init)

        h_sol = np.zeros(n)
        p_sol = np.zeros(n)
        T_sol = np.zeros(n)
        S_sol = np.zeros(n)

        for i, unit in enumerate(self.unit_operations):
            h_sol[i] = unit.control_volume.properties_out[0].enth_mol()
            p_sol[i] = unit.outlet.pressure[0].value
            T_sol[i] = unit.control_volume.properties_out[0].temperature()
            S_sol[i] = unit.control_volume.properties_out[0].entr_mol()

        def add_warning():
            if self.optimization_converged == None:
                # Have not run the optimization yet
                pass
            elif not self.optimization_converged:
                xlim = plt.gca().get_xlim()
                ylim = plt.gca().get_ylim()
                x = xlim[1] - (xlim[1] - xlim[0]) * 0.1
                y = ylim[0] + (ylim[1] - ylim[0]) * 0.1
                plt.text(x, y, "Warning: did not converge", color="red", fontsize=12, bbox=dict(facecolor='white', alpha=0.8), va='bottom', ha='right')

        # Phase 3b: Helmholtz-only diagrams removed (not on the generic PR package).

        for unit in self.unit_operations:
            unit.report()
