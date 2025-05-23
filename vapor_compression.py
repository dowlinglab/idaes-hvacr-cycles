# Import required IDAES-PSE modules
from idaes.core import FlowsheetBlock
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
)
from idaes.models.unit_models import (Heater, Turbine, Compressor, 
                                      Mixer, Separator, PressureChanger,
                                      Valve)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
from pyomo.environ import ConcreteModel, value, Objective, SolverFactory, maximize, minimize, TransformationFactory, Param, Var
from idaes.core.util.initialization import propagate_state
from pyomo.network import Arc
from pyomo.environ import units as pyunits
import numpy as np
import matplotlib.pyplot as plt

import logging
import CoolProp.CoolProp as CP

from idaes.core.util import DiagnosticsToolbox
from enum import Enum

class Mode(Enum):
    ORIGINAL_TPX = "original_TPx"
    IMPROVED_TPX = "improved_TPx"
    PH = "PH"

# Conversion factor
C_to_K = 273.15

class SimpleVaporCompressionCycle:



    def __init__(self,fluid_name, compressor_efficiency=0.85, mode=Mode.IMPROVED_TPX):
        ''' Simple Vapor Compression Cycle

        Parameters:
            fluid_name : str
                Name of the fluid
            compressor_efficiency : float
                Isentropic efficiency of the compressor
        
        '''

        self.fluid_name = fluid_name

        # Create the ConcreteModel and Flowsheet
        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        # State variables
        if mode == Mode.PH:
            sv = StateVars.PH
        else:
            sv = StateVars.TPX

        # Save the mode
        self.mode = mode

        # Add property package for a common refrigerant (e.g., R134a) using General Helmholtz model
        self.model.fs.properties = HelmholtzParameterBlock(pure_component=fluid_name, 
                                                    state_vars=sv, 
                                                    amount_basis=AmountBasis.MASS)
        
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
        self.model.fs.evaporator_to_compressor_expanded.flow_mass_equality.deactivate()

        # Let's see if this helps with convergence
        # self.model.fs.evaporator_to_compressor_expanded.pressure_equality.deactivate()
        # self.model.fs.expansion_valve_to_evaporator_expanded.pressure_equality.deactivate()
        # self.model.fs.expansion_valve_to_evaporator_expanded.temperature_equality.deactivate()

        # Set up the objective function
        '''
        self.model.fs.COP = Objective(expr=(self.model.fs.evaporator.heat_duty[0]) /
                                (self.model.fs.compressor.work_mechanical[0]), sense=maximize)
        '''
        self.model.fs.cop = Var(initialize=1, units=pyunits.dimensionless, bounds=(0.1, 10))

        @self.model.fs.Constraint(doc="COP constraint")
        def compute_cop(b):
            return b.cop * b.compressor.work_mechanical[0] == b.evaporator.heat_duty[0] 
        
        @self.model.fs.Objective(doc="Maximize COP", sense=maximize)
        def obj(b):
            return b.cop
                                
        self.model.fs.compute_cop.deactivate()
        self.model.fs.obj.deactivate()
        
        self.model.fs.evaporator.superheating = Param(initialize=0, units=pyunits.K, mutable=True)

        @self.model.fs.evaporator.Constraint(doc="Superheat evaporator outlet")
        def superheating_constraint(b):
            return b.control_volume.properties_out[0].temperature >= b.control_volume.properties_out[0].temperature_sat + b.superheating
        
        self.model.fs.evaporator.superheating_constraint.deactivate()

        self.model.fs.condenser.subcooling = Param(initialize=0, units=pyunits.K, mutable=True)

        @self.model.fs.condenser.Constraint(doc="Subcool condenser outlet") 
        def subcooling_constraint(b):
            return b.control_volume.properties_out[0].temperature <= b.control_volume.properties_out[0].temperature_sat - b.subcooling
        
        self.model.fs.condenser.subcooling_constraint.deactivate()

        @self.model.fs.compressor.Constraint(doc="Must be a vapor")
        def vapor_constraint(b):
            return b.control_volume.properties_out[0].temperature >= b.control_volume.properties_out[0].temperature_sat
            # return b.outlet.pressure[0]/1e3 >= b.control_volume.properties_out[0].pressure_sat/1e3

        self.model.fs.compressor.vapor_constraint.deactivate()

        '''
        # This constraint is redundant with another constraint
        @self.model.fs.expansion_valve.Constraint(doc="Must be two-phase")
        def two_phase_constraint(b):
            # return b.control_volume.properties_out[0].temperature[0] == b.control_volume.properties_out[0].temperature_sat
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
        self.model.fs.properties.hp_diagram()
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.show()

        self.model.fs.properties.ts_diagram()
        plt.show()

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

        self.model.fs.properties.hp_diagram()
        plt.plot(self.h_init/1000, self.p_init/1000, 'ko')
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.plot(self.T_init, self.p_init/1000, 'ko')
        plt.show()

        self.S_init = np.zeros(len(self.h_init))
        for i in range(len(self.h_init)):
            # Inputs are in J/kg and Pa, hence 1000 is needed for unit conversion
            # Output has units J/kg.K, hence 1000 is needed for unit conversion
            self.S_init[i] = CP.PropsSI('S', 'H', self.h_init[i], 'P', self.p_init[i], self.fluid_name)

        self.model.fs.properties.ts_diagram()
        plt.plot(self.S_init/1000, self.T_init, 'ko')
        plt.show()  

    def initialize(self, verbose=False):
        ''' Initialize the flowsheet '''


        # p_init has units of Pa, hence the scale factor is 1 to get Pa
        p_scale = 1

        ## Evaporator
        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)   # Example value

        self.model.fs.evaporator.inlet.pressure[0].fix(self.p_init[-1]*p_scale)

        if self.mode == Mode.PH:
            # Initialize H
            self.model.fs.evaporator.inlet.enth_mass[0].fix(self.h_init[-1])
            self.model.fs.evaporator.outlet.enth_mass[0].fix(self.h_init[0])
        else:
            self.model.fs.evaporator.inlet.temperature[0].fix(self.T_init[-1])
            self.model.fs.evaporator.outlet.temperature[0].fix(self.T_init[0])

        self.logger.info("Initializing evaporator...")
        self.model.fs.evaporator.initialize(outlvl=logging.WARNING)

        if verbose:
            self.model.fs.evaporator.report()

        propagate_state(self.model.fs.evaporator_to_compressor)

        ## Compressor
        self.model.fs.compressor.inlet.pressure[0].fix(self.p_init[0]*p_scale)

        if self.mode == Mode.PH:
            self.model.fs.compressor.inlet.enth_mass[0].fix(self.h_init[0])
        else:
            self.model.fs.compressor.inlet.temperature[0].fix(self.T_init[0])

        self.model.fs.compressor.outlet.pressure[0].fix(self.p_init[1]*p_scale) # Set to target
        # self.model.fs.compressor.outlet.vapor_frac[0].fix(1.0)  # Ensure vapor phase

        self.model.fs.compressor.efficiency_isentropic.fix(self.compressor_efficiency)

        self.logger.info("Initializing compressor...")
        self.model.fs.compressor.initialize(outlvl=logging.WARNING)

        if verbose:
            self.model.fs.compressor.report()

        propagate_state(self.model.fs.compressor_to_condenser)

        ## Condenser
        
        self.model.fs.condenser.inlet.pressure[0].fix(self.p_init[1]*p_scale)

        if self.mode == Mode.PH:
            self.model.fs.condenser.inlet.enth_mass[0].fix(self.h_init[1])
            self.model.fs.condenser.outlet.enth_mass[0].fix(self.h_init[2])
        else:
            self.model.fs.condenser.inlet.temperature[0].fix(self.T_init[1])
            self.model.fs.condenser.outlet.temperature[0].fix(self.T_init[2])

            # self.model.fs.condenser.outlet.vapor_frac[0].fix(0.0)  # Ensure liquid phase (is this needed?)

        self.logger.info("Initializing condenser...")
        self.model.fs.condenser.initialize(outlvl=logging.WARNING)
        if verbose:
            self.model.fs.condenser.report()

        propagate_state(self.model.fs.condenser_to_expansion_valve)

        ## Expansion Valve
        self.model.fs.expansion_valve.inlet.pressure[0].fix(self.p_init[2]*p_scale)
        self.model.fs.expansion_valve.outlet.pressure[0].fix(self.p_init[3]*p_scale)

        self.logger.info("Initializing expansion valve...")
        self.model.fs.expansion_valve.initialize(outlvl=logging.WARNING)

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
                           max_pressure_ratio = 4):
        
        assert superheating >= 0, "Superheating must be greater than or equal to 0"
        assert subcooling >= 0, "Subcooling must be greater than or equal to 0"
        assert max_pressure_ratio > 1.2, "Maximum pressure ratio must be greater than 1.2"
        # assert not (Tsat_constraints and bound_vapor_frac), "Cannot have both Tsat constraints and vapor fraction bounds"

        ## Unfix the conditions from initialization
        self.unit_operations = [self.model.fs.evaporator, self.model.fs.compressor, self.model.fs.condenser, self.model.fs.expansion_valve]

        for unit in self.unit_operations:
            # Unfix all variables
            unit.inlet.flow_mass[0].unfix()
            unit.outlet.flow_mass[0].unfix()

            if self.mode == Mode.PH:
                unit.inlet.enth_mass[0].unfix()
                unit.outlet.enth_mass[0].unfix()
            else:
                unit.inlet.temperature[0].unfix()
                unit.inlet.vapor_frac[0].unfix()
                unit.outlet.temperature[0].unfix()
                unit.outlet.vapor_frac[0].unfix()
            
            unit.inlet.pressure[0].unfix()
            unit.outlet.pressure[0].unfix()

            unit.T_lower_bound.deactivate()
            unit.T_upper_bound.deactivate()

            # Set bounds for the vapor fraction to ensure it is within [0,1]
            if self.mode == Mode.ORIGINAL_TPX or self.mode == Mode.IMPROVED_TPX:
                unit.inlet.vapor_frac[0].setlb(0)
                unit.inlet.vapor_frac[0].setub(1)

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

        # Convert pressures from kPa to Pa
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
        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)

        ## Evaporator

        # Evaporator pressure bounds
        if low_side_pressure_min:
            self.model.fs.evaporator.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.evaporator.inlet.pressure[0].setub(low_side_pressure_max)

        # Evaporator temperature bounds (outlet)
        ## TODO: Set these bounds on the control variable if PH mode
        if check_input(evaporator_temperature):
            if evaporator_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.evaporator.Tmin = evaporator_temperature[0] + C_to_K
                    self.model.fs.evaporator.T_lower_bound.activate()

                else:
                    self.model.fs.evaporator.outlet.temperature[0].setlb(evaporator_temperature[0] + C_to_K)
            
            if evaporator_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.evaporator.Tmax = evaporator_temperature[1] + C_to_K
                    self.model.fs.evaporator.T_upper_bound.activate()
                else:
                    self.model.fs.evaporator.outlet.temperature[0].setub(evaporator_temperature[1] + C_to_K)

        # Evaporator outlet must be a vapor
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.evaporator.outlet.vapor_frac[0].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.evaporator.outlet.vapor_frac[0].fix(1.0)

        # Deactivate the complementarity-like constraint because we fixed the phase
        if self.mode == Mode.IMPROVED_TPX:
            self.model.fs.evaporator.control_volume.properties_out[0.0].eq_complementarity.deactivate()

        # Activate superheating constraint
        if superheating > 0.1:
            self.model.fs.evaporator.superheating = superheating
            self.model.fs.evaporator.superheating_constraint.activate()
        else:
            self.model.fs.evaporator.superheating_constraint.deactivate()

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
                    self.model.fs.compressor.Tmin = compressor_temperature[0] + C_to_K
                    self.model.fs.compressor.T_lower_bound.activate()
                else:
                    self.model.fs.compressor.outlet.temperature[0].setlb(compressor_temperature[0] + C_to_K)

            if compressor_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.compressor.Tmax = compressor_temperature[1] + C_to_K
                    self.model.fs.compressor.T_upper_bound.activate()
                else:
                    self.model.fs.compressor.outlet.temperature[0].setub(compressor_temperature[1] + C_to_K)
        
        # Compressor outlet must be a vapor
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.compressor.outlet.vapor_frac[0].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.compressor.outlet.vapor_frac[0].fix(1.0)

            # Deactivate the complementarity-like constraint because we fixed the phase
            self.model.fs.compressor.control_volume.properties_out[0.0].eq_complementarity.deactivate()

            # Add inequality constraint to ensure the compressor outlet is only vapor
            self.model.fs.compressor.vapor_constraint.activate()

        # Compressor only allows input work
        # self.model.fs.compressor.work_mechanical.setlb(0)

        # Set the maximum pressure ratio
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
        if check_input(condenser_temperature):
            if condenser_temperature[0]:
                if self.mode == Mode.PH:
                    self.model.fs.condenser.Tmin = condenser_temperature[0] + C_to_K
                    self.model.fs.condenser.T_lower_bound.activate()
                else:
                    self.model.fs.condenser.outlet.temperature[0].setlb(condenser_temperature[0]+ C_to_K)
            if condenser_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.condenser.Tmax = condenser_temperature[1] + C_to_K
                    self.model.fs.condenser.T_upper_bound.activate()
                else:
                    self.model.fs.condenser.outlet.temperature[0].setub(condenser_temperature[1]+ C_to_K)

        # Condenser outlet must be a liquid
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.condenser.outlet.vapor_frac[0].setub(0.01)
        elif self.mode == Mode.IMPROVED_TPX:
            self.model.fs.condenser.outlet.vapor_frac[0].fix(0.0)

            # Deactivate the complementarity-like constraint because we fixed the phase
            self.model.fs.condenser.control_volume.properties_out[0.0].eq_complementarity.deactivate()

        # Activate subcooling constraint
        if subcooling > 0.1:
            self.model.fs.condenser.subcooling = subcooling
            self.model.fs.condenser.subcooling_constraint.activate()
        else:
            self.model.fs.condenser.subcooling_constraint.deactivate()


        ## Expansion Valve

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
                    self.model.fs.expansion_valve.Tmin = expansion_valve_temperature[0] + C_to_K
                    self.model.fs.expansion_valve.T_lower_bound.activate()
                else:
                    self.model.fs.expansion_valve.outlet.temperature[0].setlb(expansion_valve_temperature[0]+ C_to_K)
            if expansion_valve_temperature[1]:
                if self.mode == Mode.PH:
                    self.model.fs.expansion_valve.Tmax = expansion_valve_temperature[1] + C_to_K
                    self.model.fs.expansion_valve.T_upper_bound.activate()
                else:
                    self.model.fs.expansion_valve.outlet.temperature[0].setub(expansion_valve_temperature[1]+ C_to_K)

        # Expansion valve outlet must be two-phase
        if self.mode == Mode.ORIGINAL_TPX:
            self.model.fs.expansion_valve.outlet.vapor_frac[0].setlb(0.01)
            self.model.fs.expansion_valve.outlet.vapor_frac[0].setub(0.99)
        elif self.mode == Mode.IMPROVED_TPX:

            # This constraint has two smoothed max operators, might have a point singularity
            self.model.fs.expansion_valve.control_volume.properties_out[0.0].eq_complementarity.deactivate()

            # Use this constraint instead
            self.model.fs.expansion_valve.control_volume.properties_out[0.0].eq_sat.activate()

        # Calculate scaling factors
        calculate_scaling_factors(self.model)

    def optimize_COP(self, verbose, initialize=True):

        solver = get_solver()
        solver.options = {'max_iter': 1000, 'tol': 1e-6, 'linear_solver':'ma57'} # relax the tolerance
        
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
            self.model.fs.cop = self.model.fs.evaporator.heat_duty[0].value / self.model.fs.compressor.work_mechanical[0].value

        self.logger.info("Setting up the optimization problem...")

        # Activate the objective function
        self.model.fs.compute_cop.activate()
        self.model.fs.obj.activate()
            
        # Solve the optimization problem
        results = solver.solve(self.model, tee=verbose)

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
            h_sol[i] = unit.control_volume.properties_out[0].enth_mass()
            p_sol[i] = unit.outlet.pressure[0].value
            T_sol[i] = unit.control_volume.properties_out[0].temperature()
            S_sol[i] = unit.control_volume.properties_out[0].entr_mass()

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

        self.model.fs.properties.hp_diagram()
        plt.plot(h_sol/1000, p_sol/1000, 'ko')
        add_warning()
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.plot(T_sol, p_sol/1000, 'ko')
        add_warning()
        plt.show()

        self.model.fs.properties.ts_diagram()
        plt.plot(S_sol/1000, T_sol, 'ko')
        add_warning()
        plt.show()

        for unit in self.unit_operations:
            unit.report()
