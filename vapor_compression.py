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
from pyomo.environ import ConcreteModel, value, Objective, SolverFactory, maximize, minimize, TransformationFactory
from idaes.core.util.initialization import propagate_state
from pyomo.network import Arc
from pyomo.environ import units as pyunits
import numpy as np
import matplotlib.pyplot as plt

import logging
import CoolProp.CoolProp as CP

from idaes.core.util import DiagnosticsToolbox

class SimpleVaporCompressionCycle:

    def __init__(self,fluid_name, compressor_efficiency=0.85):
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

        # Add property package for a common refrigerant (e.g., R134a) using General Helmholtz model
        self.model.fs.properties = HelmholtzParameterBlock(pure_component="R134a", 
                                                    state_vars=StateVars.TPX, 
                                                    amount_basis=AmountBasis.MASS)
        
        # Save the compressor efficiency
        assert 0 < compressor_efficiency < 1, "Compressor efficiency must be between 0 and 1"
        self.compressor_efficiency = compressor_efficiency
        
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

        # Set up the objective function
        self.model.fs.COP = Objective(expr=(self.model.fs.evaporator.heat_duty[0]) /
                                (self.model.fs.compressor.work_mechanical[0]), sense=maximize)
        
        self.model.fs.COP.deactivate()

    def draw_thermodynamic_diagrams(self):
        self.model.fs.properties.hp_diagram()
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.show()

        self.model.fs.properties.ts_diagram()
        plt.show()

    # TODO: Make this a private method
    # Instead, the user should just specify the high and low temperatures
    # From the PT diagram, we can determine the high and low pressures
    # Then from the PH diagram, we can determine the enthalpy values
    # And call this function
    def specify_initial_conditions(self,
                                   enthalpy = [410, 430, 256, 260], # kJ/kg
                                   pressure = [300, 1100, 1100, 300]): # kPa
    
        ''' Specify initial conditions for the flowsheet

        Components:
            0: Evaporator
            1: Compressor
            2: Condenser
            3: Expansion Valve

        Parameters:
            enthalpy : list
                Initial enthalpy values for each component in the cycle
            pressure : list
                Initial pressure values for each component in the cycle
        
        '''

        self.h_init = enthalpy
        self.p_init = pressure

        self.T_init = np.zeros(len(self.h_init))
        for i in range(len(self.h_init)):
            # Inputs are in J/kg and Pa, hence 1000 is needed for unit conversion
            self.T_init[i] = CP.PropsSI('T', 'H', self.h_init[i]*1000, 'P', self.p_init[i]*1000, self.fluid_name)

        self.model.fs.properties.hp_diagram()
        plt.plot(self.h_init, self.p_init, 'ko')
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.plot(self.T_init, self.p_init, 'ko')
        plt.show()

        self.S_init = np.zeros(len(self.h_init))
        for i in range(len(self.h_init)):
            # Inputs are in J/kg and Pa, hence 1000 is needed for unit conversion
            # Output has units J/kg.K, hence 1000 is needed for unit conversion
            self.S_init[i] = CP.PropsSI('S', 'H', self.h_init[i]*1000, 'P', self.p_init[i]*1000, self.fluid_name)/1000

        self.model.fs.properties.ts_diagram()
        plt.plot(self.S_init, self.T_init, 'ko')
        plt.show()  

    def initialize(self, verbose=False):
        ''' Initialize the flowsheet '''


        # p_init has units of kPa, hence the scale factor is 1e3 to get Pa
        p_scale = 1e3

        ## Evaporator
        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)   # Example value
        self.model.fs.evaporator.inlet.temperature[0].fix(self.T_init[-1])
        self.model.fs.evaporator.inlet.pressure[0].fix(self.p_init[-1]*p_scale)

        self.model.fs.evaporator.outlet.temperature[0].fix(self.T_init[0])

        self.logger.info("Initializing evaporator...")
        self.model.fs.evaporator.initialize(outlvl=logging.WARNING)

        if verbose:
            self.model.fs.evaporator.report()

        propagate_state(self.model.fs.evaporator_to_compressor)

        ## Compressor
        self.model.fs.compressor.inlet.temperature[0].fix(self.T_init[0])
        self.model.fs.compressor.inlet.pressure[0].fix(self.p_init[0]*p_scale)

        self.model.fs.compressor.outlet.pressure[0].fix(self.p_init[1]*p_scale) # Set to target
        # self.model.fs.compressor.outlet.vapor_frac[0].fix(1.0)  # Ensure vapor phase

        self.model.fs.compressor.efficiency_isentropic.fix(self.compressor_efficiency)

        self.logger.info("Initializing compressor...")
        self.model.fs.compressor.initialize(outlvl=logging.WARNING)

        if verbose:
            self.model.fs.compressor.report()

        propagate_state(self.model.fs.compressor_to_condenser)

        self.model.fs.condenser.inlet.temperature[0].fix(self.T_init[1])
        self.model.fs.condenser.inlet.pressure[0].fix(self.p_init[1]*p_scale)

        self.model.fs.condenser.outlet.vapor_frac[0].fix(0.0)  # Ensure liquid phase

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
                           low_side_pressure = (2E5, 5E5),
                           high_side_pressure = (1E6, 3E6),
                           evaporator_temperature = (-20, 0),
                           compressor_temperature = None,
                           condenser_temperature = (30, 50),
                           expansion_valve_temperature = None,
                           subcooling = 3, # degC
                           superheating = 3 # degC
                           ):
        
        C_to_K = 273.15

        ## Unfix the conditions from initialization
        self.unit_operations = [self.model.fs.evaporator, self.model.fs.compressor, self.model.fs.condenser, self.model.fs.expansion_valve]

        for unit in self.unit_operations:
            # Unfix all variables
            unit.inlet.flow_mass[0].unfix()
            unit.inlet.temperature[0].unfix()
            unit.inlet.pressure[0].unfix()
            unit.outlet.flow_mass[0].unfix()
            unit.outlet.temperature[0].unfix()
            unit.outlet.pressure[0].unfix()

            # Set bounds for the vapor fraction to ensure it is within [0,1]
            unit.inlet.vapor_frac[0].setlb(0)
            unit.inlet.vapor_frac[0].setub(1)

        # Convert pressures from kPa to Pa
        if not None and len(low_side_pressure) == 2:
            low_side_pressure_min = low_side_pressure[0]*1000
            low_side_pressure_max = low_side_pressure[1]*1000
        else:
            low_side_pressure_min = None
            low_side_pressure_max = None

        if not None and len(high_side_pressure) == 2:
            high_side_pressure_min = high_side_pressure[0]*1000
            high_side_pressure_max = high_side_pressure[1]*1000
        else:
            high_side_pressure_min = None
            high_side_pressure_max = None

        # Set mass flowrate to 1 kg/s because we only care about thermodynamic efficiency
        self.model.fs.evaporator.inlet.flow_mass[0].fix(1)

        def check_input(bounds):
            if bounds is not None and len(bounds) == 2:
                return True
            else:
                return False

        ## Evaporator

        # Evaporator pressure bounds
        if low_side_pressure_min:
            self.model.fs.evaporator.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            self.model.fs.evaporator.inlet.pressure[0].setub(low_side_pressure_max)

        # Evaporator temperature bounds (outlet)
        if check_input(evaporator_temperature):
            if evaporator_temperature[0]:
                self.model.fs.evaporator.outlet.temperature[0].setlb(evaporator_temperature[0] + C_to_K)
            
            if evaporator_temperature[1]:
                self.model.fs.evaporator.outlet.temperature[0].setub(evaporator_temperature[1] + C_to_K)

        # Evaporator outlet must be a vapor
        self.model.fs.evaporator.outlet.vapor_frac[0].setlb(0.999)

        
        if superheating > 0.1:
            @self.model.fs.evaporator.Constraint(doc="Superheat evaporator outlet")
            def subcooling_constraint(b):
                return b.outlet.temperature[0] >= b.control_volume.properties_out[0].t_sat_func + superheating

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
        if check_input(compressor_temperature):
            if compressor_temperature[0]:
                self.model.fs.compressor.outlet.temperature[0].setlb(compressor_temperature[0] + C_to_K)
            if compressor_temperature[1]:
                self.model.fs.compressor.outlet.temperature[0].setub(compressor_temperature[1] + C_to_K)
        
        # Compressor outlet must be a vapor
        self.model.fs.compressor.outlet.vapor_frac[0].setlb(0.999)

        # Compressor only allows input work
        self.model.fs.compressor.work_mechanical.setlb(0)

        ## Condenser

        # Condenser pressure bounds
        if high_side_pressure_min:
            self.model.fs.condenser.inlet.pressure[0].setlb(high_side_pressure_min)
            self.model.fs.condenser.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            self.model.fs.condenser.inlet.pressure[0].setub(high_side_pressure_max)
            self.model.fs.condenser.outlet.pressure[0].setub(high_side_pressure_max)

        # Condenser temperature bounds (outlet)
        if check_input(condenser_temperature):
            if condenser_temperature[0]:
                self.model.fs.condenser.outlet.temperature[0].setlb(condenser_temperature[0]+ C_to_K)
            if condenser_temperature[1]:
                self.model.fs.condenser.outlet.temperature[0].setub(condenser_temperature[1]+ C_to_K)

        # Condenser outlet must be a liquid
        self.model.fs.condenser.outlet.vapor_frac[0].setub(0.001)

        if subcooling > 0.1:
            @self.model.fs.condenser.Constraint(doc="Subcool condenser outlet") 
            def subcooling_constraint(b):
                return b.outlet.temperature[0] <= b.control_volume.properties_out[0].t_sat_func() - subcooling


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
        if check_input(expansion_valve_temperature):
            if expansion_valve_temperature[0]:
                self.model.fs.expansion_valve.outlet.temperature[0].setlb(expansion_valve_temperature[0]+ C_to_K)
            if expansion_valve_temperature[1]:
                self.model.fs.expansion_valve.outlet.temperature[0].setub(expansion_valve_temperature[1]+ C_to_K)

        # Expansion valve outlet must be two-phase
        self.model.fs.expansion_valve.outlet.vapor_frac[0].setlb(0.001)
        self.model.fs.expansion_valve.outlet.vapor_frac[0].setub(0.999)

        # Calculate scaling factors
        calculate_scaling_factors(self.model)

    def optimize_COP(self, verbose, initialize=True):

        solver = get_solver()
        solver.options = {'max_iter': 500}
        
        if initialize:
            self.logger.info("Initializing the flowsheet by solving with no objective...")

            # Disable the objective function
            self.model.fs.COP.deactivate()

            results = solver.solve(self.model, tee=verbose)

            if results.solver.termination_condition == "optimal":
                self.logger.info("Initialization successful")
            else:
                self.logger.error("Initialization failed")

            if verbose:
                self.model.fs.report()

        self.logger.info("Setting up the optimization problem...")

        # Activate the objective function
        self.model.fs.COP.activate()
            
        # Solve the optimization problem
        results = solver.solve(self.model, tee=verbose)

        # Resolve just in case the optimizer got stuck
        results = solver.solve(self.model, tee=verbose)

        # Resolve if the optimizer got stuck
        if results.solver.termination_condition != "optimal":
            results = solver.solve(self.model, tee=verbose)

        # Check the solver status
        if results.solver.termination_condition == "optimal":    
            self.logger.info("Optimization successful")
            self.logger.info("COP: {:.2f}".format(value(self.model.fs.COP)))
        else:
            self.logger.error("Optimization failed")

            # Create a diagnostics toolbox instance
            diag = DiagnosticsToolbox(self.model, constraint_residual_tolerance=1e-6)

            diag.display_constraints_with_large_residuals()

        if verbose:
            self.model.fs.report()

        return value(self.model.fs.COP)

    def report_solution(self):

        print("Optimized COP:", round(value(self.model.fs.COP),3))

        n = len(self.h_init)

        h_sol = np.zeros(n)
        p_sol = np.zeros(n)
        T_sol = np.zeros(n)
        S_sol = np.zeros(n)

        for i, unit in enumerate(self.unit_operations):
            h_sol[i] = unit.control_volume.properties_out[0].enth_mass()
            p_sol[i] = unit.outlet.pressure[0].value
            T_sol[i] = unit.outlet.temperature[0].value
            S_sol[i] = unit.control_volume.properties_out[0].entr_mass()

        self.model.fs.properties.hp_diagram()
        plt.plot(h_sol/1000, p_sol/1000, 'ko')
        plt.show()

        self.model.fs.properties.pt_diagram()
        plt.plot(T_sol, p_sol/1000, 'ko')
        plt.show()

        self.model.fs.properties.ts_diagram()
        plt.plot(S_sol/1000, T_sol, 'ko')
        plt.show()

        for unit in self.unit_operations:
            unit.report()

