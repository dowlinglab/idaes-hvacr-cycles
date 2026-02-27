import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from pyomo.environ import value

from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode


def _disable_plots():
    try:
        plt.show = lambda *args, **kwargs: None
    except Exception:
        pass


def main():
    _disable_plots()

    design_kwargs = dict(
        ambient_temperature=35,
        condenser_approach=5,
        evap_sat_temperature=-10,
        superheating=3,
        subcooling=3,
        max_pressure_ratio=4,
        debug_disable_arc_pressure_eq=True,
    )

    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        mode=Mode.IMPROVED_TPX,
    )

    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    Q_rated = cycle.rate_capacity_from_design(design_kwargs, verbose=False)
    print(f"Q_rated: {Q_rated}")

    cycle.set_specifications(plr=0.5, Q_cool_rated=Q_rated, **design_kwargs)
    cycle.optimize_COP(verbose=False, initialize=True, optimize=False)

    cop = SimpleVaporCompressionCyclePLR.compute_model_cop(cycle.model)
    mdot = value(cycle.model.fs.evaporator.inlet.flow_mass[0])
    print(f"PLR=0.5 COP_model: {cop}")
    print(f"PLR=0.5 mdot: {mdot}")

    TL_K = -10 + 273.15
    TH_K = 35 + 5 + 273.15
    cop_carnot = SimpleVaporCompressionCyclePLR.compute_carnot_cop_ref(TL_K, TH_K)
    print(f"Carnot COP: {cop_carnot}")
    print(f"COP_model / COP_carnot: {cop / cop_carnot}")


if __name__ == "__main__":
    main()
