# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Cycle configuration for the HX-ideal variant.

Context
-------
This module belongs to the `cycle_dx_hx_ideal` variant. Relative to baseline
intent, this variant preserves idealized heat exchanger behavior and only
provides a standardized configuration object and API-compatible parameters.
"""

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class CycleConfig:
    """Container for cycle-level and heat-exchanger parameters.

    Attributes
    ----------
    eta_isentropic:
        Compressor isentropic efficiency (fraction).
    m_dot_ref:
        Refrigerant mass flow rate [kg/s].
    m_dot_air_evap:
        Evaporator-side air mass flow [kg/s].
    m_dot_air_cond:
        Condenser-side air mass flow [kg/s].
    cp_air_evap:
        Evaporator-side air heat capacity [J/(kg-K)].
    cp_air_cond:
        Condenser-side air heat capacity [J/(kg-K)].
    superheat_target_K:
        Evaporator outlet superheat target above saturation [K].
    subcool_target_K:
        Condenser outlet subcool target below saturation [K].
    UA_evap_total:
        Lumped evaporator UA [W/K]. Retained for interface consistency.
    UA_cond_total:
        Lumped condenser UA [W/K]. Retained for interface consistency.
    UA_evap_tp:
        Evaporator two-phase zone UA [W/K]. Retained for interface consistency.
    UA_evap_sh:
        Evaporator superheat zone UA [W/K]. Retained for interface consistency.
    UA_cond_ds:
        Condenser desuperheat zone UA [W/K]. Retained for interface consistency.
    UA_cond_tp:
        Condenser two-phase zone UA [W/K]. Retained for interface consistency.
    UA_cond_sc:
        Condenser subcool zone UA [W/K]. Retained for interface consistency.
    hx_arrangement:
        HX flow arrangement string used by epsilon-NTU helpers.
    """

    eta_isentropic: float = 0.75
    m_dot_ref: float = 0.02
    plr: float = 0.75
    cd: float = 0.13

    m_dot_air_evap: float = 1.2
    m_dot_air_cond: float = 1.5
    cp_air_evap: float = 1006.0
    cp_air_cond: float = 1006.0

    superheat_target_K: float = 5.0
    subcool_target_K: float = 5.0

    UA_evap_total: float = 2500.0
    UA_cond_total: float = 3000.0

    UA_evap_tp: float = 1800.0
    UA_evap_sh: float = 700.0

    UA_cond_ds: float = 800.0
    UA_cond_tp: float = 1700.0
    UA_cond_sc: float = 500.0

    hx_arrangement: str = "counterflow"

    def scaled_ua(self, factor: float) -> "CycleConfig":
        """Return a new config with all UA values scaled by ``factor``."""
        return replace(
            self,
            UA_evap_total=self.UA_evap_total * factor,
            UA_cond_total=self.UA_cond_total * factor,
            UA_evap_tp=self.UA_evap_tp * factor,
            UA_evap_sh=self.UA_evap_sh * factor,
            UA_cond_ds=self.UA_cond_ds * factor,
            UA_cond_tp=self.UA_cond_tp * factor,
            UA_cond_sc=self.UA_cond_sc * factor,
        )
