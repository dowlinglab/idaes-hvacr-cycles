"""
PLR + Epsilon-NTU Lumped HX Model (Copy)

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex 
QA/Testing: Shilpa Narasimhan

Description:
    Copy-only PLR model wrapper that reproduces the non-IDAES lumped epsilon-NTU
    cycle equations using the existing `cycle_dx_hx_lumpedUA` solver path.

Context Breadcrumb:
    This module intentionally aligns with the non-IDAES lumped-UA formulation
    so baseline runs can be compared one-to-one against legacy results before
    attempting tighter IDAES unit-model coupling.

Notes:
    - All QA and validation are the responsibility of Shilpa Narasimhan.
    - This file was generated with Codex assistance.
"""

from __future__ import annotations

from dataclasses import replace
from enum import Enum
from typing import Optional

from cycle_dx_hx_lumpedUA.config import CycleConfig
from cycle_dx_hx_lumpedUA.cycle_model import solve_cycle_point


class Mode(Enum):
    """Mode placeholder kept for compatibility with existing runner imports."""

    PH = "PH"
    IMPROVED_TPX = "improved_TPx"


class SimpleVaporCompressionCyclePLREpsNTU:
    """PLR wrapper using the non-IDAES lumped epsilon-NTU cycle model.

    The public API mirrors `vapor_compression_plr.py` usage in existing run
    scripts (`set_specifications`, `optimize_COP`, `get_*_cop`) so this can be
    swapped in as a copy-only baseline without editing original PLR code.
    """

    def __init__(
        self,
        fluid_name: str,
        compressor_efficiency: float = 0.75,
        PLR: float = 0.75,
        CD: float = 0.13,
        mode: Mode = Mode.PH,
        UA_evap_total: float = 1500.0,
        UA_cond_total: float = 1800.0,
        UA_scale: float = 1.0,
        m_dot_ref: float = 0.02,
        m_dot_air_evap: float = 1.2,
        m_dot_air_cond: float = 1.5,
        cp_air_evap: float = 1006.0,
        cp_air_cond: float = 1006.0,
    ):
        """Create an epsilon-NTU PLR cycle wrapper.

        Args:
            fluid_name: Refrigerant fluid string (CoolProp-compatible).
            compressor_efficiency: Isentropic efficiency used in compressor step.
            PLR: Part-load ratio in [0, 1].
            CD: Cycling degradation coefficient in [0, 1].
            mode: Placeholder enum for script compatibility.
            UA_evap_total: Baseline evaporator UA [W/K] from non-IDAES model.
            UA_cond_total: Baseline condenser UA [W/K] from non-IDAES model.
            UA_scale: Optional scalar applied to both UA values.
            m_dot_ref: Refrigerant mass flow [kg/s].
            m_dot_air_evap: Evaporator air mass flow [kg/s].
            m_dot_air_cond: Condenser air mass flow [kg/s].
            cp_air_evap: Evaporator-side air heat capacity [J/(kg-K)].
            cp_air_cond: Condenser-side air heat capacity [J/(kg-K)].
        """
        assert 0.0 < compressor_efficiency < 1.0
        assert 0.0 <= PLR <= 1.0
        assert 0.0 <= CD <= 1.0
        assert UA_scale > 0.0

        self.fluid_name = fluid_name
        self.mode = mode

        self._base_cfg = CycleConfig(
            eta_isentropic=compressor_efficiency,
            m_dot_ref=m_dot_ref,
            plr=PLR,
            cd=CD,
            m_dot_air_evap=m_dot_air_evap,
            m_dot_air_cond=m_dot_air_cond,
            cp_air_evap=cp_air_evap,
            cp_air_cond=cp_air_cond,
            UA_evap_total=UA_evap_total,
            UA_cond_total=UA_cond_total,
        )
        self._ua_scale = UA_scale
        self._cfg = self._base_cfg.scaled_ua(UA_scale)

        self._t_evap_sat_c: float = -30.0
        self._t_cond_sat_c: float = 35.0
        self._t_air_evap_in_c: float = -20.0
        self._t_air_cond_in_c: float = 25.0

        self._last_result = None
        self.optimization_converged = None

    def specify_initial_conditions(self, low_side_temperature: float = -20.0, high_side_temperature: float = 30.0):
        """Compatibility no-op for runners expecting IDAES initialization API."""
        _ = (low_side_temperature, high_side_temperature)

    def initialize(self, verbose: bool = False):
        """Compatibility no-op for runners expecting IDAES initialization API."""
        _ = verbose

    def set_specifications(
        self,
        low_side_pressure=(60.0, 200.0),
        high_side_pressure=(500.0, 4000.0),
        evaporator_temperature=(-55.0, -20.0),
        cold_storage_setpoint: Optional[float] = None,
        evap_offset_bounds=(-35.0, 0.0),
        condenser_temperature=(15.0, 80.0),
        subcooling=3.0,
        superheating=3.0,
        max_pressure_ratio=20.0,
        ambient_temperature: Optional[float] = None,
        condenser_approach: Optional[float] = None,
        evap_sat_temperature: Optional[float] = None,
        debug_disable_arc_pressure_eq=False,
        plr: Optional[float] = None,
        cd: Optional[float] = None,
        UA_evap_total: Optional[float] = None,
        UA_cond_total: Optional[float] = None,
        UA_scale: Optional[float] = None,
        **kwargs,
    ):
        """Store operating specs and map them to non-IDAES lumped solver inputs.

        Mapping rules:
            - `t_evap_sat_c`: explicit `evap_sat_temperature` if provided; else
              midpoint of `evaporator_temperature` bounds (to mirror `_plr` style
              bound-driven specification when no explicit saturation target is set).
            - `t_cond_sat_c`: `ambient_temperature + condenser_approach` when both
              are provided; else midpoint of `condenser_temperature`.
            - Air inlets:
                evaporator air inlet uses cold-storage setpoint when provided,
                otherwise `t_evap_sat_c + 10`;
                condenser air inlet uses ambient when provided, otherwise
                `t_cond_sat_c - 10`.
        """
        _ = (
            low_side_pressure,
            high_side_pressure,
            subcooling,
            superheating,
            max_pressure_ratio,
            debug_disable_arc_pressure_eq,
            kwargs,
            evap_offset_bounds,
        )

        if UA_scale is not None:
            assert UA_scale > 0.0
            self._ua_scale = float(UA_scale)

        cfg = self._base_cfg
        if plr is not None:
            assert 0.0 <= plr <= 1.0
            cfg = replace(cfg, plr=float(plr))
        if cd is not None:
            assert 0.0 <= cd <= 1.0
            cfg = replace(cfg, cd=float(cd))
        if UA_evap_total is not None:
            cfg = replace(cfg, UA_evap_total=float(UA_evap_total))
        if UA_cond_total is not None:
            cfg = replace(cfg, UA_cond_total=float(UA_cond_total))
        self._base_cfg = cfg
        self._cfg = cfg.scaled_ua(self._ua_scale)

        if evap_sat_temperature is not None:
            self._t_evap_sat_c = float(evap_sat_temperature)
        else:
            self._t_evap_sat_c = 0.5 * (float(evaporator_temperature[0]) + float(evaporator_temperature[1]))

        if (ambient_temperature is not None) and (condenser_approach is not None):
            self._t_cond_sat_c = float(ambient_temperature + condenser_approach)
        else:
            self._t_cond_sat_c = 0.5 * (float(condenser_temperature[0]) + float(condenser_temperature[1]))

        if cold_storage_setpoint is not None:
            self._t_air_evap_in_c = float(cold_storage_setpoint)
        else:
            self._t_air_evap_in_c = float(self._t_evap_sat_c + 10.0)

        if ambient_temperature is not None:
            self._t_air_cond_in_c = float(ambient_temperature)
        else:
            self._t_air_cond_in_c = float(self._t_cond_sat_c - 10.0)

    def optimize_COP(self, verbose=False, initialize=True, optimize=False):
        """Run non-IDAES lumped epsilon-NTU solve and return full-load COP."""
        _ = (verbose, initialize, optimize)
        try:
            self._last_result = solve_cycle_point(
                fluid=self.fluid_name,
                t_evap_sat_c=self._t_evap_sat_c,
                t_cond_sat_c=self._t_cond_sat_c,
                cfg=self._cfg,
                t_air_evap_in_c=self._t_air_evap_in_c,
                t_air_cond_in_c=self._t_air_cond_in_c,
            )
            self.optimization_converged = True
            cop_full = float(self._last_result.diagnostics.get("COP_full", 0.0))
            return cop_full, True
        except Exception:
            self.optimization_converged = False
            self._last_result = None
            return float("nan"), False

    def get_full_load_cop(self):
        """Return full-load COP from latest solve."""
        if self._last_result is None:
            return float("nan")
        return float(self._last_result.diagnostics.get("COP_full", float("nan")))

    def get_part_load_cop(self):
        """Return part-load COP from latest solve."""
        if self._last_result is None:
            return float("nan")
        return float(self._last_result.cop)
