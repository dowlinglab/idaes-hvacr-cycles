# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Shared data structures used by cycle and HX model modules.

Context
-------
Shared types are isolated here to avoid circular imports while keeping the
same callable interface across HX variants.
"""

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class RefrigerantState:
    """Refrigerant thermodynamic state at a stream point.

    Attributes
    ----------
    p:
        Pressure [Pa].
    h:
        Specific enthalpy [J/kg].
    m_dot:
        Mass flow rate [kg/s].
    """

    p: float
    h: float
    m_dot: float


@dataclass(frozen=True)
class AirStream:
    """Secondary-side stream representation for HX calls.

    Attributes
    ----------
    t_in:
        Inlet temperature [K].
    m_dot:
        Mass flow rate [kg/s].
    cp:
        Heat capacity [J/(kg-K)].
    """

    t_in: float
    m_dot: float
    cp: float


@dataclass(frozen=True)
class CycleResult:
    """Cycle solution summary for one operating point."""

    fluid: str
    t_evap_sat_c: float
    t_cond_sat_c: float
    p_evap_pa: float
    p_cond_pa: float
    m_dot_ref: float
    q_evap_w: float
    q_cond_w: float
    w_comp_w: float
    cop: float
    diagnostics: Dict[str, float]
