#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Path B core, Relator-Locked Tension Model, RLTM.

Baseline physics implemented here
---------------------------------
1. The Path B vector channel uses finite Ward powers, not the linearized
   O(rho_star^2) truncation.
2. Near-edge geometry enters through delta_loc = rho_star^2 ZigmaRing.
3. Bulk vector completion enters through delta_bulk = rho_star^2 ζ_B.
4. No Path A Gaussian softening factor is inserted into Path B.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, pi

from .alpha_sector import AlphaSectorBlocks
from .common import PhysicalConstants, CONSTANTS
from .geometry import GeometryBlocks, K_EM_raw, K_T_raw_B


π = pi


@dataclass(frozen=True)
class RLTMResult:
    """Complete Path B output at fixed constants and shared geometry."""

    K_EM_raw: float
    K_T_raw_B: float
    delta_loc: float
    delta_bulk: float
    ward_T_loc: float
    ward_T_bulk: float
    ward_EM_loc: float
    ward_EM_bulk: float
    K_T_eff_B: float
    K_EM_eff_B: float
    fM1: float
    Gten: float
    a_B: float
    T: float
    R_B: float
    rlock_B: float
    m_B_kg: float
    m_B_MeV: float
    C_B: float


def finite_ward_factors(delta: float) -> tuple[float, float]:
    """Return finite Ward factors for one scalar source.

    The exact finite factors are
        tension factor = (1+delta)^(-1/2),
        EM factor      = (1+delta)^(1/4).
    They obey tension_factor * EM_factor^2 = 1 exactly up to floating arithmetic.
    """
    if delta <= -1.0:
        raise ValueError("finite Ward source must satisfy delta > -1")
    return (1.0 + delta) ** (-0.5), (1.0 + delta) ** 0.25


def compute_rltm(
    geometry: GeometryBlocks,
    alpha_blocks: AlphaSectorBlocks,
    constants: PhysicalConstants = CONSTANTS,
) -> RLTMResult:
    """Compute the baseline RLTM path with finite Ward completion."""
    α = constants.α
    ρ = geometry.ρstar
    x = geometry.x

    Kraw = K_EM_raw(x)
    KTraw = K_T_raw_B(x)

    δloc = ρ**2 * geometry.ZigmaRing
    δbulk = ρ**2 * alpha_blocks.ζB
    ward_T_loc, ward_EM_loc = finite_ward_factors(δloc)
    ward_T_bulk, ward_EM_bulk = finite_ward_factors(δbulk)

    KTeff = KTraw * ward_T_loc * ward_T_bulk
    KEMeff = Kraw * ward_EM_loc * ward_EM_bulk

    fM1 = 0.5 * α * geometry.ystar**2 * exp(-0.5 * geometry.ystar**2)
    σB = geometry.σB
    Gten = exp(-σB / α)

    a_B = (4.0 * π * α * constants.ħ * constants.c / ρ) * KEMeff
    T = (constants.ħ * constants.c / constants.ℓP**2) * KTeff * Gten * (1.0 + fM1)

    R_B = constants.ℓP * (
        (2.0 * α * KEMeff) / (ρ * KTeff * Gten * (1.0 + fM1))
    ) ** 0.5
    m_B_kg = geometry.D * constants.ħ / (constants.c * R_B)
    m_B_MeV = m_B_kg * constants.c**2 / constants.e / 1.0e6
    rlock_B = R_B / geometry.D
    C_B = m_B_kg / constants.mP

    # Exact Ward-neutral identity check for both finite factors.
    ward_check = (ward_T_loc * ward_EM_loc**2) * (ward_T_bulk * ward_EM_bulk**2)
    if abs(ward_check - 1.0) > 2.0e-14:
        raise RuntimeError("finite Ward identity check failed")

    return RLTMResult(
        K_EM_raw=Kraw,
        K_T_raw_B=KTraw,
        delta_loc=δloc,
        delta_bulk=δbulk,
        ward_T_loc=ward_T_loc,
        ward_T_bulk=ward_T_bulk,
        ward_EM_loc=ward_EM_loc,
        ward_EM_bulk=ward_EM_bulk,
        K_T_eff_B=KTeff,
        K_EM_eff_B=KEMeff,
        fM1=fM1,
        Gten=Gten,
        a_B=a_B,
        T=T,
        R_B=R_B,
        rlock_B=rlock_B,
        m_B_kg=m_B_kg,
        m_B_MeV=m_B_MeV,
        C_B=C_B,
    )
