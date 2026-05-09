#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Path A core, Relator-Locked Vacuum-in-Volume Model, RLVM.

Baseline physics implemented here
---------------------------------
1. The Path A scalar channel uses DC_new(alpha), the updated scalar mother-law
   Coulombic branch from the alpha-sector update.
2. The Path A local electromagnetic kernel contains only the ring-local induced
   scalar zeta_soft built from Lambda_ind.
3. The completed vector-shell factor zeta_B is not inserted into K_EM^eff,A.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, pi

from .alpha_sector import AlphaSectorBlocks
from .common import PhysicalConstants, CONSTANTS
from .geometry import GeometryBlocks, K_EM_raw, N_eff_raw_A, sG_softening


π = pi


@dataclass(frozen=True)
class RLVMResult:
    """Complete Path A output at fixed constants and shared geometry."""

    K_EM_raw: float
    sG: float
    K_EM_mid_A: float
    K_EM_eff_A: float
    N_eff_raw_A: float
    N_eff_eff_A: float
    a_A: float
    b_A: float
    R_A: float
    rlock_A: float
    m_A_kg: float
    m_A_MeV: float
    C_A: float


def compute_rlvm(
    geometry: GeometryBlocks,
    alpha_blocks: AlphaSectorBlocks,
    constants: PhysicalConstants = CONSTANTS,
) -> RLVMResult:
    """Compute the RLVM path with the updated scalar DC_new branch."""
    α = constants.α
    ρ = geometry.ρstar
    x = geometry.x
    κ = geometry.κTens

    Kraw = K_EM_raw(x)
    sG = sG_softening(κ)
    Kmid = Kraw * sG
    K_eff_A = Kmid * (1.0 - 0.5 * ρ**2 * alpha_blocks.ζsoft)

    Nraw = N_eff_raw_A(x)
    Neff = Nraw * (1.0 + ρ**2 * alpha_blocks.DC_new)

    ζC = π
    σA = geometry.σA
    lock = exp(-σA / α)

    a_A = (4.0 * π * α * constants.ħ * constants.c / ρ) * K_eff_A
    b_A = (constants.ħ * constants.c / constants.ℓP**4) * (Neff * ζC**4 * ρ**2 / 8.0)

    R_A = constants.ℓP * exp(σA / (4.0 * α)) * (
        (32.0 * π * α * K_eff_A) / (3.0 * Neff * ζC**4 * ρ**3)
    ) ** 0.25

    m_A_kg = geometry.D * constants.ħ / (constants.c * R_A)
    m_A_MeV = m_A_kg * constants.c**2 / constants.e / 1.0e6
    rlock_A = R_A / geometry.D
    C_A = m_A_kg / constants.mP

    # Algebraic cross-check against the direct mass formula.
    C_A_direct = geometry.D * (
        3.0 * Neff * ζC**4 * ρ**3 * lock / (32.0 * π * α * K_eff_A)
    ) ** 0.25
    if abs(C_A - C_A_direct) > 5.0e-12 * max(1.0, abs(C_A)):
        raise RuntimeError("RLVM algebraic consistency check failed")

    return RLVMResult(
        K_EM_raw=Kraw,
        sG=sG,
        K_EM_mid_A=Kmid,
        K_EM_eff_A=K_eff_A,
        N_eff_raw_A=Nraw,
        N_eff_eff_A=Neff,
        a_A=a_A,
        b_A=b_A,
        R_A=R_A,
        rlock_A=rlock_A,
        m_A_kg=m_A_kg,
        m_A_MeV=m_A_MeV,
        C_A=C_A,
    )
