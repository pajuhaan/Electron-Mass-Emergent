#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Mehrdad Pajuhaan

Updated alpha-sector blocks used by the electron-mass calculations.

This module implements only the scalar Coulomb branch and the vector-shell
representative that are used by the corrected electron-mass paths.

Path A, RLVM, uses the scalar Coulomb branch DC_new(alpha) in the scalar count
N_eff^eff,A.  Path B, RLTM, uses the vector overlap zeta_B in the finite Ward
factors.  No combined scalar-vector bridge is computed in this codebase.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import log, pi, sqrt

from scipy.optimize import brentq


π = pi

# ---------------------------------------------------------------------------
# Current reduced vector-shell representative from the updated alpha paper.
# ---------------------------------------------------------------------------
# Completed shell-chain representative used only in the Path B finite Ward
# vector overlap.  It is dimensionless.
Λgeom_th = 0.691_683_146_106_991_026_716_159_707_634_574_473_9664

# Central ALP coefficient, denoted K_ov or LockCoeff in the paper.
Kov = (150.0 * π**2 - 8.0 * π**4 - 315.0) / (180.0 * π**6)

# Vector overlap strength used in Path B finite Ward factors.
ζB_default = Kov * Λgeom_th / (2.0 * π**2)

# Path A local inductive logarithm for the ring-local soft correction only.
# This is not the completed vector shell chain.
Λind = log(8.0 * sqrt(π)) - 2.0
ζsoft_default = Kov * Λind / (2.0 * π**2)

# ---------------------------------------------------------------------------
# Rank-5 scalar mother-law realization.
# ---------------------------------------------------------------------------
# The scalar branch solves
#       D = (alpha/pi) sqrt(R_moth(D)),
#       R_moth(D) = 1 - Theta_1_eff D + D^2 Phi_dyn(D).
# The constants below are dimensionless values of the current finite-rank
# scalar evaluator used by the alpha-sector update.
Θ1_eff = 2.746_368_406_272_133_363_756_18
s_uv = log(2.0)
s_ir = 1.0 / (8.0 * π**2)
A_uv_5 = 2.863_721_013_654_457_602_270_00
A_ir_5 = 2.850_875_150_489_833_584_190_00
ρ_dyn_5 = 0.000_440_300_042_572_573_590_000_000
χ_5 = ρ_dyn_5 * sqrt(s_uv * s_ir)

# Reference values reported by the updated alpha-sector calculation.  They are
# included as provenance checks and are not used as fitted inputs.
Dstar_th = 0.002_315_457_831_961_859_388_055_180_838_225_655_037
alpha_star_th = 0.007_297_352_565_050_600_333_38


@dataclass(frozen=True)
class AlphaSectorBlocks:
    """Dimensionless alpha-sector blocks used by Path A and Path B."""

    α: float
    DC_new: float
    ζB: float
    ζsoft: float
    Λgeom: float
    Kov: float
    R_moth: float
    Φ_dyn: float


def Phi_dyn(D: float) -> float:
    """Visible dynamic scalar response Phi_dyn(D) for the current rank-5 branch."""
    numerator = (
        A_uv_5 * (1.0 + s_ir * D)
        + A_ir_5 * (1.0 + s_uv * D)
        - 2.0 * χ_5 * D * sqrt(A_uv_5 * A_ir_5)
    )
    denominator = (1.0 + s_uv * D) * (1.0 + s_ir * D) - (χ_5**2) * (D**2)
    return numerator / denominator


def R_moth(D: float) -> float:
    """Scalar mother radicand on the admissible positive branch."""
    return 1.0 - Θ1_eff * D + D**2 * Phi_dyn(D)


@lru_cache(maxsize=4096)
def DC_new(alpha: float) -> float:
    """Return the updated scalar Coulomb branch DC_new(alpha).

    The equation solved is
        D - (alpha/pi) sqrt(R_moth(D)) = 0.
    The positive root near alpha/pi is selected.  This is the scalar branch used
    in the Path A scalar count.
    """
    if alpha <= 0.0:
        raise ValueError("alpha must be positive")
    x = alpha / π

    def residual(D: float) -> float:
        radicand = R_moth(D)
        if radicand <= 0.0:
            return D
        return D - x * sqrt(radicand)

    lo = 0.0
    hi = max(4.0 * x, 1.0e-2)
    while residual(hi) <= 0.0:
        hi *= 2.0
        if hi > 1.0:
            raise RuntimeError("failed to bracket DC_new root")
    return brentq(residual, lo, hi, xtol=1.0e-17, rtol=1.0e-15, maxiter=200)


def build_alpha_blocks(alpha: float, *, ζB: float = ζB_default, ζsoft: float = ζsoft_default) -> AlphaSectorBlocks:
    """Build the dimensionless alpha-sector block bundle at the supplied alpha."""
    D = DC_new(alpha)
    return AlphaSectorBlocks(
        α=alpha,
        DC_new=D,
        ζB=ζB,
        ζsoft=ζsoft,
        Λgeom=Λgeom_th,
        Kov=Kov,
        R_moth=R_moth(D),
        Φ_dyn=Phi_dyn(D),
    )
