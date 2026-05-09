#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Pipeline assembly for the Relator electron-mass calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import exp, pi

from scipy.optimize import brentq

from .alpha_sector import build_alpha_blocks, DC
from .common import PhysicalConstants, CONSTANTS
from .geometry import build_geometry
from .rlvm_core import RLVMResult, compute_rlvm
from .rltm_core import RLTMResult, compute_rltm, finite_ward_factors


@dataclass(frozen=True)
class FullRun:
    """Bundle returned by the complete electron-mass run."""

    constants: PhysicalConstants
    geometry: object
    alpha_blocks: object
    rlvm: RLVMResult
    rltm: RLTMResult


@lru_cache(maxsize=256)
def run_baseline(alpha: float = CONSTANTS.α, G: float = CONSTANTS.G) -> FullRun:
    """Run the complete baseline calculation at the supplied alpha and G."""
    constants = CONSTANTS.with_updates(α=alpha, G=G)
    geometry = build_geometry(constants.α)
    alpha_blocks = build_alpha_blocks(constants.α)
    rlvm = compute_rlvm(geometry, alpha_blocks, constants)
    rltm = compute_rltm(geometry, alpha_blocks, constants)
    return FullRun(constants=constants, geometry=geometry, alpha_blocks=alpha_blocks, rlvm=rlvm, rltm=rltm)


def alpha_fixed_point_map(alpha_unknown: float, run: FullRun) -> float:
    """Path A / Path B equality map at the fixed shared locked geometry.

    The scalar entry is always the current Relator DC(alpha).  The vector-shell
    representative ζ_B is a fixed Path B Ward input of the selected run.
    """
    g = run.geometry
    A = run.rlvm
    B = run.rltm
    ab = run.alpha_blocks
    ρ = g.ρstar

    D_scalar = DC(alpha_unknown)
    Neff = A.N_eff_raw_A * (1.0 + ρ**2 * D_scalar)
    fM1 = 0.5 * alpha_unknown * g.ystar**2 * exp(-0.5 * g.ystar**2)

    # ζ_B is the vector-shell representative selected by the run.  It enters
    # only through Path B finite Ward factors.
    δloc = B.delta_loc
    δbulk = ρ**2 * ab.ζB
    wTloc, wEMloc = finite_ward_factors(δloc)
    wTbulk, wEMbulk = finite_ward_factors(δbulk)
    KTeff = B.K_T_raw_B * wTloc * wTbulk
    KEMeff = B.K_EM_raw * wEMloc * wEMbulk

    return (
        8.0
        / (3.0 * pi**3)
        * A.K_EM_eff_A
        / (Neff * ρ)
        * (KTeff / KEMeff) ** 2
        * (1.0 + fM1) ** 2
    )


def solve_alpha_crosscheck(run: FullRun) -> float:
    """Solve alpha = A_AB(alpha) for the Path A / Path B equality diagnostic."""

    def residual(alpha_value: float) -> float:
        return alpha_value - alpha_fixed_point_map(alpha_value, run)

    lo = 0.006
    hi = 0.0085
    if residual(lo) * residual(hi) > 0.0:
        lo = 0.001
        hi = 0.02
    return brentq(residual, lo, hi, xtol=2.0e-16, rtol=1.0e-14, maxiter=200)
