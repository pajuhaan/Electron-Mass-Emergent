#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Mehrdad Pajuhaan

Reference-only QED+ALP scalar branch for an isolated diagnostic report.

This module implements the retained M=5 QED-induced scalar branch described in
Alpha.tex.  It is deliberately separate from the primary Relator scalar branch
`DC(α)` used by the baseline electron-mass code.

The diagnostic branch is

    DCQEDALP(alpha) = sum_{n=1}^5 c_n^QED (alpha/pi)^n,

where the coefficients c_n^QED are the retained scalar coefficients obtained
from the pure-photonic QED benchmark through the ALP bridge in the alpha paper.
Only the scalar audit branch DCQEDALP is exported here.

Nothing in this module is used by the baseline Path A or Path B calculations.
The report that imports it is a separate one-way diagnostic.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

from scipy.optimize import brentq

π = pi

# Retained pure-photonic benchmark coefficients for a_e^QED(x), x = alpha/pi.
a1_QED = 0.5
a2_QED = -0.328_478_965_579_193
a3_QED = 1.181_241_456_587
a4_QED = -1.912_245_764_926_445_574_152_647_153_1265
a5_QED = 5.891

# Retained ALP-bridge scalar coefficients c_n^QED through M = 5.
c1_QED = 1.0
c2_QED = -1.373_184_203_277_606_742_852_040_178_93
c3_QED = 3.800_116_429_431_190_284_498_232_364_10
c4_QED = -8.857_883_408_219_170_032_455_881_235_19
c5_QED = 24.735_181_450_506_472_598_325_174_337_7

C_QED: tuple[float, ...] = (c1_QED, c2_QED, c3_QED, c4_QED, c5_QED)

# Current reduced target reported in the alpha paper.  It is included here only
# to reproduce the inverse-alpha diagnostic of the QED+ALP branch.
Dstar_QED_ALP_target = 0.002_315_457_831_961_859_388_055
alpha_QED_to_D_reference = 0.007_297_352_564_484_953_942_37


@dataclass(frozen=True)
class QEDALPBranch:
    """Reference-only QED+ALP scalar branch values at one alpha."""

    α: float
    x: float
    DC_qed_alp: float
    target_residual: float
    truncation_order: int


def DC_qed_alp(alpha: float, order: int = 5) -> float:
    """Evaluate the retained QED+ALP scalar branch DCQEDALP_M(alpha).

    Parameters
    ----------
    alpha
        Positive dimensionless fine-structure coupling.
    order
        Retained truncation order, 1 <= order <= 5.  The manuscript diagnostic
        uses order = 5.
    """
    if alpha <= 0.0:
        raise ValueError("alpha must be positive")
    if not 1 <= order <= len(C_QED):
        raise ValueError("order must satisfy 1 <= order <= 5")
    x = alpha / π
    value = 0.0
    power = x
    for coeff in C_QED[:order]:
        value += coeff * power
        power *= x
    return value


def build_qed_alp_branch(alpha: float, order: int = 5) -> QEDALPBranch:
    """Return the QED+ALP diagnostic scalar branch bundle."""
    D = DC_qed_alp(alpha, order=order)
    return QEDALPBranch(
        α=alpha,
        x=alpha / π,
        DC_qed_alp=D,
        target_residual=D - Dstar_QED_ALP_target,
        truncation_order=order,
    )


def solve_alpha_qed_alp_target(order: int = 5) -> float:
    """Solve DCQEDALP_M(alpha) = Dstar_QED_ALP_target.

    This reproduces the reference-only inverse-alpha diagnostic in the alpha
    paper.  It is not used in the electron-mass baseline.
    """
    def residual(alpha: float) -> float:
        return DC_qed_alp(alpha, order=order) - Dstar_QED_ALP_target

    return brentq(residual, 0.006, 0.0085, xtol=2.0e-16, rtol=1.0e-14, maxiter=200)
