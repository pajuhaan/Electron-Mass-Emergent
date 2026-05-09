#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Finite-width S1 ring geometry used by both electron-mass paths.

The root y_star is obtained from the shape stationarity equation
    3 F_struct'(y)/F_struct(y) + 2/y - y = 0.
All quantities in this module are dimensionless except where explicitly stated.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import cos, exp, pi, sin, sqrt

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erfcx
from scipy.special import ellipk


π = pi


@dataclass(frozen=True)
class GeometryBlocks:
    """Shared locked S1 geometry."""

    α: float
    ystar: float
    Kring: float
    Lring: float
    Gring: float
    ZigmaRing: float
    Fstruct: float
    D: float
    ρstar: float
    x: float
    κTens: float
    βrel: float
    σA: float
    σB: float


def chord_S1(Δ: float) -> float:
    """Dimensionless chord on S1 for angular separation Delta."""
    return 2.0 * sin(abs(Δ) / 2.0)


def J_kernel(a: float, y: float) -> float:
    """Gaussian near-field kernel J(a,y) = sqrt(pi)/(2y) erfcx(a/(2y))."""
    return sqrt(π) / (2.0 * y) * erfcx(a / (2.0 * y))


@lru_cache(maxsize=2048)
def Kring_of(y: float) -> float:
    """K0(y), the scalar near-field ring harmonic."""
    value, _ = quad(lambda Δ: J_kernel(chord_S1(Δ), y), 0.0, π, epsabs=1e-13, epsrel=1e-13, limit=300)
    return value / π


@lru_cache(maxsize=2048)
def Lring_of(y: float) -> float:
    """K1(y), the first near-field ring harmonic."""
    value, _ = quad(lambda Δ: cos(Δ) * J_kernel(chord_S1(Δ), y), 0.0, π, epsabs=1e-13, epsrel=1e-13, limit=300)
    return value / π


def Fstruct_of(y: float, alpha: float) -> float:
    """Structure factor F_struct(y)."""
    K = Kring_of(y)
    L = Lring_of(y)
    Z = K * L / (2.0 * π**2)
    return 2.0 + Z + alpha * K / (8.0 * π**2)


def dFstruct_dy(y: float, alpha: float) -> float:
    """Central finite-difference derivative of F_struct with respect to y."""
    h = 5.0e-6 * (1.0 + abs(y))
    return (Fstruct_of(y + h, alpha) - Fstruct_of(y - h, alpha)) / (2.0 * h)


def shape_stationarity(y: float, alpha: float) -> float:
    """Shape stationarity residual."""
    F = Fstruct_of(y, alpha)
    return 3.0 * dFstruct_dy(y, alpha) / F + 2.0 / y - y


def solve_ystar(alpha: float) -> float:
    """Solve the shared S1 shape-lock equation for y_star."""
    y_lo = 0.35
    step = 0.02
    f_lo = shape_stationarity(y_lo, alpha)
    y = y_lo + step
    while y <= 3.0:
        f_y = shape_stationarity(y, alpha)
        if f_lo * f_y < 0.0:
            return brentq(lambda yy: shape_stationarity(yy, alpha), y - step, y, xtol=1.0e-13, rtol=1.0e-13, maxiter=200)
        y_lo = y
        f_lo = f_y
        y += step
    raise RuntimeError("no y_star bracket found in [0.35, 3.0]")


def kappa_tension(x: float) -> float:
    """Locked second-moment fraction kappa = x^2/(1+x^2)."""
    return x * x / (1.0 + x * x)


def K_EM_raw(x: float) -> float:
    """Raw along-ring elliptic electromagnetic kernel."""
    return 2.0 * x / π * ellipk(-4.0 * x * x)


def sG_softening(κ: float) -> float:
    """Gaussian Maxwell softening factor for Path A local EM kernel."""
    if κ <= 0.0:
        return 1.0
    if κ < 1.0e-12:
        return 1.0 - κ / 2.0 + 0.75 * κ * κ
    return sqrt(π / κ) * erfcx(1.0 / sqrt(κ))


def N_eff_raw_A(x: float) -> float:
    """Raw scalar channel count used in Path A."""
    κ = kappa_tension(x)
    return 2.0 * (1.0 + κ) * (1.0 - 0.5 * κ**2)


def K_T_raw_B(x: float) -> float:
    """Raw tangential tension kernel used in Path B."""
    return 2.0 * kappa_tension(x)


def build_geometry(alpha: float) -> GeometryBlocks:
    """Construct all shared locked S1 geometry blocks."""
    ystar = solve_ystar(alpha)
    K = Kring_of(ystar)
    L = Lring_of(ystar)
    Ggap = K - L
    Z = K * L / (2.0 * π**2)
    Fstruct = 2.0 + Z + alpha * K / (8.0 * π**2)
    D = 4.0 * Fstruct / 3.0
    ρstar = 1.0 / D
    x = ρstar / ystar
    κ = kappa_tension(x)
    βrel = (6.0 / π) * (D / Fstruct)
    σA = 4.0 / βrel
    σB = 2.0 / βrel
    return GeometryBlocks(
        α=alpha,
        ystar=ystar,
        Kring=K,
        Lring=L,
        Gring=Ggap,
        ZigmaRing=Z,
        Fstruct=Fstruct,
        D=D,
        ρstar=ρstar,
        x=x,
        κTens=κ,
        βrel=βrel,
        σA=σA,
        σB=σB,
    )
