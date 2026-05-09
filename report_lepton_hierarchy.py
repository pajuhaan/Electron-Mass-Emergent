#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

RLVM electron anchor lifted through the imported charged-lepton hierarchy.

The ladder is dimensionless and is applied only after the Path A electron anchor
has been formed.  It does not modify the RLVM electron kernel, DC_new scalar
count, Planck prefactor, or sigma_A.
"""

from __future__ import annotations

from functools import lru_cache
from math import atan, atanh, exp, log, pi, sin, sinh, sqrt
from scipy.integrate import quad

from relator_electron.alpha_sector import DC_new, Kov
from relator_electron.common import CONSTANTS, constants_table, format_float, make_console, ppm, print_header, rich_table
from relator_electron.pipeline import run_baseline

π = pi
C0_GAUSS = 0.5 * (log(2.0) + 0.577215664901532860606512090082402431)
I_TOT = 1.0 / 6.0 - 1.0 / (4.0 * π**2)


def η_n(n: int) -> float:
    return 1.0 / (n * π)


def ℓ_n(n: int) -> float:
    return 1.0 / (n * π * sqrt(π))


def ΔΛ_out(η: float) -> float:
    return -π * (log(1.0 - η**4) / (2.0 * η) + atanh(η) - atan(η))


@lru_cache(maxsize=32)
def P_IR(ℓ: float) -> float:
    def integrand(u: float) -> float:
        weight = u**2 * sin(π * u) ** 2
        gate = 1.0 - (1.0 / 3.0) * (1.0 - u) ** 2 / ((1.0 - u) ** 2 + ℓ**2)
        return weight * gate * exp(-((1.0 - u) / ℓ) ** 2)

    value, _ = quad(integrand, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13, limit=300)
    return value / I_TOT


CORE_LOG_CONVENTION = "absolute_inverse_eta"


def L_core(n: int) -> float:
    """Imported hierarchy core log.

    The numerical charged-lepton ladder uses the absolute shell inverse-radius
    term ln(1/eta_n).  This is the convention that reproduces the imported
    hierarchy table.  It differs from an electron-anchored ln(eta_1/eta_n)
    core and should be kept explicit in the manuscript.
    """
    return log(1.0 / η_n(n)) + log(P_IR(ℓ_n(1)) / P_IR(ℓ_n(n))) + log(abs(ΔΛ_out(η_n(1))) / abs(ΔΛ_out(η_n(n))))


def k_n(n: int) -> int:
    return n * (2 * n - 1)


def w_n(n: int) -> float:
    if n == 1:
        return 1.0
    if n == 2:
        return 11.0 / 2.0
    if n == 3:
        return 41.0 / 16.0
    raise ValueError("only n=1,2,3 are used")


def gamma_geom(n: int) -> float:
    η = η_n(n)
    return 0.5 * sinh(η) / η


def kappa_curv(n: int) -> float:
    η = η_n(n)
    return sinh(η) / η - 1.0


def X_n(n: int, DC: float) -> float:
    return Kov / (2.0 * DC) * C0_GAUSS * P_IR(ℓ_n(n))


def gamma_lad(n: int, DC: float) -> tuple[float, float]:
    X = X_n(n, DC)
    Γgeom = gamma_geom(n)
    Γmap = w_n(n) * k_n(n) * X / (1.0 + kappa_curv(n) * w_n(n) * k_n(n) * X)
    γ = Γmap / (Γgeom + Γmap)
    N = (Γgeom + Γmap) * P_IR(ℓ_n(n))
    return γ, N


def ladder_logs(alpha: float) -> dict[str, float]:
    DC = DC_new(alpha)
    γ1, N1 = gamma_lad(1, DC)
    γ2, N2 = gamma_lad(2, DC)
    γ3, N3 = gamma_lad(3, DC)

    L1 = 0.0
    L2 = L_core(2) - γ2 * log(N1 / N2)
    L3 = L_core(3) - γ3 * log(N1 / N3)
    ΔL1 = DC * (0.0 * log(2.0) - log(1.0))
    ΔL2 = DC * (1.0 * log(2.0) - log(2.0))
    ΔL3 = DC * (2.0 * log(2.0) - log(3.0))
    return {
        "DC_new": DC,
        "L_e_cp1": L1,
        "L_mu_cp1": L2,
        "L_tau_cp1": L3,
        "L_e_dc": L1 + ΔL1,
        "L_mu_dc": L2 + ΔL2,
        "L_tau_dc": L3 + ΔL3,
        "DeltaL_tau_DC": ΔL3,
        "N1": N1,
        "N2": N2,
        "N3": N3,
        "gamma2": γ2,
        "gamma3": γ3,
    }


def main() -> None:
    console = make_console()
    print_header(console, "RLVM charged-lepton hierarchy report", "Dimensionless ladder applied to the updated Path A anchor. Author: Mehrdad Pajuhaan")
    console.print(constants_table(CONSTANTS))

    run = run_baseline()
    logs = ladder_logs(CONSTANTS.α)
    m_e_anchor = run.rlvm.m_A_MeV

    ladder_rows = [
        ("m_e anchor", format_float(m_e_anchor, 16), "MeV", "updated RLVM Path A"),
        ("core-log convention", CORE_LOG_CONVENTION, "—", "imported hierarchy numerical convention"),
        ("DC_new", format_float(logs["DC_new"], 16), "dimensionless", "scalar mother branch"),
        ("Delta L_tau^DC", format_float(logs["DeltaL_tau_DC"], 16), "dimensionless", "DC_new ln(4/3)"),
        ("N1,N2,N3", f"{format_float(logs['N1'], 12)}, {format_float(logs['N2'], 12)}, {format_float(logs['N3'], 12)}", "dimensionless", "sync normalizers"),
        ("gamma2,gamma3", f"{format_float(logs['gamma2'], 12)}, {format_float(logs['gamma3'], 12)}", "dimensionless", "TT map shares"),
    ]
    console.print(rich_table("Imported ladder checks", ("Quantity", "Value", "Unit", "Role"), ladder_rows))

    refs = {"e": CONSTANTS.m_e_ref_MeV, "mu": CONSTANTS.m_μ_ref_MeV, "tau": CONSTANTS.m_τ_ref_MeV}
    labels = [
        ("e", "L_e_cp1", "L_e_dc"),
        ("mu", "L_mu_cp1", "L_mu_dc"),
        ("tau", "L_tau_cp1", "L_tau_dc"),
    ]

    for variant, suffix in (("CP-1 baseline", "cp1"), ("Master log plus DC", "dc")):
        rows = []
        for lep, key_cp1, key_dc in labels:
            key = key_cp1 if suffix == "cp1" else key_dc
            L = logs[key]
            m = m_e_anchor * exp(L)
            Δ = m - refs[lep]
            rows.append((lep, format_float(L, 16), format_float(m, 16), format_float(Δ, 14), format_float(ppm(Δ / refs[lep]), 14)))
        console.print(rich_table(variant, ("Lepton", "L_n", "m c² [MeV]", "Δ [MeV]", "Δ [ppm]"), rows))


if __name__ == "__main__":
    main()
