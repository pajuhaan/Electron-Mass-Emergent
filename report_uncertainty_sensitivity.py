#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Uncertainty and sensitivity report for the Relator electron-mass code.

The report propagates physical-constant uncertainties through the closed RLVM
and RLTM masses.  For the alpha sensitivity, the shared locked geometry is held
at the reported operating point and the alpha-dependent scalar and vector blocks
are recomputed.  This matches the local sensitivity convention used in the
manuscript and avoids re-solving the shape equation for an infinitesimal
constant-propagation step.
"""

from __future__ import annotations

from math import log

from relator_electron.alpha_sector import build_alpha_blocks
from relator_electron.common import CONSTANTS, constants_table, eV_from_MeV, format_float, make_console, ppm, print_header, rich_table
from relator_electron.pipeline import run_baseline
from relator_electron.rlvm_core import compute_rlvm
from relator_electron.rltm_core import compute_rltm


def masses_at_alpha_fixed_geometry(alpha: float, base_run) -> tuple[float, float]:
    """Return Path A and Path B masses in MeV at fixed shared geometry."""
    constants = CONSTANTS.with_updates(α=alpha)
    alpha_blocks = build_alpha_blocks(alpha)
    rlvm = compute_rlvm(base_run.geometry, alpha_blocks, constants)
    rltm = compute_rltm(base_run.geometry, alpha_blocks, constants)
    return rlvm.m_A_MeV, rltm.m_B_MeV


def central_log_sensitivity_alpha(base_run, rel_step: float = 1.0e-6) -> tuple[float, float]:
    """Return local d ln m_A / d ln alpha and d ln m_B / d ln alpha."""
    α0 = CONSTANTS.α
    mA_plus, mB_plus = masses_at_alpha_fixed_geometry(α0 * (1.0 + rel_step), base_run)
    mA_minus, mB_minus = masses_at_alpha_fixed_geometry(α0 * (1.0 - rel_step), base_run)
    sA = (log(mA_plus) - log(mA_minus)) / (2.0 * rel_step)
    sB = (log(mB_plus) - log(mB_minus)) / (2.0 * rel_step)
    return sA, sB


def main() -> None:
    console = make_console()
    print_header(console, "Uncertainty and sensitivity report", "Current DC(α) and finite Ward baseline. Author: Mehrdad Pajuhaan")
    console.print(constants_table(CONSTANTS))

    base = run_baseline()
    mA = base.rlvm.m_A_MeV
    mB = base.rltm.m_B_MeV
    sA_alpha, sB_alpha = central_log_sensitivity_alpha(base)

    constants_info = [
        ("ħ", CONSTANTS.u_rel_ħ, 0.0, 0.0, "SI exact"),
        ("c", CONSTANTS.u_rel_c, 0.0, 0.0, "SI exact"),
        ("e", CONSTANTS.u_rel_e, 0.0, 0.0, "SI exact"),
        ("G", CONSTANTS.u_rel_G, -0.5, -0.5, "analytic m ∝ G^(-1/2)"),
        ("α", CONSTANTS.u_rel_α, sA_alpha, sB_alpha, "local finite difference"),
    ]

    rows_A = []
    rows_B = []
    quad_A = 0.0
    quad_B = 0.0
    for symbol, u_rel, sA, sB, note in constants_info:
        dA_MeV = abs(sA) * u_rel * mA
        dB_MeV = abs(sB) * u_rel * mB
        dA_eV = eV_from_MeV(dA_MeV)
        dB_eV = eV_from_MeV(dB_MeV)
        quad_A += dA_eV**2
        quad_B += dB_eV**2
        rows_A.append((symbol, format_float(u_rel, 12), format_float(sA, 12), format_float(dA_eV, 12), format_float(ppm(dA_MeV / mA) if mA else 0.0, 12), note))
        rows_B.append((symbol, format_float(u_rel, 12), format_float(sB, 12), format_float(dB_eV, 12), format_float(ppm(dB_MeV / mB) if mB else 0.0, 12), note))

    console.print(rich_table("Path A one-sigma uncertainty", ("Const", "rel σ", "s_A", "δA [eV]", "δA [ppm]", "Method"), rows_A))
    console.print(rich_table("Path B one-sigma uncertainty", ("Const", "rel σ", "s_B", "δB [eV]", "δB [ppm]", "Method"), rows_B))

    total_rows = [
        ("Path A", format_float(quad_A ** 0.5, 14), format_float(ppm((quad_A ** 0.5) / eV_from_MeV(mA)), 14)),
        ("Path B", format_float(quad_B ** 0.5, 14), format_float(ppm((quad_B ** 0.5) / eV_from_MeV(mB)), 14)),
    ]
    console.print(rich_table("Combined one-sigma uncertainty", ("Path", "RSS [eV]", "RSS [ppm]"), total_rows))

    ρ = base.geometry.ρstar
    DC = base.alpha_blocks.DC
    dlnNeff_dDC = ρ**2 / (1.0 + ρ**2 * DC)
    local_rows = [
        ("Path A", "d ln m_A / d ln K_EM_eff,A", "-1/4", "aggregate kernel"),
        ("Path A", "d ln m_A / d ln N_eff_eff,A", "+1/4", "aggregate scalar count"),
        ("Path A", "d ln m_A / d DC(α)", format_float(0.25 * dlnNeff_dDC, 14), "fixed N_raw and rho*"),
        ("Path A", "d ln m_A / d sigma_A", format_float(-1.0 / (4.0 * CONSTANTS.α), 14), "explicit exponential"),
        ("Path B", "d ln m_B / d ln K_T_eff,B", "+1/2", "aggregate tension kernel"),
        ("Path B", "d ln m_B / d ln K_EM_eff,B", "-1/2", "aggregate EM kernel"),
        ("Path B", "d ln m_B / d sigma_B", format_float(-1.0 / (2.0 * CONSTANTS.α), 14), "explicit exponential"),
    ]
    console.print(rich_table("Local logarithmic sensitivities", ("Path", "Derivative", "Value", "Convention"), local_rows))


if __name__ == "__main__":
    main()
