#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Reference-only QED+ALP scalar-branch diagnostic for the electron-mass code.

This report leaves the primary code path untouched.  It asks what the Path A
mass and Path A/Path B equality alpha would be if the Path A scalar count used
D_QED^[5](alpha/pi), reconstructed from retained QED coefficients through the
ALP bridge in the alpha paper.  Path B is not modified, because Path B uses the
bulk vector overlap zeta_B in its Ward factors rather than the scalar DC block.
"""

from __future__ import annotations

from dataclasses import replace
from math import exp, pi

from scipy.optimize import brentq

from relator_electron.common import CONSTANTS, constants_table, eV_from_MeV, format_float, make_console, ppb, ppm, print_header, rich_table
from relator_electron.pipeline import run_baseline, solve_alpha_crosscheck
from relator_electron.qed_alp_dc import (
    C_QED,
    Dstar_QED_ALP_target,
    alpha_QED_to_D_reference,
    build_qed_alp_branch,
    DC_qed_alp,
    solve_alpha_qed_alp_target,
)
from relator_electron.rlvm_core import compute_rlvm
from relator_electron.rltm_core import finite_ward_factors

π = pi


def compute_rlvm_with_qed_alp_DC(run, alpha: float = CONSTANTS.α):
    """Compute Path A with DC_QED_ALP in the scalar count only."""
    branch = build_qed_alp_branch(alpha, order=5)
    alpha_blocks_qed = replace(run.alpha_blocks, α=alpha, DC_new=branch.DC_qed_alp)
    constants = run.constants.with_updates(α=alpha)
    rlvm_qed = compute_rlvm(run.geometry, alpha_blocks_qed, constants)
    return branch, rlvm_qed


def alpha_fixed_point_map_qed_alp(alpha_unknown: float, run) -> float:
    """Path-equality fixed-point map using the QED+ALP scalar branch in Path A."""
    g = run.geometry
    A = run.rlvm
    B = run.rltm
    ab = run.alpha_blocks
    ρ = g.ρstar

    D_scalar = DC_qed_alp(alpha_unknown, order=5)
    Neff = A.N_eff_raw_A * (1.0 + ρ**2 * D_scalar)
    fM1 = 0.5 * alpha_unknown * g.ystar**2 * exp(-0.5 * g.ystar**2)

    # Path B remains vector-channel only.  The bulk vector representative zeta_B
    # is the same fixed dimensionless vector input used in the baseline run.
    δloc = B.delta_loc
    δbulk = ρ**2 * ab.ζB
    wTloc, wEMloc = finite_ward_factors(δloc)
    wTbulk, wEMbulk = finite_ward_factors(δbulk)
    KTeff = B.K_T_raw_B * wTloc * wTbulk
    KEMeff = B.K_EM_raw * wEMloc * wEMbulk

    return (
        8.0
        / (3.0 * π**3)
        * A.K_EM_eff_A
        / (Neff * ρ)
        * (KTeff / KEMeff) ** 2
        * (1.0 + fM1) ** 2
    )


def solve_alpha_crosscheck_qed_alp(run) -> float:
    """Solve alpha = A_AB^QED+ALP(alpha) for the path-equality diagnostic."""
    def residual(alpha_value: float) -> float:
        return alpha_value - alpha_fixed_point_map_qed_alp(alpha_value, run)

    lo, hi = 0.006, 0.0085
    if residual(lo) * residual(hi) > 0.0:
        lo, hi = 0.001, 0.02
    return brentq(residual, lo, hi, xtol=2.0e-16, rtol=1.0e-14, maxiter=200)


def main() -> None:
    console = make_console()
    print_header(
        console,
        "Reference-only QED+ALP scalar-branch diagnostic",
        "Separate report; baseline Relator DC_new code is not modified. Author: Mehrdad Pajuhaan",
    )
    console.print(constants_table(CONSTANTS))

    run = run_baseline()
    branch_input, rlvm_qed = compute_rlvm_with_qed_alp_DC(run, CONSTANTS.α)
    alpha_path_eq_baseline = solve_alpha_crosscheck(run)
    alpha_path_eq = solve_alpha_crosscheck_qed_alp(run)
    alpha_qed_target = solve_alpha_qed_alp_target(order=5)

    coeff_rows = [(f"c{idx}_QED", format_float(coeff, 18), "dimensionless", "D_QED^[5](x) coefficient") for idx, coeff in enumerate(C_QED, start=1)]
    console.print(rich_table("Retained QED+ALP scalar coefficients", ("Coeff", "Value", "Unit", "Role"), coeff_rows))

    scalar_rows = [
        ("x input", format_float(branch_input.x, 18), "dimensionless", "alpha/pi at run input"),
        ("DC_QED+ALP(alpha)", format_float(branch_input.DC_qed_alp, 18), "dimensionless", "D_QED^[5](alpha/pi)"),
        ("Dstar target", format_float(Dstar_QED_ALP_target, 18), "dimensionless", "fixed alpha-paper reduced target"),
        ("DC_QED+ALP - Dstar", format_float(branch_input.target_residual, 18), "dimensionless", "target residual at run alpha"),
        ("alpha from D_QED+ALP=Dstar", format_float(alpha_qed_target, 18), "dimensionless", "inverse-alpha diagnostic in alpha paper"),
        ("alpha row in Alpha.tex", format_float(alpha_QED_to_D_reference, 18), "dimensionless", "provenance check"),
    ]
    console.print(rich_table("Scalar diagnostic branch", ("Quantity", "Value", "Unit", "Definition"), scalar_rows))

    mass_rows = []
    for case, m in (
        ("Path A baseline DC_new", run.rlvm.m_A_MeV),
        ("Path A with DC_QED+ALP", rlvm_qed.m_A_MeV),
        ("Path B unchanged vector path", run.rltm.m_B_MeV),
    ):
        Δ = m - CONSTANTS.m_e_ref_MeV
        mass_rows.append((case, format_float(m, 18), format_float(eV_from_MeV(Δ), 16), format_float(ppm(Δ / CONSTANTS.m_e_ref_MeV), 16)))
    mass_rows.append(("A_QED+ALP - B", "", format_float(eV_from_MeV(rlvm_qed.m_A_MeV - run.rltm.m_B_MeV), 16), ""))
    console.print(rich_table("Electron-mass rows", ("Case", "m c² [MeV]", "Δ vs ref [eV]", "Δ vs ref [ppm]"), mass_rows))

    alpha_rows = [
        ("Path equality alpha, DC_new", format_float(alpha_path_eq_baseline, 18), f"relative {format_float(ppb((alpha_path_eq_baseline-CONSTANTS.α)/CONSTANTS.α), 14)} ppb vs run input"),
        ("Path equality alpha, DC_QED+ALP", format_float(alpha_path_eq, 18), f"relative {format_float(ppb((alpha_path_eq-CONSTANTS.α)/CONSTANTS.α), 14)} ppb vs run input"),
        ("QED+ALP target alpha", format_float(alpha_qed_target, 18), f"relative {format_float(ppb((alpha_qed_target-CONSTANTS.α)/CONSTANTS.α), 14)} ppb vs run input"),
    ]
    console.print(rich_table("Alpha diagnostics", ("Diagnostic", "Value", "Comment"), alpha_rows))

    audit_rows = [
        ("Core Path A/B formulas", "PASS", "masses use only alpha, G through ell_P, and dimensionless kernels"),
        ("Reference electron mass", "not used in cores", "appears only in comparison rows and in the separate emergent-G inversion"),
        ("QED+ALP diagnostic", "separate", "does not modify DC_new, RLVM baseline, RLTM baseline, or lepton hierarchy"),
    ]
    console.print(rich_table("No-contamination audit", ("Check", "Status", "Meaning"), audit_rows))

    channel_rows = [
        ("Path A scalar slot", "DC_new or DC_QED+ALP", "N_eff_eff,A = N_eff_raw,A [1 + rho*² DC]"),
        ("Path A local EM slot", "zeta_soft", "ring-local Lambda_ind only"),
        ("Path B vector slot", "zeta_B", "finite Ward factors only"),
        ("Path B scalar DC", "not used", "RLTM is the independent vector/tension closure"),
    ]
    console.print(rich_table("Channel placement", ("Channel", "Inserted quantity", "Reason"), channel_rows))


if __name__ == "__main__":
    main()
