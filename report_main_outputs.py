#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Main numerical report for the updated Relator electron-mass calculation.

Constants used by this report are imported from relator_electron.common.CONSTANTS
and printed at the top of the output.  The calculation uses the updated scalar
mother-law DC_new in Path A and finite Ward powers in Path B.
"""

from __future__ import annotations

from relator_electron.alpha_sector import Λind
from relator_electron.common import CONSTANTS, constants_table, eV_from_MeV, format_float, make_console, ppb, ppm, print_header, rich_table
from relator_electron.pipeline import run_baseline, solve_alpha_crosscheck


def main() -> None:
    console = make_console()
    run = run_baseline()
    c = run.constants
    g = run.geometry
    a = run.alpha_blocks
    A = run.rlvm
    B = run.rltm

    print_header(
        console,
        "Updated Relator electron-mass run",
        "Path A uses DC_new; Path B uses finite Ward powers. Author: Mehrdad Pajuhaan",
    )
    console.print(constants_table(c))

    shared_rows = [
        ("α", format_float(c.α, 16), "dimensionless", "run input"),
        ("y*", format_float(g.ystar, 16), "dimensionless", "shape-stationarity root"),
        ("K_ring(y*)", format_float(g.Kring, 16), "dimensionless", "near-field K0"),
        ("L_ring(y*)", format_float(g.Lring, 16), "dimensionless", "near-field K1"),
        ("G_ring(y*)", format_float(g.Gring, 16), "dimensionless", "K0 - K1, not used in ZigmaRing"),
        ("ZigmaRing(y*)", format_float(g.ZigmaRing, 16), "dimensionless", "K_ring L_ring /(2π²)"),
        ("Fstruct(y*)", format_float(g.Fstruct, 16), "dimensionless", "2 + ZigmaRing + α K_ring/(8π²)"),
        ("D", format_float(g.D, 16), "dimensionless", "4 Fstruct/3"),
        ("rho*", format_float(g.ρstar, 16), "dimensionless", "1/D"),
        ("x", format_float(g.x, 16), "dimensionless", "rho*/y*"),
        ("kappa_T", format_float(g.κTens, 16), "dimensionless", "x²/(1+x²)"),
        ("beta_rel", format_float(g.βrel, 16), "dimensionless", "internal Relator slope"),
        ("sigma_A", format_float(g.σA, 16), "dimensionless", "Path A exponent"),
        ("sigma_B", format_float(g.σB, 16), "dimensionless", "Path B exponent"),
    ]
    console.print(rich_table("Shared locked geometry", ("Quantity", "Value", "Unit", "Definition"), shared_rows))

    alpha_rows = [
        ("Kov", format_float(a.Kov, 16), "dimensionless", "ALP angular overlap coefficient"),
        ("Lambda_geom^(th)", format_float(a.Λgeom, 16), "dimensionless", "current reduced vector representative"),
        ("zeta_B", format_float(a.ζB, 16), "dimensionless", "Kov Lambda_geom /(2π²), Path B only"),
        ("Lambda_ind", format_float(Λind, 16), "dimensionless", "ring-local inductive logarithm"),
        ("zeta_soft", format_float(a.ζsoft, 16), "dimensionless", "Kov Lambda_ind /(2π²), Path A local EM only"),
        ("DC_new(alpha)", format_float(a.DC_new, 16), "dimensionless", "updated scalar mother-law branch"),
        ("R_moth(DC_new)", format_float(a.R_moth, 16), "dimensionless", "scalar mother radicand"),
        ("Phi_dyn(DC_new)", format_float(a.Φ_dyn, 16), "dimensionless", "rank-5 visible dynamic response"),
    ]
    console.print(rich_table("Updated alpha-sector blocks", ("Quantity", "Value", "Unit", "Role"), alpha_rows))

    pathA_rows = [
        ("K_EM_raw", format_float(A.K_EM_raw, 16), "dimensionless", "raw along-ring elliptic kernel"),
        ("s_G", format_float(A.sG, 16), "dimensionless", "Gaussian Maxwell softening"),
        ("K_EM_mid,A", format_float(A.K_EM_mid_A, 16), "dimensionless", "K_EM_raw s_G"),
        ("K_EM_eff,A", format_float(A.K_EM_eff_A, 16), "dimensionless", "mid × (1 - 1/2 rho*² zeta_soft)"),
        ("N_eff_raw,A", format_float(A.N_eff_raw_A, 16), "dimensionless", "2(1+kappa)(1-kappa²/2)"),
        ("N_eff_eff,A", format_float(A.N_eff_eff_A, 16), "dimensionless", "raw × (1 + rho*² DC_new)"),
        ("a_A", format_float(A.a_A, 16), "J m", "1/R coefficient"),
        ("b_A", format_float(A.b_A, 16), "J m⁻³", "R³ coefficient before exponential"),
        ("R_A", format_float(A.R_A, 16), "m", "stationary closure radius"),
        ("rlock_A", format_float(A.rlock_A, 16), "m", "R_A/D"),
        ("m_A c²", format_float(A.m_A_MeV, 16), "MeV", "baseline RLVM"),
    ]
    console.print(rich_table("Path A, RLVM core outputs", ("Quantity", "Value", "Unit", "Definition"), pathA_rows))

    pathB_rows = [
        ("K_EM_raw", format_float(B.K_EM_raw, 16), "dimensionless", "raw along-ring elliptic kernel"),
        ("K_T_raw,B", format_float(B.K_T_raw_B, 16), "dimensionless", "2 kappa_T"),
        ("delta_loc", format_float(B.delta_loc, 16), "dimensionless", "rho*² ZigmaRing"),
        ("delta_bulk", format_float(B.delta_bulk, 16), "dimensionless", "rho*² zeta_B"),
        ("Ward T loc", format_float(B.ward_T_loc, 16), "dimensionless", "(1+delta_loc)^(-1/2)"),
        ("Ward T bulk", format_float(B.ward_T_bulk, 16), "dimensionless", "(1+delta_bulk)^(-1/2)"),
        ("Ward EM loc", format_float(B.ward_EM_loc, 16), "dimensionless", "(1+delta_loc)^(1/4)"),
        ("Ward EM bulk", format_float(B.ward_EM_bulk, 16), "dimensionless", "(1+delta_bulk)^(1/4)"),
        ("K_T_eff,B", format_float(B.K_T_eff_B, 16), "dimensionless", "finite Ward tension kernel"),
        ("K_EM_eff,B", format_float(B.K_EM_eff_B, 16), "dimensionless", "finite Ward EM kernel"),
        ("f_M1", format_float(B.fM1, 16), "dimensionless", "mode-selective tangential remainder"),
        ("G_ten", format_float(B.Gten, 16), "dimensionless", "exp(-sigma_B/alpha)"),
        ("a_B", format_float(B.a_B, 16), "J m", "1/R coefficient"),
        ("T", format_float(B.T, 16), "J m⁻¹", "string tension"),
        ("R_B", format_float(B.R_B, 16), "m", "stationary closure radius"),
        ("rlock_B", format_float(B.rlock_B, 16), "m", "R_B/D"),
        ("m_B c²", format_float(B.m_B_MeV, 16), "MeV", "baseline RLTM"),
    ]
    console.print(rich_table("Path B, RLTM core outputs", ("Quantity", "Value", "Unit", "Definition"), pathB_rows))

    compare_rows = []
    for name, m in (("Path A", A.m_A_MeV), ("Path B", B.m_B_MeV)):
        Δ = m - c.m_e_ref_MeV
        compare_rows.append((name, format_float(m, 16), format_float(eV_from_MeV(Δ), 16), format_float(ppm(Δ / c.m_e_ref_MeV), 16)))
    compare_rows.append(("A - B", "", format_float(eV_from_MeV(A.m_A_MeV - B.m_B_MeV), 16), ""))
    console.print(rich_table("Mass comparison", ("Case", "m c² [MeV]", "Δ [eV]", "Δ [ppm]"), compare_rows))

    αhat = solve_alpha_crosscheck(run)
    alpha_rows = [
        ("alpha input", format_float(c.α, 18), "baseline"),
        ("alpha_hat", format_float(αhat, 18), f"relative {format_float(ppb((αhat-c.α)/c.α), 12)} ppb"),
    ]
    console.print(rich_table("Path equality alpha diagnostic", ("Quantity", "Value", "Comment"), alpha_rows))


if __name__ == "__main__":
    main()
