#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Emergent-G inversion report.

This script inverts the closed mass formula using the laboratory electron mass.
It does not use the experimental value of G in the inversion.  CODATA G is used
only in the final comparison table.
"""

from __future__ import annotations

from math import sqrt

from relator_electron.common import CONSTANTS, constants_table, format_float, make_console, ppm, print_header, rich_table
from relator_electron.pipeline import run_baseline


def infer_G(C_path: float, m_ref_kg: float) -> float:
    """Infer G from m = sqrt(hbar c/G) C_path at fixed dimensionless C_path."""
    return CONSTANTS.ħ * CONSTANTS.c * (C_path / m_ref_kg) ** 2


def mass_from_G(C_path: float, G_value: float) -> float:
    """Return MeV mass from a dimensionless path coefficient and supplied G."""
    mP = sqrt(CONSTANTS.ħ * CONSTANTS.c / G_value)
    mkg = mP * C_path
    return mkg * CONSTANTS.c**2 / CONSTANTS.e / 1.0e6


def main() -> None:
    console = make_console()
    print_header(console, "Emergent G inversion", "Uses updated DC_new and finite Ward path coefficients. Author: Mehrdad Pajuhaan")
    console.print(constants_table(CONSTANTS))

    run = run_baseline()
    C_A = run.rlvm.C_A
    C_B = run.rltm.C_B
    m_ref = CONSTANTS.m_e_ref_kg

    G_A = infer_G(C_A, m_ref)
    G_B = infer_G(C_B, m_ref)
    G_geo = sqrt(G_A * G_B)

    rows = []
    for name, G_value in (("G_A from Path A", G_A), ("G_B from Path B", G_B), ("G_geo", G_geo)):
        rows.append((
            name,
            format_float(G_value, 16),
            format_float(G_value - CONSTANTS.G, 14),
            format_float(ppm((G_value - CONSTANTS.G) / CONSTANTS.G), 14),
        ))
    console.print(rich_table("Inferred Newton coupling", ("Quantity", "Value [m³ kg⁻¹ s⁻²]", "Δ vs CODATA", "Δ [ppm]"), rows))

    back_rows = []
    for name, C_path in (("Path A with G_geo", C_A), ("Path B with G_geo", C_B)):
        m = mass_from_G(C_path, G_geo)
        Δ = m - CONSTANTS.m_e_ref_MeV
        back_rows.append((name, format_float(m, 16), format_float(Δ * 1.0e6, 14), format_float(ppm(Δ / CONSTANTS.m_e_ref_MeV), 14)))
    console.print(rich_table("Back-predicted masses using G_geo", ("Case", "m c² [MeV]", "Δ [eV]", "Δ [ppm]"), back_rows))

    δAB = 2.0 * abs(G_A - G_B) / (G_A + G_B)
    console.print(rich_table("Internal G consistency", ("Diagnostic", "Value"), (("2|G_A-G_B|/(G_A+G_B)", f"{format_float(δAB * 1e9, 12)} ppb"),)))


if __name__ == "__main__":
    main()
