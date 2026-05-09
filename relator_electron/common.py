#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Shared constants, units, and report-formatting utilities for the Relator
emergent electron-mass numerical package.

All dimensional constants used by the numerical scripts are declared here.  The
core modules import these constants rather than redefining them, which keeps the
calculation auditable and avoids silent inconsistencies between reports.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import sqrt
from typing import Iterable, Sequence

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text


@dataclass(frozen=True)
class PhysicalConstants:
    """Physical constants used by the electron-mass closures.

    The SI values below follow the manuscript convention.  In the revised SI,
    c, e, and h are exact by definition.  The numerical value of ħ is displayed
    at the same precision as the manuscript pipeline.  G is the CODATA 2022
    central value used for the baseline run, and alpha is the run input.
    """

    c: float = 299_792_458.0
    e: float = 1.602_176_634e-19
    ħ: float = 1.054_571_817e-34
    G: float = 6.674_30e-11
    α: float = 7.297_352_5643e-3
    m_e_ref_MeV: float = 0.510_998_950_69
    m_μ_ref_MeV: float = 105.658_3745
    m_τ_ref_MeV: float = 1776.86

    # Relative one-sigma uncertainties used in the uncertainty report.
    u_rel_c: float = 0.0
    u_rel_e: float = 0.0
    u_rel_ħ: float = 0.0
    u_rel_G: float = 2.2e-5
    u_rel_α: float = 1.5e-10

    @property
    def ℓP(self) -> float:
        """Planck length, ℓ_P = sqrt(ħ G / c^3), in metres."""
        return sqrt(self.ħ * self.G / self.c**3)

    @property
    def tP(self) -> float:
        """Planck time, t_P = ℓ_P / c, in seconds."""
        return self.ℓP / self.c

    @property
    def mP(self) -> float:
        """Planck mass in the convention m_P = ħ / (c ℓ_P)."""
        return self.ħ / (self.c * self.ℓP)

    @property
    def m_e_ref_kg(self) -> float:
        """Reference electron mass in kg from the displayed MeV rest energy."""
        return self.m_e_ref_MeV * 1.0e6 * self.e / self.c**2

    def with_updates(self, **kwargs: float) -> "PhysicalConstants":
        """Return a copy with updated numeric constants."""
        return replace(self, **kwargs)


CONSTANTS = PhysicalConstants()


def format_float(value: float, digits: int = 15) -> str:
    """Return a compact scientific/decimal representation for report tables."""
    if value == 0:
        return "0"
    abs_v = abs(value)
    if 1e-4 <= abs_v < 1e6:
        return f"{value:.{digits}g}"
    return f"{value:.{digits}e}"


def eV_from_MeV(value_MeV: float) -> float:
    """Convert MeV to eV."""
    return 1.0e6 * value_MeV


def ppm(relative_value: float) -> float:
    """Convert a relative value to parts per million."""
    return 1.0e6 * relative_value


def ppb(relative_value: float) -> float:
    """Convert a relative value to parts per billion."""
    return 1.0e9 * relative_value


def make_console() -> Console:
    """Create a Rich console used by all report scripts."""
    return Console(width=132)


def print_header(console: Console, title: str, subtitle: str | None = None) -> None:
    """Print a professional report header."""
    text = Text(title, style="bold cyan")
    if subtitle:
        text.append("\n" + subtitle, style="white")
    console.print(Panel(text, expand=False, border_style="cyan"))


def rich_table(
    title: str,
    columns: Sequence[str],
    rows: Iterable[Sequence[object]],
    *,
    caption: str | None = None,
) -> Table:
    """Build a Rich table with consistent styling."""
    table = Table(title=title, caption=caption, show_lines=False, header_style="bold magenta")
    for col in columns:
        table.add_column(str(col), overflow="fold")
    for row in rows:
        table.add_row(*[str(item) for item in row])
    return table


def constants_table(constants: PhysicalConstants = CONSTANTS) -> Table:
    """Table of constants used in the run."""
    rows = [
        ("c", format_float(constants.c, 16), "m s⁻¹", "SI exact"),
        ("e", format_float(constants.e, 16), "C", "SI exact"),
        ("ħ", format_float(constants.ħ, 16), "J s", "SI via fixed h, displayed precision"),
        ("G", format_float(constants.G, 16), "m³ kg⁻¹ s⁻²", "CODATA 2022 baseline used in manuscript"),
        ("α", format_float(constants.α, 16), "dimensionless", "fine-structure input"),
        ("m_e c²", format_float(constants.m_e_ref_MeV, 16), "MeV", "reference rest energy"),
        ("ℓ_P", format_float(constants.ℓP, 16), "m", "sqrt(ħG/c³)"),
        ("m_P", format_float(constants.mP, 16), "kg", "ħ/(cℓ_P)"),
    ]
    return rich_table("Constants used by the run", ("Symbol", "Value", "Unit", "Meaning"), rows)
