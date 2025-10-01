#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2025 by Mehrdad Pajuhaan (pajuhaan@gmail.com)
"""
Numeric Uncertainty Propagation (standalone)
- Uses Relator_Electron_Emergent_Mass.py in the same folder
- Central finite-difference at ±1σ (relative) for each constant
- Reports δm_A, δm_B in eV and ppm (relative to baseline masses), plus log-sensitivities
"""

import os, sys, importlib
from math import isfinite
try:
    import mpmath as mp
except Exception:
    raise SystemExit("mpmath is required.")

# ---------------------------------------------------------------------------
# Ensure we can import the sibling module
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

try:
    rel = importlib.import_module("Relator_Electron_Emergent_Mass")
except Exception as e:
    raise SystemExit(f"Cannot import Relator_Electron_Emergent_Mass.py from {HERE}.\n{e}")

# Tiny helpers mirroring the main module’s style
def nstr(x, n=28): return mp.nstr(x, n)
def ppm(x): return mp.mpf('1e6') * x

def _format_table(headers, rows, title=None):
    """Pure-Python table (no external deps)."""
    cols = len(headers)
    widths = [len(str(h)) for h in headers]
    for r in rows:
        for j in range(cols):
            widths[j] = max(widths[j], len(str(r[j])))
    def line(vals): return "  ".join(str(vals[j]).ljust(widths[j]) for j in range(cols))
    out = []
    if title: out.append(title)
    out.append(line(headers))
    out.append(line(["-"*w for w in widths]))
    for r in rows: out.append(line(r))
    return "\n".join(out)

# ---------------------------------------------------------------------------
# Core numeric propagation
# ---------------------------------------------------------------------------
def propagate_one(const_key: str, u_rel: mp.mpf, baseline, alpha0):
    """
    Central finite difference at ±1σ for the given constant.
    Returns dict with deltas for A/B in eV and ppm, and log-sensitivities.
    """
    mA0 = baseline["mA"]; mB0 = baseline["mB"]

    # If uncertainty is zero, return zeros without recomputing
    if u_rel == 0:
        return {
            "const": const_key, "u_rel": u_rel,
            "dA_eV": mp.mpf('0'), "dA_ppm": mp.mpf('0'), "sA": mp.mpf('0'),
            "dB_eV": mp.mpf('0'), "dB_ppm": mp.mpf('0'), "sB": mp.mpf('0'),
        }

    # Helper to compute masses with current rel.PHYS (and alpha argument)
    def _run(alpha_value):
        R = rel.compute_all(alpha_value)
        return R["mA"], R["mB"]

    # Alpha is passed as an argument; others live in rel.PHYS
    if const_key == "alpha_in":
        a_plus  = alpha0 * (1 + u_rel)
        a_minus = alpha0 * (1 - u_rel)
        mA_p, mB_p = _run(a_plus)
        mA_m, mB_m = _run(a_minus)
    else:
        old = rel.PHYS[const_key]["value"]
        try:
            rel.PHYS[const_key]["value"] = old * (1 + u_rel)
            mA_p, mB_p = _run(alpha0)
            rel.PHYS[const_key]["value"] = old * (1 - u_rel)
            mA_m, mB_m = _run(alpha0)
        finally:
            rel.PHYS[const_key]["value"] = old  # restore

    # Central 1σ effect (absolute) in MeV, then to eV
    dA_MeV = abs(mA_p - mA_m) / 2
    dB_MeV = abs(mB_p - mB_m) / 2
    dA_eV  = rel.eV_from_MeV(dA_MeV)
    dB_eV  = rel.eV_from_MeV(dB_MeV)

    # ppm relative to baseline predictions
    dA_ppm = ppm(dA_MeV / mA0)
    dB_ppm = ppm(dB_MeV / mB0)

    # Log-sensitivity d ln m / d ln const (dimensionless)
    sA = (mp.log(mA_p) - mp.log(mA_m)) / (2 * u_rel)
    sB = (mp.log(mB_p) - mp.log(mB_m)) / (2 * u_rel)

    return {
        "const": const_key, "u_rel": u_rel,
        "dA_eV": dA_eV, "dA_ppm": dA_ppm, "sA": sA,
        "dB_eV": dB_eV, "dB_ppm": dB_ppm, "sB": sB,
    }

def main():
    # Baseline run at PHYS["alpha_in"]["value"]
    alpha0 = rel.PHYS["alpha_in"]["value"]
    base   = rel.compute_all(alpha0)

    mA0 = base["mA"]; mB0 = base["mB"]

    print(_format_table(
        ["Quantity", "Value", "Note"],
        [
            ["alpha (input)", nstr(alpha0, 20), rel.PHYS["alpha_in"]["source"]],
            ["m_A [MeV]",     nstr(mA0, 18), ""],
            ["m_B [MeV]",     nstr(mB0, 18), ""],
        ],
        title="=== BASELINE (from Relator_Electron_Emergent_Mass) ==="
    ))
    print()

    # Constants to propagate (present in rel.PHYS)
    keys = ["hbar", "c", "e", "G", "alpha_in"]
    rows = []
    quad_A = mp.mpf('0')
    quad_B = mp.mpf('0')

    for k in keys:
        if k not in rel.PHYS:
            continue
        u_rel = mp.mpf(rel.PHYS[k].get("u_rel", 0))
        res = propagate_one(k, u_rel, base, alpha0)

        # accumulate quadrature (only finite numbers)
        if isfinite(res["dA_eV"]): quad_A += res["dA_eV"]**2
        if isfinite(res["dB_eV"]): quad_B += res["dB_eV"]**2

        rows.append([
            k,
            nstr(res["u_rel"], 12),
            nstr(res["dA_eV"], 12), nstr(res["dA_ppm"], 12), nstr(res["sA"], 10),
            nstr(res["dB_eV"], 12), nstr(res["dB_ppm"], 12), nstr(res["sB"], 10),
        ])

    print(_format_table(
        ["Constant", "rel σ (1σ)",
         "δm_A [eV]", "δm_A [ppm]", "∂ln m_A/∂ln const",
         "δm_B [eV]", "δm_B [ppm]", "∂ln m_B/∂ln const"],
        rows,
        title="=== NUMERIC UNCERTAINTY PROPAGATION (±1σ central diff) ==="
    ))
    print()

    # Totals (quadrature)
    dA_tot_eV = mp.sqrt(quad_A)
    dB_tot_eV = mp.sqrt(quad_B)
    dA_tot_ppm = ppm(dA_tot_eV / rel.eV_from_MeV(mA0))
    dB_tot_ppm = ppm(dB_tot_eV / rel.eV_from_MeV(mB0))

    print(_format_table(
        ["Path", "Total δm [eV] (RSS)", "Total δm [ppm] (RSS)"],
        [
            ["A", nstr(dA_tot_eV, 12), nstr(dA_tot_ppm, 12)],
            ["B", nstr(dB_tot_eV, 12), nstr(dB_tot_ppm, 12)],
        ],
        title="=== COMBINED (QUADRATURE) 1σ UNCERTAINTY ==="
    ))
    print()

    # Reference: constants table (value & source & u_rel)
    cref_rows = []
    for k in keys:
        if k not in rel.PHYS: continue
        d = rel.PHYS[k]
        cref_rows.append([k, nstr(d["value"], 20), d.get("source",""), nstr(d.get("u_rel",0), 12)])
    print(_format_table(
        ["Constant", "Value", "Source", "rel σ"],
        cref_rows,
        title="=== CONSTANTS USED (from module PHYS) ==="
    ))

if __name__ == "__main__":
    main()
