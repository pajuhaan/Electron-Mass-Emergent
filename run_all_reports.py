#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mehrdad Pajuhaan

Run all Relator electron-mass reports in sequence.

Each report is launched as a fresh Python process.  The PYTHONUNBUFFERED setting
keeps progress and Rich tables visible in deterministic order.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


REPORT_FILES: tuple[str, ...] = (
    "report_main_outputs.py",
    "report_uncertainty_sensitivity.py",
    "report_emergent_G.py",
    "report_lepton_hierarchy.py",
    "report_qed_alp_dc.py",
)


def main() -> None:
    """Run all report scripts in deterministic order."""
    root = Path(__file__).resolve().parent
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    for filename in REPORT_FILES:
        print("\n" + "=" * 120, flush=True)
        print(f"Running {filename}", flush=True)
        print("=" * 120, flush=True)
        subprocess.run([sys.executable, str(root / filename)], cwd=root, env=env, check=True)


if __name__ == "__main__":
    main()
