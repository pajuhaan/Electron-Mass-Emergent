# Relator electron-mass numerical code

Author: Mehrdad Pajuhaan

This directory contains a cleaned, modular numerical implementation of the current electron-mass calculations.

## Physics update implemented

- Path A, RLVM, uses the current Relator scalar Coulomb block `DC(alpha)` in `N_eff_eff_A`.
- Path A does not compute or insert any combined scalar-vector bridge.  The scalar channel remains scalar.
- Path B, RLTM, uses finite Ward factors
  
  `K_T_eff = K_T_raw (1+delta_loc)^(-1/2) (1+delta_bulk)^(-1/2)`
  
  `K_EM_eff = K_EM_raw (1+delta_loc)^(1/4) (1+delta_bulk)^(1/4)`.

## Files

- `relator_electron/common.py` contains constants, units, and Rich-table helpers.
- `relator_electron/alpha_sector.py` contains `DC`, `ζB`, and `ζsoft`.
- `relator_electron/geometry.py` contains the finite-width S1 geometry and shared shape root.
- `relator_electron/rlvm_core.py` contains the Path A core calculation.
- `relator_electron/rltm_core.py` contains the Path B core calculation.
- `relator_electron/pipeline.py` assembles both paths and the alpha-equality diagnostic.
- `report_main_outputs.py` prints all main intermediate values and final masses.
- `report_uncertainty_sensitivity.py` prints uncertainty and sensitivity tables.
- `report_emergent_G.py` prints the inferred-G cross-check.
- `report_lepton_hierarchy.py` prints the charged-lepton hierarchy lift.
- `report_qed_alp_dc.py` prints the separate reference-only QED+ALP scalar-branch diagnostic.
- `run_all_reports.py` prints a compact all-in-one summary from the same core calculations, including the separate QED+ALP diagnostic report.

## Run

```bash
cd relator_electron_code
python3 report_main_outputs.py
python3 report_uncertainty_sensitivity.py
python3 report_emergent_G.py
python3 report_lepton_hierarchy.py
python3 report_qed_alp_dc.py
python3 run_all_reports.py
```

The scripts require `scipy` and `rich`, both of which are available in the current execution environment.
