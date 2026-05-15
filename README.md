# Relator electron-mass numerical code

Author: Mehrdad Pajuhaan

Repository: https://github.com/pajuhaan/Electron-Mass-Emergent  
Zenodo: https://doi.org/10.5281/zenodo.17219278  
Paper: https://www.researchgate.net/publication/396387730_Emergent_Electron_Mass_from_Two-Space_Boundary

This repository contains the cleaned numerical implementation for the current Relator electron-mass calculation.

## Abstract

We formulate a candidate rest-frame closure for the electron branch as the stationary equilibrium of a locked Gaussian ring satisfying $R\Omega=c$ in a two-space kinematics, with generator space $\mathcal C$ separated from propagation space $\mathbb R^3$. Two independent closures are evaluated on the same locked shape data. Path A is the scalar vacuum-in-volume model, RLVM, while Path B is the tangential vector-tension model, RLTM. The only microscopic length inserted into the closure prefactors is

$$
\ell_P=\sqrt{\frac{\hbar G}{c^3}},
$$

and all remaining factors are dimensionless functions of the locked geometry, $\alpha$, and the assigned scalar or vector completion. For each path the stationary radius is obtained in closed form, and the branch mass is read from

$$
m_\bullet
=
\frac{D\hbar}{cR_\bullet}
=
m_P\,\mathcal C_\bullet(\alpha)\,
\exp\!\left(-\frac{\pi}{8\alpha}\right),
\qquad
m_P=\sqrt{\frac{\hbar c}{G}},
\qquad
\bullet\in\{\mathrm A,\mathrm B\}.
$$

In the current numerical run the two rest-energy equivalents are

$$
m_Ac^2
=
0.5109977017083316\ \mathrm{MeV},
\qquad
m_Bc^2
=
0.5109977017184913\ \mathrm{MeV}.
$$

Relative to the CODATA electron rest-energy reference, these undershoot by

$$
-2.44419615096268\ \mathrm{ppm}
\qquad\text{and}\qquad
-2.444176268804592\ \mathrm{ppm},
$$

respectively, below the $11\ \mathrm{ppm}$ one-sigma uncertainty propagated from the present uncertainty of $G$. Equating the two closed masses gives the scale-free diagnostic

$$
\widehat{\alpha}
=
0.00729735256488327037,
\qquad
\frac{\widehat{\alpha}-\alpha}{\alpha}
=
0.0799290435808\ \mathrm{ppb}.
$$

Since the Planck prefactor gives $m_\bullet\propto(\hbar c/G)^{1/2}$ at fixed dimensionless brackets, the same formulas can be algebraically inverted to define an internal $G$-diagnostic from $(m_e,\alpha)$,

$$
\widehat{G}_{\mathrm{geo}}^{(e)}
=
6.6742673735758264\times10^{-11}\ 
\mathrm{m^3\,kg^{-1}\,s^{-2}}.
$$

This inversion is a consistency check, not a measurement of $G$. The construction therefore provides a reproducible rest-frame audit of electron inertia within the stated Relator closure assumptions.

## Implemented physics

- Path A, RLVM, uses the current Relator scalar Coulomb block `DC(alpha)` only in `N_eff_eff_A`.
- Path A keeps the scalar channel scalar and does not insert a combined scalar-vector bridge.
- Path B, RLTM, uses finite Ward factors

$$
K_{T,\mathrm{eff}}
=
K_{T,\mathrm{raw}}
(1+\delta_{\mathrm{loc}})^{-1/2}
(1+\delta_{\mathrm{bulk}})^{-1/2},
$$

$$
K_{\mathrm{EM},\mathrm{eff}}
=
K_{\mathrm{EM},\mathrm{raw}}
(1+\delta_{\mathrm{loc}})^{1/4}
(1+\delta_{\mathrm{bulk}})^{1/4}.
$$

## Files

- `relator_electron/common.py`: constants, units, and table helpers.
- `relator_electron/alpha_sector.py`: `DC`, `zetaB`, and `zetaSoft`.
- `relator_electron/geometry.py`: finite-width `S1` geometry and shared shape root.
- `relator_electron/rlvm_core.py`: Path A, RLVM, core calculation.
- `relator_electron/rltm_core.py`: Path B, RLTM, core calculation.
- `relator_electron/pipeline.py`: assembled two-path pipeline and alpha-equality diagnostic.
- `report_main_outputs.py`: main intermediate values and final masses.
- `report_uncertainty_sensitivity.py`: uncertainty and sensitivity tables.
- `report_emergent_G.py`: inferred-$G$ cross-check.
- `report_lepton_hierarchy.py`: charged-lepton ladder lift.
- `report_qed_alp_dc.py`: reference-only QED+ALP scalar diagnostic.
- `run_all_reports.py`: compact all-in-one summary.

## Run

```bash
cd relator_electron_code

python3 report_main_outputs.py
python3 report_uncertainty_sensitivity.py
python3 report_emergent_G.py
python3 report_lepton_hierarchy.py
python3 report_qed_alp_dc.py
python3 run_all_reports.py