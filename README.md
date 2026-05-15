# Emergent Electron Mass from Two-Space Boundary -- Numerical Code

Author: Mehrdad Pajuhaan

Repository: https://github.com/pajuhaan/Electron-Mass-Emergent  
Zenodo concept DOI: https://doi.org/10.5281/zenodo.17219278  
Zenodo release DOI: https://doi.org/10.5281/zenodo.17219279  
Paper: https://www.researchgate.net/publication/396387730_Emergent_Electron_Mass_from_Two-Space_Boundary

This repository contains the cleaned numerical implementation for the current Relator electron-mass calculation. The code evaluates the two closed rest-frame electron-mass paths used in the manuscript, with the current Relator scalar Coulomb branch in Path A and the finite Ward vector/tension completion in Path B.

All formulas in this README are written with explicit mathematical notation rather than the custom LaTeX macros used in the paper. For example, the paper macro for the Planck mass is written here explicitly as $m_P$, the speed of light as $c$, the reduced Planck constant as $\hbar$, the fine-structure constant as $\alpha$, and the Planck length as $\ell_P$.

## Abstract

We formulate a candidate rest-frame closure for the electron branch as the stationary equilibrium of a locked Gaussian ring satisfying $R\Omega=c$ in a two-space kinematics, with generator space $\mathcal C$ separated from propagation space $\mathbb R^3$. Two independent closures are evaluated on the same locked shape data. Path A is the scalar vacuum-in-volume model, RLVM, while Path B is the tangential vector-tension model, RLTM.

The only microscopic length inserted into the closure prefactors is $\ell_P=\sqrt{\hbar G/c^3}$, and all remaining factors are dimensionless functions of the locked geometry, $\alpha$, and the assigned scalar or vector completion. The corresponding Planck mass is $m_P=\sqrt{\hbar c/G}=\hbar/(c\ell_P)$.

For each path the stationary radius is obtained in closed form, and the branch mass is read from $m_\bullet=D\hbar/(cR_\bullet)=m_P\,\mathcal C_\bullet(\alpha)\,\exp[-\pi/(8\alpha)]$, with $\bullet\in\{\mathrm A,\mathrm B\}$.

In the current numerical run the two rest-energy equivalents are $m_Ac^2=0.5109977017083316\ \mathrm{MeV}$ and $m_Bc^2=0.5109977017184913\ \mathrm{MeV}$.

Relative to the CODATA electron rest-energy reference $m_ec^2=0.51099895069\ \mathrm{MeV}$, these undershoot by $-2.44419615096268\ \mathrm{ppm}$ and $-2.444176268804592\ \mathrm{ppm}$, respectively. These offsets are below the $11\ \mathrm{ppm}$ one-sigma uncertainty propagated from the present uncertainty of $G$.

Equating the two closed masses gives the scale-free diagnostic $\widehat{\alpha}=0.00729735256488327037$, with $(\widehat{\alpha}-\alpha)/\alpha=0.0799290435808\ \mathrm{ppb}$ at the run input $\alpha=0.0072973525643$.

Since the Planck prefactor gives $m_\bullet\propto(\hbar c/G)^{1/2}$ at fixed dimensionless brackets, the same formulas can be algebraically inverted to define an internal $G$-diagnostic from $(m_e,\alpha)$. The current geometric mean is $\widehat{G}_{\mathrm{geo}}^{(e)}=6.6742673735758264\times10^{-11}\ \mathrm{m^3\,kg^{-1}\,s^{-2}}$.

This inversion is a consistency check, not a measurement of $G$. The construction therefore provides a reproducible rest-frame audit of electron inertia within the stated Relator closure assumptions.

## Implemented physics

### Shared Relator lock

The shared locked geometry is built from the two-space boundary condition $R\Omega=c$. The resolved closure radius is denoted by $R$, while the observable rest-clock radius is $r_{\mathrm{lock}}=R/D$. For a branch mass $m_\bullet$, the readout relation is $m_\bullet=D\hbar/(cR_\bullet)$.

The geometry root is $y_\star$, the transverse scale is $\rho_\star=1/D$, and the locked ratio is $x=\rho_\star/y_\star$. The current run uses $D=2.668305228471116$, $\rho_\star=0.3747697187450251$, $y_\star=1.412858872914209$, and $x=0.2652563011987268$.

### Path A, RLVM

Path A is the Relator-Locked Vacuum-in-Volume Model. Its closure energy is $U_A(R)=a_A/R+b_A\exp(-\sigma_A/\alpha)R^3$.

The electromagnetic coefficient is $a_A=(4\pi\alpha\hbar c/\rho_\star)K_{\mathrm{EM,eff}}^{A}$.

The volume coefficient before the exponential is $b_A=(\hbar c/\ell_P^4)N_{\mathrm{eff,eff}}^{A}\zeta_C^4\rho_\star^2/8$, with $\zeta_C=\pi$.

The scalar count is $N_{\mathrm{eff,eff}}^{A}=N_{\mathrm{eff,raw}}^{A}[1+\rho_\star^2D_C(\alpha)]$.

The Path A local electromagnetic kernel is $K_{\mathrm{EM,eff}}^{A}=K_{\mathrm{EM,raw}}\,s_G(\kappa_T)[1-\rho_\star^2\zeta_{\mathrm{soft}}/2]$.

The current Path A dilation exponent is $\sigma_A=\pi/2$. The stationary radius is $R_A=\ell_P\exp[\sigma_A/(4\alpha)]\{32\pi\alpha K_{\mathrm{EM,eff}}^{A}/[3N_{\mathrm{eff,eff}}^{A}\zeta_C^4\rho_\star^3]\}^{1/4}$.

Equivalently, the mass can be written as $m_A=m_P\,D\{3N_{\mathrm{eff,eff}}^{A}\zeta_C^4\rho_\star^3\exp(-\sigma_A/\alpha)/[32\pi\alpha K_{\mathrm{EM,eff}}^{A}]\}^{1/4}$.

Path A uses the current Relator scalar Coulomb block `DC(alpha)` only in `N_eff_eff_A`. It does not insert the completed Path B vector representative into the Path A scalar slot.

### Path B, RLTM

Path B is the Relator-Locked Tension Model. Its closure energy is $U_B(R)=a_B/R+2\pi T R$.

The electromagnetic coefficient is $a_B=(4\pi\alpha\hbar c/\rho_\star)K_{\mathrm{EM,eff}}^{B}$.

The microscopic tension is $T=(\hbar c/\ell_P^2)K_{T,\mathrm{eff}}^{B}\exp(-\sigma_B/\alpha)[1+f_{M1}(y_\star,\alpha)]$.

The current Path B dilation exponent is $\sigma_B=\pi/4$. The tension weight is $G_{\mathrm{ten}}(\alpha)=\exp(-\sigma_B/\alpha)$.

The finite Ward sources are $\delta_{\mathrm{loc}}=\rho_\star^2\mathcal Z_{\mathrm{ring}}(y_\star)$ and $\delta_{\mathrm{bulk}}=\rho_\star^2\zeta_B(\alpha)$.

The finite Ward tension kernel is $K_{T,\mathrm{eff}}^{B}=K_{T,\mathrm{raw}}^{B}(1+\delta_{\mathrm{loc}})^{-1/2}(1+\delta_{\mathrm{bulk}})^{-1/2}$.

The finite Ward electromagnetic kernel is $K_{\mathrm{EM,eff}}^{B}=K_{\mathrm{EM,raw}}(1+\delta_{\mathrm{loc}})^{1/4}(1+\delta_{\mathrm{bulk}})^{1/4}$.

The exact Ward-neutral identity checked by the code is $(1+\delta)^{-1/2}[(1+\delta)^{1/4}]^2=1$ for each finite Ward source.

The stationary radius is $R_B=\ell_P\{2\alpha K_{\mathrm{EM,eff}}^{B}/[\rho_\star K_{T,\mathrm{eff}}^{B}G_{\mathrm{ten}}(\alpha)(1+f_{M1})]\}^{1/2}$.

Equivalently, the mass can be written as $m_B=m_P\,D\{\rho_\star K_{T,\mathrm{eff}}^{B}G_{\mathrm{ten}}(\alpha)(1+f_{M1})/[2\alpha K_{\mathrm{EM,eff}}^{B}]\}^{1/2}$.

Path B uses finite Ward factors only. It does not use the Path A scalar Coulomb block `DC(alpha)` in its baseline vector/tension closure.

## Notation map

| Code name | Explicit mathematical notation | Meaning |
|---|---:|---|
| `alpha` | $\alpha$ | Fine-structure input used by the electron-mass run |
| `c` | $c$ | Speed of light |
| `e` | $e$ | Elementary charge magnitude |
| `hbar` | $\hbar$ | Reduced Planck constant |
| `G` | $G$ | Newton coupling used in $\ell_P$ |
| `ell_P` | $\ell_P=\sqrt{\hbar G/c^3}$ | Planck length |
| `m_P` | $m_P=\sqrt{\hbar c/G}=\hbar/(c\ell_P)$ | Planck mass |
| `m_e` | $m_e$ | Reference electron mass used only for comparison and $G$ inversion |
| `D` | $D=4F_{\mathrm{struct}}(y_\star)/3$ | Closure-to-rest-clock scale factor |
| `ystar` | $y_\star$ | Shared shape-stationarity root |
| `rho_star` | $\rho_\star=1/D$ | Locked transverse ratio |
| `x` | $x=\rho_\star/y_\star$ | Shared shape ratio used in kernels |
| `kappa_T` | $\kappa_T=x^2/(1+x^2)$ | Tangential second-moment fraction |
| `Kring` | $K_{\mathrm{ring}}(y_\star)$ | Near-field scalar ring harmonic |
| `Lring` | $L_{\mathrm{ring}}(y_\star)$ | First near-field ring harmonic |
| `Gring` | $G_{\mathrm{ring}}=K_{\mathrm{ring}}-L_{\mathrm{ring}}$ | Ring gap diagnostic, not used in `ZigmaRing` |
| `ZigmaRing` | $\mathcal Z_{\mathrm{ring}}=K_{\mathrm{ring}}L_{\mathrm{ring}}/(2\pi^2)$ | Near-edge geometric Ward source |
| `Fstruct` | $F_{\mathrm{struct}}=2+\mathcal Z_{\mathrm{ring}}+\alpha K_{\mathrm{ring}}/(8\pi^2)$ | Shared structure factor |
| `beta_rel` | $\beta_{\mathrm{rel}}=8/\pi$ | Internal Relator dilation slope in the current geometry |
| `sigma_A` | $\sigma_A=4/\beta_{\mathrm{rel}}=\pi/2$ | Path A dilation exponent |
| `sigma_B` | $\sigma_B=2/\beta_{\mathrm{rel}}=\pi/4$ | Path B dilation exponent |
| `K_EM_raw` | $K_{\mathrm{EM,raw}}(x)$ | Raw along-ring electromagnetic kernel |
| `s_G` | $s_G(\kappa_T)$ | Path A Gaussian Maxwell softening |
| `K_EM_mid_A` | $K_{\mathrm{EM,mid}}^{A}=K_{\mathrm{EM,raw}}s_G$ | Path A intermediate electromagnetic kernel |
| `K_EM_eff_A` | $K_{\mathrm{EM,eff}}^{A}$ | Path A effective electromagnetic kernel |
| `N_eff_raw_A` | $N_{\mathrm{eff,raw}}^{A}$ | Path A raw scalar channel count |
| `N_eff_eff_A` | $N_{\mathrm{eff,eff}}^{A}$ | Path A scalar count corrected by $D_C(\alpha)$ |
| `a_A` | $a_A$ | Path A $1/R$ coefficient, units $\mathrm{J\,m}$ |
| `b_A` | $b_A$ | Path A $R^3$ coefficient before the dilation exponential, units $\mathrm{J\,m^{-3}}$ |
| `R_A` | $R_A$ | Path A stationary closure radius |
| `rlock_A` | $r_{\mathrm{lock},A}=R_A/D$ | Path A observable rest-clock radius |
| `m_A` | $m_A$ | Path A branch mass |
| `K_T_raw_B` | $K_{T,\mathrm{raw}}^{B}=2\kappa_T$ | Path B raw tension kernel |
| `delta_loc` | $\delta_{\mathrm{loc}}=\rho_\star^2\mathcal Z_{\mathrm{ring}}$ | Path B local finite Ward source |
| `delta_bulk` | $\delta_{\mathrm{bulk}}=\rho_\star^2\zeta_B$ | Path B bulk vector finite Ward source |
| `K_T_eff_B` | $K_{T,\mathrm{eff}}^{B}$ | Path B finite Ward tension kernel |
| `K_EM_eff_B` | $K_{\mathrm{EM,eff}}^{B}$ | Path B finite Ward electromagnetic kernel |
| `fM1` | $f_{M1}=(\alpha/2)y_\star^2\exp(-y_\star^2/2)$ | Path B mode-selective tangential remainder |
| `Gten` | $G_{\mathrm{ten}}=\exp(-\sigma_B/\alpha)$ | Path B tension dilation weight |
| `a_B` | $a_B$ | Path B $1/R$ coefficient, units $\mathrm{J\,m}$ |
| `T` | $T$ | Path B string tension, units $\mathrm{J\,m^{-1}}$ |
| `R_B` | $R_B$ | Path B stationary closure radius |
| `rlock_B` | $r_{\mathrm{lock},B}=R_B/D$ | Path B observable rest-clock radius |
| `m_B` | $m_B$ | Path B branch mass |
| `DC(alpha)` | $D_C(\alpha)$ | Current Relator scalar Coulomb branch used only in Path A baseline |
| `DCQEDALP(alpha)` | $D_{C}^{\mathrm{QED+ALP}}(\alpha)$ | Reference-only scalar diagnostic, not part of baseline Path A/B |
| `Kov` | $K_{\mathrm{ov}}=(150\pi^2-8\pi^4-315)/(180\pi^6)$ | Current angular overlap coefficient from the alpha-sector module |
| `Lambda_B_th` | $\Lambda_B^{(\mathrm{th})}$ | Current reduced vector representative used for Path B |
| `zeta_B` | $\zeta_B=K_{\mathrm{ov}}\Lambda_B^{(\mathrm{th})}/(2\pi^2)$ | Path B bulk vector overlap |
| `Lambda_ind` | $\Lambda_{\mathrm{ind}}=\ln(8\sqrt{\pi})-2$ | Ring-local inductive logarithm |
| `zeta_soft` | $\zeta_{\mathrm{soft}}=K_{\mathrm{ov}}\Lambda_{\mathrm{ind}}/(2\pi^2)$ | Path A local soft electromagnetic correction only |
| `R_moth` | $R_{\mathrm{moth}}(D)$ | Scalar mother-law radicand used in $D_C(\alpha)$ |
| `Phi_dyn` | $\Phi_{\mathrm{dyn}}(D)$ | Rank-5 scalar dynamic response |

## Current numerical benchmark

### Constants used by the run

| Quantity | Value | Unit |
|---|---:|---|
| $c$ | $2.9979245800000000\times10^8$ | $\mathrm{m\,s^{-1}}$ |
| $e$ | $1.6021766339999999\times10^{-19}$ | $\mathrm C$ |
| $\hbar$ | $1.0545718170000000\times10^{-34}$ | $\mathrm{J\,s}$ |
| $G$ | $6.6742999999999994\times10^{-11}$ | $\mathrm{m^3\,kg^{-1}\,s^{-2}}$ |
| $\alpha$ | $0.0072973525643$ | dimensionless |
| $m_ec^2$ | $0.51099895069$ | $\mathrm{MeV}$ |
| $\ell_P$ | $1.6162550239285500\times10^{-35}$ | $\mathrm m$ |
| $m_P$ | $2.1764343420511269\times10^{-8}$ | $\mathrm{kg}$ |

### Shared locked geometry

| Quantity | Value |
|---|---:|
| $y_\star$ | $1.412858872914209$ |
| $K_{\mathrm{ring}}(y_\star)$ | $0.4156324192936439$ |
| $L_{\mathrm{ring}}(y_\star)$ | $0.05653957686872045$ |
| $G_{\mathrm{ring}}(y_\star)$ | $0.3590928424249234$ |
| $\mathcal Z_{\mathrm{ring}}(y_\star)$ | $0.001190507753137068$ |
| $F_{\mathrm{struct}}(y_\star)$ | $2.001228921353337$ |
| $D$ | $2.668305228471116$ |
| $\rho_\star$ | $0.3747697187450251$ |
| $x$ | $0.2652563011987268$ |
| $\kappa_T$ | $0.06573568314719434$ |
| $\beta_{\mathrm{rel}}$ | $2.546479089470325$ |
| $\sigma_A$ | $1.570796326794897$ |
| $\sigma_B$ | $0.7853981633974484$ |

### Alpha-sector blocks used by the current electron-mass run

| Quantity | Value | Role |
|---|---:|---|
| $K_{\mathrm{ov}}$ | $0.002231538916531971$ | angular overlap coefficient |
| $\Lambda_B^{(\mathrm{th})}$ | $0.691683146106991$ | reduced vector representative |
| $\zeta_B$ | $7.8195528195469194\times10^{-5}$ | Path B bulk vector overlap |
| $\Lambda_{\mathrm{ind}}$ | $0.6518064846045357$ | ring-local inductive logarithm |
| $\zeta_{\mathrm{soft}}$ | $7.3687428458752324\times10^{-5}$ | Path A local EM correction |
| $D_C(\alpha)$ | $0.002315457831882793$ | current scalar Coulomb branch |
| $R_{\mathrm{moth}}[D_C(\alpha)]$ | $0.9936715126327368$ | scalar mother-law radicand |
| $\Phi_{\mathrm{dyn}}[D_C(\alpha)]$ | $5.70992324488732$ | rank-5 scalar dynamic response |

### Path A, RLVM outputs

| Quantity | Value | Unit |
|---|---:|---|
| $K_{\mathrm{EM,raw}}$ | $0.249071299300973$ | dimensionless |
| $s_G$ | $0.9699360322089555$ | dimensionless |
| $K_{\mathrm{EM,mid}}^{A}$ | $0.2415832277811149$ | dimensionless |
| $K_{\mathrm{EM,eff}}^{A}$ | $0.2415819776396209$ | dimensionless |
| $N_{\mathrm{eff,raw}}^{A}$ | $2.126866130533706$ | dimensionless |
| $N_{\mathrm{eff,eff}}^{A}$ | $2.127557811806196$ | dimensionless |
| $a_A$ | $1.8688398940229415\times10^{-27}$ | $\mathrm{J\,m}$ |
| $b_A$ | $1.6856879087891495\times10^{114}$ | $\mathrm{J\,m^{-3}}$ |
| $R_A$ | $1.0303933101755780\times10^{-12}$ | $\mathrm m$ |
| $r_{\mathrm{lock},A}$ | $3.8616021105125673\times10^{-13}$ | $\mathrm m$ |
| $m_Ac^2$ | $0.5109977017083316$ | $\mathrm{MeV}$ |

### Path B, RLTM outputs

| Quantity | Value | Unit |
|---|---:|---|
| $K_{\mathrm{EM,raw}}$ | $0.249071299300973$ | dimensionless |
| $K_{T,\mathrm{raw}}^{B}$ | $0.1314713662943887$ | dimensionless |
| $\delta_{\mathrm{loc}}$ | $0.0001672096022022919$ | dimensionless |
| $\delta_{\mathrm{bulk}}$ | $1.0982745075879497\times10^{-5}$ | dimensionless |
| $(1+\delta_{\mathrm{loc}})^{-1/2}$ | $0.9999164056820823$ | dimensionless |
| $(1+\delta_{\mathrm{bulk}})^{-1/2}$ | $0.9999945086726945$ | dimensionless |
| $(1+\delta_{\mathrm{loc}})^{1/4}$ | $1.000041799779645$ | dimensionless |
| $(1+\delta_{\mathrm{bulk}})^{1/4}$ | $1.000002745674961$ | dimensionless |
| $K_{T,\mathrm{eff}}^{B}$ | $0.1314596541432451$ | dimensionless |
| $K_{\mathrm{EM,eff}}^{B}$ | $0.2490823943238152$ | dimensionless |
| $f_{M1}$ | $0.002684541055181403$ | dimensionless |
| $G_{\mathrm{ten}}$ | $1.8106344293342899\times10^{-47}$ | dimensionless |
| $a_B$ | $1.9268619288542303\times10^{-27}$ | $\mathrm{J\,m}$ |
| $T$ | $0.0002888448754525408$ | $\mathrm{J\,m^{-1}}$ |
| $R_B$ | $1.0303933101550915\times10^{-12}$ | $\mathrm m$ |
| $r_{\mathrm{lock},B}$ | $3.8616021104357901\times10^{-13}$ | $\mathrm m$ |
| $m_Bc^2$ | $0.5109977017184913$ | $\mathrm{MeV}$ |

### Mass and diagnostic comparison

| Case | $mc^2$ | Deviation from $m_ec^2$ | Relative deviation |
|---|---:|---:|---:|
| Path A, RLVM | $0.5109977017083316\ \mathrm{MeV}$ | $-1.248981668422466\ \mathrm{eV}$ | $-2.44419615096268\ \mathrm{ppm}$ |
| Path B, RLTM | $0.5109977017184913\ \mathrm{MeV}$ | $-1.248971508660546\ \mathrm{eV}$ | $-2.444176268804592\ \mathrm{ppm}$ |
| Path A minus Path B |  | $-1.0159761920647270\times10^{-5}\ \mathrm{eV}$ |  |

The path-equality diagnostic gives $\widehat{\alpha}=0.00729735256488327037$, with $(\widehat{\alpha}-\alpha)/\alpha=0.0799290435808\ \mathrm{ppb}$.

The internal geometric mean from the $G$ inversion is $\widehat{G}_{\mathrm{geo}}^{(e)}=6.6742673735758264\times10^{-11}\ \mathrm{m^3\,kg^{-1}\,s^{-2}}$. Using this value in the same formulas gives $m_Ac^2=0.5109989506849203\ \mathrm{MeV}$ and $m_Bc^2=0.5109989506950799\ \mathrm{MeV}$.

## Diagnostic branches and channel separation

The baseline code keeps the channels separated.

- Path A scalar slot uses $D_C(\alpha)$ through $N_{\mathrm{eff,eff}}^{A}=N_{\mathrm{eff,raw}}^{A}[1+\rho_\star^2D_C(\alpha)]$.
- Path A local EM softening uses $\zeta_{\mathrm{soft}}=K_{\mathrm{ov}}\Lambda_{\mathrm{ind}}/(2\pi^2)$ only in $K_{\mathrm{EM,eff}}^{A}$.
- Path B uses $\zeta_B=K_{\mathrm{ov}}\Lambda_B^{(\mathrm{th})}/(2\pi^2)$ only through finite Ward factors.
- Path B does not use $D_C(\alpha)$ in the baseline vector/tension closure.
- The QED+ALP scalar diagnostic $D_C^{\mathrm{QED+ALP}}(\alpha)$ is implemented in a separate report only. It does not modify the baseline RLVM mass, the baseline RLTM mass, or the charged-lepton hierarchy report.

## Files

- `relator_electron/common.py`: constants, SI units, and table helpers.
- `relator_electron/alpha_sector.py`: current `DC(alpha)`, `zetaB`, `zetaSoft`, and alpha-sector provenance constants.
- `relator_electron/geometry.py`: finite-width `S1` geometry, shared shape root, ring harmonics, and raw kernels.
- `relator_electron/rlvm_core.py`: Path A, RLVM, core calculation.
- `relator_electron/rltm_core.py`: Path B, RLTM, core calculation.
- `relator_electron/pipeline.py`: assembled two-path pipeline and alpha-equality diagnostic.
- `relator_electron/qed_alp_dc.py`: reference-only QED+ALP scalar diagnostic branch.
- `report_main_outputs.py`: main constants, geometry, alpha-sector blocks, Path A outputs, Path B outputs, and path-equality diagnostic.
- `report_uncertainty_sensitivity.py`: uncertainty propagation and local sensitivity tables.
- `report_emergent_G.py`: internal inferred-$G$ consistency check.
- `report_lepton_hierarchy.py`: charged-lepton ladder lift from the current Path A anchor.
- `report_qed_alp_dc.py`: separate reference-only QED+ALP scalar-branch diagnostic.
- `run_all_reports.py`: compact all-in-one summary from the same core calculations.

## Run

Install the Python dependencies if needed.

```bash
python3 -m pip install scipy rich
```

Run the reports from the repository root.

```bash
cd relator_electron_code

python3 report_main_outputs.py
python3 report_uncertainty_sensitivity.py
python3 report_emergent_G.py
python3 report_lepton_hierarchy.py
python3 report_qed_alp_dc.py
python3 run_all_reports.py
```

The calculations are deterministic. There is no stochastic sampling and no fitted microscopic parameter in the two baseline mass closures.

## Reproducibility notes

The baseline electron-mass pipeline uses the run input $\alpha=0.0072973525643$ and $G=6.6742999999999994\times10^{-11}\ \mathrm{m^3\,kg^{-1}\,s^{-2}}$. The constants $c$, $e$, and $h$ are exact in SI-2019. The reference electron rest energy $m_ec^2=0.51099895069\ \mathrm{MeV}$ is used only for comparison rows and for the separate algebraic $G$ inversion.

The core mass formulas use only $\alpha$, $G$ through $\ell_P$, and dimensionless kernels fixed by the locked geometry and the specified scalar or vector completion. The measured electron mass is not used to compute the baseline Path A or Path B masses.