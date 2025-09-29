# Electron Mass from Relator Kinematics — Reproducible Code

Minimal, clean code to reproduce the numerics and figures for **Electron Mass from Relator Kinematics**.  
Planck-only scaling, two independent closures (Path A/B), and a geometry-only fine-structure identity — **no fits**.

**Manuscript (PDF):** https://doi.org/10.5281/zenodo.17219278

---

## Paths at a glance

- **Path A — RLVM (Vacuum-in-Volume closure)**  
  Local EM softening only; scalar S4 completion; dilation weight \(e^{-\sigma_A/\alpha}\) with \(\sigma_A=\pi/2\).  
  **Mass formula**  
  \[
  m_A \;=\; \frac{D\,\hbar}{c\,\ell_P}\,
  \Bigg[\frac{3\,N_{\mathrm{eff}}^{\mathrm{eff}}\;\zeta_C^{\,4}\;\rho^{3}\;e^{-\sigma_A/\alpha}}
  {32\pi\,\alpha\;K_{\mathrm{EM}}^{\mathrm{eff}}}\Bigg]^{\!1/4}.
  \]

- **Path B — RLTM (Tangential Tension closure)**  
  Ward-preserving product kernels; includes the mode-\(M{=}1\) remainder \(f_{M1}\); weight \(e^{-\sigma_B/\alpha}\) with \(\sigma_B=\pi/4\).  
  **Mass formula**  
  \[
  m_B \;=\; \frac{D\,\hbar}{c\,\ell_P}\,
  \sqrt{\frac{\rho\;K_T^{\mathrm{eff}}\;e^{-\sigma_B/\alpha}\,(1+f_{M1})}
  {2\,\alpha\;K_{\mathrm{EM}}^{\mathrm{raw,eff}}}},
  \qquad
  f_{M1}=\frac{\alpha}{2}\,e^{-y^{*2}/2}\,y^{*2}.
  \]

**Shared definitions (used by both paths)**  
\(D=\tfrac{4}{3}S(y^*)\), \(\rho=1/D\), \(x=\rho/y^*\);  
\(S(y)=2+\zeta_{\mathrm{geom}}(y)+\alpha K_0(y)/(8\pi^2)\), \(\zeta_{\mathrm{geom}}(y)=K_0(y)\Lambda(y)/(2\pi^2)\);  
\(K_{\mathrm{EM}}^{\mathrm{raw}}(x)=(2x/\pi)\,K(m=-4x^2)\) (complete elliptic integral), \(K_T=2x^2/(1+x^2)\);  
for Path B,
\(K_T^{\mathrm{eff}}=K_T\big(1-\tfrac12\rho^2\zeta_{\mathrm{geom}}\big)\big(1-\tfrac12\rho^2\zeta(\alpha)\big)\),  
\(K_{\mathrm{EM}}^{\mathrm{raw,eff}}=K_{\mathrm{EM}}^{\mathrm{raw}}\big(1+\tfrac14\rho^2\zeta_{\mathrm{geom}}\big)\big(1+\tfrac14\rho^2\zeta(\alpha)\big)\).

---

## Key results (typical, no tuning)

> Reproduced by the scripts here; exact digits depend only on the chosen constants set and high-precision arithmetic.

| Quantity | Value / Deviation | Notes |
|---|---:|---|
| Emergent electron mass — **Path A** | \(\delta m_A = -2.45046970191\) **ppm** | Undershoot vs. \(m_e\); within the \(\sim 11\) ppm uncertainty set by \(G\). |
| Emergent electron mass — **Path B** | \(\delta m_B = -2.45076007903\) **ppm** | Same envelope. |
| Fine-structure (geometry-only identity) | \(\alpha^{-1} \approx 137.0359991769\ldots\) | From \(q^2=4\pi\varepsilon_0\,\alpha\,\hbar c\); sub-ppb agreement. |
| Gravitational constant (Relator, paper) | \(G_{\rm emergent}\approx 6.67426728776\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}\) | ≈ **−4.9 ppm** vs. CODATA-2022 \(6.67430\times10^{-11}\); within CODATA σ. |

---

## Install and run

**Requirements:** Python ≥ 3.10, [mpmath](https://pypi.org/project/mpmath/). `tabulate` is optional (pretty tables).

```bash
git clone https://github.com/pajuhaan/Electron-Mass-Emergent
cd Electron-Mass-Emergent

# Install dependencies
python -m pip install --upgrade pip
python -m pip install mpmath tabulate

# Main run: two paths (A/B), tables, checksum.txt
python Relator_Electron_Emergent_Mass.py

# Uncertainty propagation (±1σ central differences for constants)
python "Uncertainty budget from physical constants.py"
````

* The main script prints shared geometry, Λ-chain, Path A/B masses, per-term effects, and saves a `checksum.txt` with a SHA-256 of the full report. 
* The uncertainty script perturbs (\hbar,c,e,G,\alpha) and reports (\delta m) in eV/ppm plus log-sensitivities. 

---

## Cite this work

Please cite the Zenodo record:

**Pajuhaan, M. (2025). *Electron Mass from Relator Kinematics*. Zenodo.**
[https://doi.org/10.5281/zenodo.17219279](https://doi.org/10.5281/zenodo.17219279)

**How to cite (BibTeX):**

```bibtex
@misc{Pajuhaan2025ElectronMassRelator,
  author       = {Pajuhaan, Mehrdad},
  title        = {Electron Mass from Relator Kinematics},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17219279},
  url          = {https://doi.org/10.5281/zenodo.17219279}
}
```

**Manuscript (PDF):** [https://doi.org/10.5281/zenodo.17219278](https://doi.org/10.5281/zenodo.17219278)

---

## Reproducibility

* No fitting or post-matching; the only dimensionful input is (\ell_P=\sqrt{\hbar G/c^3}).
* All outputs are deterministic given the repo commit and constants table inside the scripts.

