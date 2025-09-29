# Electron Mass Emergent — Relator Kinematics [Reproducible Code]

Minimal, clean code to reproduce the numerics and figures for **Electron Mass from Relator Kinematics**.  
Planck-only scaling, two independent closures (Path A/B), and a geometry-only fine-structure identity — **no fits**.

**Manuscript (PDF):** https://doi.org/10.5281/zenodo.17219278

---

## Paths at a glance

- **Path A — RLVM (Vacuum-in-Volume closure)**  
  Local EM softening only (from $\Lambda_{\rm ind}$); scalar S4 completion; dilation weight $e^{-\sigma_A/\alpha}$ with $\sigma_A=\pi/2$.

  $$m_A = \frac{D\,\hbar}{c\,\ell_P}\,
  \left[\frac{3\,N_{\rm eff}^{\rm eff}\,\zeta_C^{4}\,\rho^{3}\,e^{-\sigma_A/\alpha}}
  {32\pi\,\alpha\,K_{\rm EM}^{\rm eff}}\right]^{1/4}.$$

- **Path B — RLTM (Tangential Tension closure)**  
  Ward-preserving product kernels; includes the mode-$M{=}1$ remainder $f_{M1}$; weight $e^{-\sigma_B/\alpha}$ with $\sigma_B=\pi/4$.

  $$m_B = \frac{D\,\hbar}{c\,\ell_P}\,
  \sqrt{\frac{\rho\,K_T^{\rm eff}\,e^{-\sigma_B/\alpha}\,(1+f_{M1})}
  {2\,\alpha\,K_{\rm EM}^{\rm raw,eff}}}, \qquad
  f_{M1}=\frac{\alpha}{2}\,e^{-y^{*2}/2}\,y^{*2}.$$


---

## Key results (typical, no tuning)

> Reproduced by the scripts in this repo; digits depend only on the constants set and high-precision arithmetic.

| Quantity | Value / Deviation | Notes |
|---|---:|---|
| Emergent electron mass — **Path A** | $\delta m_A = -2.45046970191$ ppm | Undershoot vs $m_e$; within the $\sim 11$ ppm uncertainty set by $G$. |
| Emergent electron mass — **Path B** | $\delta m_B = -2.45076007903$ ppm | Same envelope. |
| Fine-structure (geometry-only identity) | $\alpha^{-1} \approx 137.0359991769\ldots$ | From $q^2=4\pi\varepsilon_0\,\alpha\,\hbar c$; sub-ppb agreement. |
| Gravitational constant (from $\ell_P$) | $G \approx 6.67426728776\times10^{-11}\ \mathrm{m^3\,kg^{-1}\,s^{-2}}$ | ≈ −4.9 ppm vs CODATA-2022; $G$ limits mass precision. |

---

## Install and run

**Requirements:** Python ≥ 3.10, [`mpmath`](https://pypi.org/project/mpmath/). [`tabulate`](https://pypi.org/project/tabulate/) is optional (pretty tables).

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

Typical outputs are saved in `outputs/` (e.g., `mass_pathA.txt`, `mass_pathB.txt`, `alpha_identity.txt`).

---

## Cite this work

**Pajuhaan, M. (2025). *Electron Mass from Relator Kinematics*. Zenodo.**
[https://doi.org/10.5281/zenodo.17219279](https://doi.org/10.5281/zenodo.17219279)

```bibtex
@misc{Pajuhaan2025ElectronMassRelator,
  author    = {Pajuhaan, Mehrdad},
  title     = {Electron Mass from Relator Kinematics},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.17219279},
  url       = {https://doi.org/10.5281/zenodo.17219279}
}
```

**Manuscript (PDF):** [https://doi.org/10.5281/zenodo.17219278](https://doi.org/10.5281/zenodo.17219278)

---

## Reproducibility

* No fitting or post-matching; the only dimensionful input is $\ell_P=\sqrt{\hbar G/c^3}$.
* All outputs are deterministic for a given commit and constants table.


