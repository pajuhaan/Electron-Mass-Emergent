#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2025 by Mehrdad Pajuhaan
# =============================================================================
#  Relator — Emergent Electron Mass (PAPER-ALIGNED NAMING: Γ_str, K_ring, ζ_ring)
#  --------------------------------------------------------------------
#  This program reproduces the paper's numerics with notation matched to the
#  manuscript output:
#    • Γ_str(y)  (was S(y))
#    • K_ring(y) (was K0(y)),  Λ_ring(y) (was Λ(y)),  ζ_ring(y) (was ζ_geom)
#    • Printed labels now use Γ_str, K_ring, Λ_ring, ζ_ring across reports
#
#  Path A (scalar):  K_EM^{mid,A} → K_EM^{eff,A}
#    - mid  = raw × s_G(κ)                 (Gaussian softening on S¹)
#    - effA = mid × [1 − ½ ρ² ζ_soft],     (ζ_soft from Λ_ind ONLY)
#    - N_eff^{eff,A} = N_eff^{raw,A} × [1 + ρ² 𝓓_phys^S4(α, ζ(α))]
#
#  Path B (vector), Ward-preserving products:
#    K_T^eff = K_T·(1 − ½ρ² ζ_ring)·(1 − ½ρ² ζ(α))
#    K_EM^{raw,eff} = K_EM^{raw}·(1 + ¼ρ² ζ_ring)·(1 + ¼ρ² ζ(α))
#
#  The code prints all named intermediates with paper-style symbols, performs
#  per-term toggle audits, propagates constant uncertainties, and saves a
#  checksumed report for reproducibility.
#
#  IMPORTANT (No double counting):
#    - Path A's local EM kernel never uses ζ(α); it uses ζ_soft from Λ_ind only.
#    - Bulk/vector effects enter Path A only via N_eff^{eff,A} through 𝓓_C and ζ(α).
# =============================================================================

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import mpmath as mp
import os, sys, platform, datetime, hashlib

# ============================= Global precision ==============================
mp.mp.dps = 140  # ppb-class stability

# ========================== Algorithm configuration ==========================
OUT_MODE   = 'exact'   # 'exact'|'series'|'dipole'
OUT_LBASE  = 23        # odd multipole base for OUT extrapolation
GL_NODES   = 512       # Gauss–Legendre nodes on S²

SPEC_M_MAX = 30        # spectral cutoff for L_{2m} in 𝓓_C
SPEC_TOL   = mp.mpf('1e-40')

CURV_SERIES_ON    = True
CURV_SERIES_ORDER = 20

CHI_LADDER_ON = True
SELF_LADDER_ON = True

USE_PARALLEL = True
MAX_WORKERS  = None
MPMATH_DPS_RUNTIME = mp.mp.dps

# ========================= Report file (enabled) =============================
SAVE_REPORT = "checksum.txt"

# ========================= Physical constants (SI) ===========================
PHYS = {
    "c":     {"value": mp.mpf('299792458.0'), "source": "SI (exact)",               "u_rel": mp.mpf('0')},
    "e":     {"value": mp.mpf('1.602176634e-19'), "source": "SI (exact)",            "u_rel": mp.mpf('0')},
    "hbar":  {"value": mp.mpf('1.054571817e-34'), "source": "SI (fixed h → ħ)",      "u_rel": mp.mpf('0')},
    "G":     {"value": mp.mpf('6.67430e-11'),      "source": "CODATA 2022",          "u_rel": mp.mpf('2.2e-5')},
    "alpha_in": {"value": mp.mpf('7.2973525643e-3'), "source": "CODATA 2022 (input)", "u_rel": mp.mpf('1.5e-10')},
    "m_e_exp_MeV": {"value": mp.mpf('0.51099895069'), "source": "CODATA 2022 (reporting only)", "u_rel": mp.mpf('0')}
}

# ================================ Utilities ==================================
def nstr(x, n=28): return mp.nstr(x, n)
def eV_from_MeV(m): return m * mp.mpf('1e6')
def ppm(rel): return mp.mpf('1e6') * rel

def _format_table(headers, rows, title=None):
    """Minimal pretty table; falls back if tabulate is unavailable."""
    try:
        from tabulate import tabulate
        out = []
        if title: out.append(title)
        out.append(tabulate(rows, headers=headers, tablefmt="github", floatfmt=".12g"))
        return "\n".join(out)
    except Exception:
        cols = len(headers)
        widths = [len(str(h)) for h in headers]
        for r in rows:
            for j in range(cols):
                widths[j] = max(widths[j], len(str(r[j])))
        def line(vals):
            return "  ".join(str(vals[j]).ljust(widths[j]) for j in range(cols))
        out = []
        if title: out.append(title)
        out.append(line(headers))
        out.append(line(["-"*w for w in widths]))
        for r in rows: out.append(line(r))
        return "\n".join(out)

# ====================== Paper notation: shared constants =====================
pi      = mp.pi
C0_UNI  = (1/pi) * (mp.mpf('4')/3 + 1/(4*pi**2))  # C_{uni}
C0_GAUSS= mp.mpf('0.5') * (mp.log(2) + mp.euler)  # ½(ln2+γ)
EPSILON = 1/mp.sqrt(pi)                           # ε = 1/√π
ETA0    = 1/pi                                    # η = 1/π
ELL0    = EPSILON*ETA0                            # ℓ = ε η
ITOT    = mp.mpf('1')/6 - mp.mpf('1')/(4*pi**2)   # ∫_0^1 x² sin²(πx) dx

# ============================= Ring integrals (S¹) ===========================
def erfcx(z): return mp.e**(z**2) * mp.erfc(z)  # scaled erfc
def a_of_delta(Delta): return 2*mp.sin(Delta/2)

def J_of(a, y):
    r"""J(a,y)=√π/(2y)·erfcx(a/(2y))  (Gaussian overlap kernel; >0 for y>0)."""
    return mp.sqrt(pi)/(2*y) * erfcx(a/(2*y))

def K_ring_of(y):
    r"""K_ring(y)=(1/π)∫_0^π J(a(Δ),y)dΔ — local EM azimuthal average on the ring."""
    return mp.quad(lambda D: J_of(a_of_delta(D), y), [0, pi]) / pi

def Lambda_ring_of(y):
    r"""Λ_ring(y)=(1/π)∫_0^π cosΔ·J(a(Δ),y)dΔ — local vector overlap on the ring."""
    return mp.quad(lambda D: mp.cos(D) * J_of(a_of_delta(D), y), [0, pi]) / pi

def Gamma_str_of(y, alpha):
    r"""Γ_str(y;α)=2 + ζ_ring(y) + α K_ring(y)/(8π²), where ζ_ring(y)=[K_ring/(2π²)]·Λ_ring."""
    Kring = K_ring_of(y); Lring = Lambda_ring_of(y)
    zeta_ring = (Kring/(2*pi**2))*Lring
    Gamma = (2 + zeta_ring) + alpha * Kring / (8*pi**2)
    return Gamma, zeta_ring, Kring, Lring

def dGamma_dy(y, alpha):
    h = mp.mpf('5e-7')*(1+abs(y))
    Gp,_,_,_ = Gamma_str_of(y+h, alpha); Gm,_,_,_ = Gamma_str_of(y-h, alpha)
    return (Gp - Gm)/(2*h)

def F_shape(y, alpha):
    r"""Stationarity on S¹:  3 Γ'/Γ + (2/y − y) = 0."""
    Gamma,_,_,_ = Gamma_str_of(y, alpha)
    return 3*dGamma_dy(y, alpha)/Gamma + (2/y - y)

def solve_y_star(alpha, y_lo=mp.mpf('0.35'), y_hi=mp.mpf('3.0'), step=mp.mpf('0.02')):
    """Bracket + bisection; fallback to argmin |F| if no sign change is found."""
    yL = y_lo; fL = F_shape(yL, alpha); y = y_lo + step
    while y <= y_hi:
        fR = F_shape(y, alpha)
        if fL*fR < 0:
            a, b = yL, y
            for _ in range(260):
                m  = (a+b)/2
                fm = F_shape(m, alpha)
                if abs(fm) < mp.mpf('1e-32') or abs(b-a) < mp.mpf('1e-32'): return m
                if fL*fm <= 0: b = m
                else: a = m; fL = fm
            return (a+b)/2
        yL, fL = y, fR; y += step
    # Fallback: argmin |F|
    argmin, best = y_lo, abs(F_shape(y_lo, alpha)); yy = y_lo + step
    while yy <= y_hi:
        val = abs(F_shape(yy, alpha))
        if val < best: best, argmin = val, yy
        yy += step
    return argmin

# ================= Spectrum {K, L_{2m}} and Coulomb block 𝓓_C(α) ============
def In_m1(n: int) -> mp.mpf:
    n = mp.mpf(n)
    return (-1)**(n - 1) / (((n - 1) * pi)**2) + (-1)**n / (((n + 1) * pi)**2)

def series_K(tol: mp.mpf = mp.mpf('1e-40')) -> mp.mpf:
    S, n = mp.mpf('0'), 2
    while True:
        t = (2 * In_m1(n))**2 / (n**2 - 1)
        S += t
        if n > 120 and abs(t) < tol: break
        n += 1
    return (2 / pi**2) * S

def _CS_pair(m: int, a: mp.mpf):
    if a == 0: return mp.mpf('1')/(m+1), mp.mpf('0')
    C = mp.sin(a)/a; S = (1 - mp.cos(a))/a
    if m == 0: return C, S
    for k in range(1, m+1):
        C, S = mp.sin(a)/a - (k/a)*S, (1 - mp.cos(a))/a + (k/a)*C
    return C, S

def I_nm(n: int, m2: int) -> mp.mpf:
    C1, _ = _CS_pair(m2, (n - 1)*pi); C2, _ = _CS_pair(m2, (n + 1)*pi)
    return mp.mpf('0.5') * (C1 - C2)

def L_2m(m: int, tol: mp.mpf = mp.mpf('1e-40'), nmax: int = 1200) -> mp.mpf:
    S = mp.mpf('0')
    for n in range(2, nmax+1):
        t = (2*I_nm(n, 2*m))**2 / (n**2 - 1)
        S += t
        if n > 80 and abs(t) < tol: break
    return (2/pi**2) * S

def precompute_spectrum(M: int = SPEC_M_MAX, tol: mp.mpf = SPEC_TOL):
    K = series_K(tol=tol)
    L_list = []
    if M >= 2:
        for m in range(2, M+1):
            L_list.append((m, L_2m(m, tol=tol)))
    return K, L_list

def DC_of_alpha_fixedK(alpha: mp.mpf, K: mp.mpf, L_list, M: int = SPEC_M_MAX) -> mp.mpf:
    r"""𝓓_C(α) = (α/π)[ √(1−ξ) − (ξ/2)K − Σ_{m≥2} (ξ/2)^m L_{2m} ], with ξ=2 C0^{uni} α."""
    a  = mp.mpf(alpha)
    xi = 2 * C0_UNI * a
    D  = (a/pi)*mp.sqrt(1 - xi) - (a/pi)*(xi/2)*K
    if M >= 2:
        for (m, Lm) in L_list:
            D -= (a/pi) * ((xi/2)**m) * Lm
    return D

# ================================ OUT block =================================
_GL_CACHE = {}
def _gauss_legendre(n: int):
    if n in _GL_CACHE: return _GL_CACHE[n]
    xs, ws = [], []
    tol = mp.mpf(10)**(-(mp.mp.dps - 8))
    for k in range(1, n+1):
        x = mp.cos(pi * (k - mp.mpf('0.25')) / (n + mp.mpf('0.5')))
        for _ in range(25):
            Pn  = mp.legendre(n, x)
            dPn = n/(1 - x**2) * (mp.legendre(n-1, x) - x*Pn)
            dx  = -Pn/dPn
            x  += dx
            if abs(dx) < tol: break
        w = 2/((1 - x**2)*(dPn**2))
        xs.append(x); ws.append(w)
    _GL_CACHE[n] = (xs, ws); return xs, ws

def _B_rho(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    rc = mp.sqrt((1 + rho)**2 + z**2)
    k2 = 4*rho/((1 + rho)**2 + z**2)
    if k2 <= 0: return mp.mpf('0')
    if k2 >= 1: k2 = mp.mpf('1') - mp.mpf('1e-30')
    K_ = mp.ellipk(k2); E_ = mp.ellipe(k2)
    denom = (1 - rho)**2 + z**2
    if rho == 0: return mp.mpf('0')
    return (z / (2*pi*rho*rc)) * (-K_ + ((1 + rho**2 + z**2)/denom) * E_)

def _B_z(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    rc = mp.sqrt((1 + rho)**2 + z**2)
    k2 = 4*rho/((1 + rho)**2 + z**2)
    if k2 <= 0: return mp.mpf('1')/(2*rc**3)
    if k2 >= 1: k2 = mp.mpf('1') - mp.mpf('1e-30')
    K_ = mp.ellipk(k2); E_ = mp.ellipe(k2)
    denom = (1 - rho)**2 + z**2
    return (1/(2*pi*rc)) * (K_ + ((1 - rho**2 - z**2)/denom) * E_)

def _dPdx_leg(l: int, x: mp.mpf) -> mp.mpf:
    if l == 0: return mp.mpf('0')
    Pl  = mp.legendre(l, x); Pl1 = mp.legendre(l-1, x)
    return (l/(1 - x**2)) * (Pl1 - x*Pl)

def _Btheta_on_sphere_x(x: mp.mpf, rstar: mp.mpf) -> mp.mpf:
    """Tangential field on the unit sphere (x=cosθ), evaluated at R=r*."""
    s = mp.sqrt(1 - x**2); rho = rstar*s; z = rstar*x
    return _B_rho(rho, z) * x - _B_z(rho, z) * s

def Lambda_OUT_exact(eta: mp.mpf, lmax: int = OUT_LBASE, gl_nodes: int = GL_NODES) -> mp.mpf:
    """Exact OUT energy via multipoles up to lmax; then Λ_out = −2 U_out."""
    if OUT_MODE == 'dipole':
        return - (pi/6) * (eta**3)
    rstar = 1/eta
    xs, ws = _gauss_legendre(gl_nodes)
    Uout = mp.mpf('0')
    Lm = lmax if (lmax % 2 == 1) else (lmax - 1)
    for l in range(1, Lm+1, 2):
        Il = l*(l+1)*2/(2*l + 1)
        s  = mp.mpf('0')
        for x, w in zip(xs, ws):
            s += _Btheta_on_sphere_x(x, rstar) * (-(1 - x**2) * _dPdx_leg(l, x)) * w
        a_l  = - (rstar**(l + 2)) * s / Il
        Uout += ((l + 1)/(2*l + 1)) * (a_l**2) * (rstar**(-(2*l + 1)))
    Uout *= (2*pi)
    return -2 * Uout

def Lambda_OUT_series(eta: mp.mpf, lmax: int = OUT_LBASE) -> mp.mpf:
    s  = mp.mpf('0')
    Lm = lmax if (lmax % 2 == 1) else (lmax - 1)
    for n in range(1, Lm+1, 2):
        s += pi/((n+1)*(2*n + 1)) * (eta**(2*n + 1))
    return -s

def Lambda_OUT(eta: mp.mpf, lmax: int = OUT_LBASE, mode: str = OUT_MODE) -> mp.mpf:
    return Lambda_OUT_series(eta, lmax) if mode == 'series' else \
           Lambda_OUT_exact(eta, lmax, gl_nodes=GL_NODES)

def Lambda_OUT_extrapolated(eta: mp.mpf, lbase: int = OUT_LBASE, mode: str = OUT_MODE) -> mp.mpf:
    """Richardson-like extrapolation using three odd multipoles L,L+2,L+4."""
    S1 = Lambda_OUT(eta, lmax=lbase,     mode=mode)
    S2 = Lambda_OUT(eta, lmax=lbase + 2, mode=mode)
    S3 = Lambda_OUT(eta, lmax=lbase + 4, mode=mode)
    denom = (S3 - 2*S2 + S1)
    if denom == 0 or abs(denom) < mp.mpf('1e-30') * max(1, abs(S3)): return S3
    Sout = S1 - (S2 - S1)**2 / denom
    if abs(Sout) > 10 * max(abs(S1), abs(S2), abs(S3)): return S3
    return Sout

# ============================= Λ-chain & ζ(α) ================================
def P_IR_chi(ell: mp.mpf) -> mp.mpf:
    """IR acceptance P_IR^(χ)(ℓ) with exact S¹ weight (no small-κ expansions)."""
    w = lambda x: (x**2) * (mp.sin(pi*x)**2)
    num = mp.quad(lambda x: w(x) * (1 - mp.mpf('1')/3 * ((1 - x)**2/(((1 - x)**2) + ell**2))) *
                            mp.e**(-((1 - x)/ell)**2), [0, 1])
    return num / ITOT

def Lambda_ind(eps: mp.mpf) -> mp.mpf:
    """Λ_ind = ln(8/ε) − 2  (local inductive baseline)."""
    return mp.log(8/eps) - 2

def DeltaLambda_UVIR(P_ir: mp.mpf) -> mp.mpf:
    """ΔΛ_{UV→IR} = C0^{Gauss} · P_IR^(χ)."""
    return C0_GAUSS * P_ir

def curvature_series_eta(eta: mp.mpf, order: int = CURV_SERIES_ORDER) -> mp.mpf:
    if not CURV_SERIES_ON: return eta**2 / 6
    if order < 2: return mp.mpf('0')
    s = mp.mpf('0'); max_k = int(order // 2)
    for k in range(1, max_k + 1):
        s += (eta**(2*k)) / mp.factorial(2*k + 1)
    return s

def gamma_geom(eta: mp.mpf) -> mp.mpf:
    """Geometric gain γ_geom(η) on S² (curvature-series)."""
    return mp.mpf('0.5') * (1 + curvature_series_eta(eta, order=CURV_SERIES_ORDER))

def gamma_eff(eta: mp.mpf, K: mp.mpf, DC_lock: mp.mpf, P_ir: mp.mpf) -> mp.mpf:
    """Effective gain γ_eff = γ_geom + γ_map(𝓓_C)."""
    return gamma_geom(eta) + (K / (2 * DC_lock)) * C0_GAUSS * P_ir

def deltaLambda_chi_ladder_extra(eta: mp.mpf, K: mp.mpf, DC_lock: mp.mpf,
                                 P_ir: mp.mpf, dLambda_OUT: mp.mpf) -> mp.mpf:
    """χ-ladder closed resummation (extra term)."""
    k = mp.sinh(eta)/eta - 1
    x = (K / (2 * DC_lock)) * C0_GAUSS * P_ir
    return (-k * x**2 / (1 + k * x)) * P_ir * dLambda_OUT

def deltaLambda_self_ladder(eta: mp.mpf, K: mp.mpf, P_ir: mp.mpf,
                            Lambda_eff: mp.mpf, alpha: mp.mpf) -> mp.mpf:
    """Self-ladder (Dyson-like) closed correction."""
    ggeom = gamma_geom(eta)
    k = mp.sinh(eta)/eta - 1
    ep = (alpha / pi) * (K / (2 * pi**2)) * P_ir * Lambda_eff
    return - ggeom * P_ir * Lambda_eff * (ep / (1 + k * ep))

def build_Lambda_eff(alpha: mp.mpf, K: mp.mpf, L_list):
    """Full Λ_eff(α) with OUT extrapolation, ladders, and DC from spectrum."""
    DC = DC_of_alpha_fixedK(alpha, K, L_list, M=SPEC_M_MAX)
    P_ir  = P_IR_chi(ELL0)
    Lin   = Lambda_ind(EPSILON)
    dUV   = DeltaLambda_UVIR(P_ir)
    dOUT  = Lambda_OUT_extrapolated(ETA0, lbase=OUT_LBASE, mode=OUT_MODE)
    base  = Lin + dUV + dOUT
    g_eff = gamma_eff(ETA0, K, DC, P_ir)
    dSYNC = g_eff * P_ir * dOUT
    if CHI_LADDER_ON:
        dSYNC += deltaLambda_chi_ladder_extra(ETA0, K, DC, P_ir, dOUT)
    Lam_eff = base + dSYNC
    if SELF_LADDER_ON:
        dSYNC += deltaLambda_self_ladder(ETA0, K, P_ir, Lam_eff, alpha)
        Lam_eff = base + dSYNC
    return {
        "P_ir": P_ir, "Lambda_ind": Lin, "Delta_UVIR": dUV,
        "Delta_out": dOUT, "Delta_sync": dSYNC, "Lambda_eff": Lam_eff,
        "DC": DC
    }

# ========================= Ring-level kernels (Path A/B) =====================
def K_EM_raw_closed(x):
    """K_EM^{raw}(x) = (2x/π) K(m=-4x²)."""
    m = -4*(x**2)
    return (2*x/pi) * mp.ellipk(m)

def kappa_tension(x): return (x**2)/(1+x**2)

def softening_factor(kappa):
    """Gaussian EM-softening; stable as κ→0 via scaled erfc."""
    if kappa <= mp.mpf('0'):
        return mp.mpf('1')
    if kappa < mp.mpf('1e-20'):
        return mp.mpf('1') - kappa/2 + (3*(kappa**2))/4
    return mp.sqrt(mp.pi)/mp.sqrt(kappa) * erfcx(1/mp.sqrt(kappa))

def K_EM_mid_A_closed(x):
    """K_EM^{mid,A}(x) = K_EM^{raw}(x) × s_G(κ(x))   (purely geometric)."""
    return K_EM_raw_closed(x) * softening_factor(kappa_tension(x))

def N_eff_raw_A_closed(x):
    """N_eff^{raw,A}(x) = 2(1+κ)(1 − ½ κ²)."""
    κ = kappa_tension(x)
    return 2*(1+κ)*(1 - mp.mpf('0.5')*κ**2)

def K_T_closed(x): return 2*kappa_tension(x)

# ============================= Geometric RG locks ============================
def beta0_from_DS(D_over_Gamma): return (6/pi)*(D_over_Gamma)   # D/Γ=4/3 → β0=8/π
def sigmas_from_beta0(beta0): return 4/beta0, 2/beta0          # (σ_A, σ_B)

# =============================== Mass formulae ===============================
def mass_path_A(ELL_P, D, rho, KEM_eff_A, N_eff_eff_A, ZETA_C, sigma_A, alpha, C, HBAR, ECHG):
    r"""Path A (scalar):
       m_A = (Dħ/cℓ_P) [ 3 N_eff^{eff,A} ζ_C^4 ρ^3 e^{-σ_A/α} / (32π α K_EM^{eff,A}) ]^{1/4}.
    """
    lock = mp.e**(-sigma_A/alpha)
    num  = 3*N_eff_eff_A*(ZETA_C**4)*(rho**3)*lock
    den  = 32*pi*alpha*KEM_eff_A
    F    = (num/den)**(mp.mpf('1')/4)
    mkg  = (D*HBAR)/(C*ELL_P) * F
    return (mkg*C**2)/ECHG/1e6  # MeV

def mass_path_B(ELL_P, D, rho, K_EM_raw_eff_B, K_T_eff_B, sigma_B, y_star, alpha, C, HBAR, ECHG, use_fM1=True):
    r"""Path B (vector):
       m_B = (Dħ/cℓ_P) √[ ρ K_T^{eff} e^{-σ_B/α} (1+f_{M1}) / (2 α K_EM^{raw,eff}) ].
       f_{M1} = (α/2) e^{−y*²/2} y*²   (ring-local M1 dressing).
    """
    lock = mp.e**(-sigma_B/alpha)
    fM1  = (alpha*mp.e**(-(y_star**2)/2)*(y_star**2)/2) if use_fM1 else mp.mpf('0')
    F    = mp.sqrt( (rho*K_T_eff_B*lock*(1+fM1)) / (2*alpha*K_EM_raw_eff_B) )
    mkg  = (D*HBAR)/(C*ELL_P) * F
    return (mkg*C**2)/ECHG/1e6  # MeV

# ====================== α from equality (diagnostic only) ====================
def alpha_eq_closed(KEM_eff_A, N_eff_eff_A, rho, K_T, K_T_eff_B, K_EM_raw, K_EM_raw_eff_B, y_star, alpha):
    r"""Closed diagnostic for α̂ from A=B (micro factors explicit):
       α̂ = (8/(3π^3)) (K_EM^{eff,A}/(N_eff^{eff,A} ρ)) (K_T/K_EM^{raw})^2
            [ (K_T^{eff}/K_T)/(K_EM^{raw,eff}/K_EM^{raw}) ]^2 (1+f_{M1})^2.
    """
    fM1  = alpha*mp.e**(-(y_star**2)/2)*(y_star**2)/2
    base  = (K_T / K_EM_raw)**2
    micro = ((K_T_eff_B / K_T) / (K_EM_raw_eff_B / K_EM_raw))**2
    return (mp.mpf('8')/(3*pi**3)) * (KEM_eff_A/(N_eff_eff_A*rho)) * base * micro * (1+fM1)**2

# ============================ Parallel helpers ===============================
def _mp_init(dps): mp.mp.dps = dps
def _task_Kring(y, dps):   _mp_init(dps); return K_ring_of(y)
def _task_Lring(y, dps):   _mp_init(dps); return Lambda_ring_of(y)
def _task_Kraw(x, dps):    _mp_init(dps); return K_EM_raw_closed(x)
def _task_KmidA(x, dps):   _mp_init(dps); return K_EM_mid_A_closed(x)

def eval_weights_parallel(y_star, x_rel, dps=MPMATH_DPS_RUNTIME, max_workers=None):
    """Concurrent evaluation of heavy S¹ weights."""
    tasks = {
        "Kring": (_task_Kring, (y_star, dps)),
        "Lring": (_task_Lring, (y_star, dps)),
        "Kraw":  (_task_Kraw,  (x_rel,  dps)),
        "KmidA": (_task_KmidA, (x_rel,  dps)),
    }
    out = {}
    workers = max_workers or os.cpu_count() or 4
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = { ex.submit(func, *args): key for key,(func,args) in tasks.items() }
        for f in as_completed(futs):
            out[futs[f]] = f.result()
    return out

# =============================== Core pipeline ===============================
def compute_all(alpha_in):
    """One complete evaluation at the given α_in (no tuning)."""
    # SI constants
    HBAR  = PHYS["hbar"]["value"]
    C     = PHYS["c"]["value"]
    ECHG  = PHYS["e"]["value"]
    G     = PHYS["G"]["value"]

    # Planck scales
    ELL_P = mp.sqrt(HBAR*G/C**3)
    t_P   = mp.sqrt(HBAR*G/C**5)

    # Spectrum and Coulomb block
    K_spec, L_list = precompute_spectrum(M=SPEC_M_MAX, tol=SPEC_TOL)

    # Λ-chain at α_in (exact + L-extrap), per Alpha code
    LAMBDA = build_Lambda_eff(alpha_in, K_spec, L_list)
    zeta_alpha = (K_spec/(2*pi**2)) * LAMBDA["Lambda_eff"]

    # Ring reduction + shape
    y_star = solve_y_star(alpha_in)
    Gamma_star, zeta_ring, K_ring_y, Lambda_ring_y = Gamma_str_of(y_star, alpha_in)
    D      = (mp.mpf('4')/3)*Gamma_star
    rho    = 1/D
    x_rel  = rho / y_star

    # Heavy S¹ weights (parallel)
    if USE_PARALLEL:
        W = eval_weights_parallel(y_star, x_rel, dps=mp.mp.dps, max_workers=MAX_WORKERS)
        K_ring_y, Lambda_ring_y = W["Kring"], W["Lring"]
        zeta_ring = (K_ring_y/(2*pi**2)) * Lambda_ring_y
        K_EM_raw  = W["Kraw"]
        K_EM_mid_A= W["KmidA"]
    else:
        K_EM_raw  = K_EM_raw_closed(x_rel)
        K_EM_mid_A= K_EM_mid_A_closed(x_rel)

    # Path A: local EM (mid→eff) and scalar completion (raw→eff)
    zeta_soft   = (K_spec/(2*pi**2)) * LAMBDA["Lambda_ind"]        # Λ_ind ONLY
    K_EM_eff_A  = K_EM_mid_A * (1 - mp.mpf('0.5')*rho**2 * zeta_soft)
    N_eff_raw_A = N_eff_raw_A_closed(x_rel)
    # Scalar completion: 𝓓_phys^S4(α, ζ(α))
    D_C         = LAMBDA["DC"]
    D_phys_S4   = D_C - (alpha_in/pi)*zeta_alpha + ((alpha_in*zeta_alpha)/pi)**2 / (4*D_C)
    N_eff_eff_A = N_eff_raw_A * (1 + rho**2 * D_phys_S4)
    ZETA_C      = pi

    # Path B: Ward-preserving vector completions
    K_T     = K_T_closed(x_rel)
    K_T_eff = K_T * (1 - mp.mpf('0.5') * rho**2 * zeta_ring) * (1 - mp.mpf('0.5') * rho**2 * zeta_alpha)
    K_EM_raw_eff_B = K_EM_raw * (1 + mp.mpf('0.25') * rho**2 * zeta_ring) * (1 + mp.mpf('0.25') * rho**2 * zeta_alpha)

    # RG locks from geometry (β0≈8/π)
    beta0            = beta0_from_DS(D/Gamma_star)
    sigma_A, sigma_B = sigmas_from_beta0(beta0)

    # Masses
    mA = mass_path_A(ELL_P, D, rho, K_EM_eff_A, N_eff_eff_A, ZETA_C, sigma_A, alpha_in, C, HBAR, ECHG)
    mB = mass_path_B(ELL_P, D, rho, K_EM_raw_eff_B, K_T_eff,    sigma_B, y_star,   alpha_in, C, HBAR, ECHG, use_fM1=True)

    # Diagnostic: α̂ from equality (no scales)
    alpha_hat = alpha_eq_closed(K_EM_eff_A, N_eff_eff_A, rho, K_T, K_T_eff, K_EM_raw, K_EM_raw_eff_B, y_star, alpha_in)

    # Reduced Compton radii
    rbarA = (HBAR*C) / (mA*1e6*ECHG)
    rbarB = (HBAR*C) / (mB*1e6*ECHG)

    return {
        # shared
        "alpha_in": alpha_in, "ELL_P": ELL_P, "t_P": t_P,
        "y_star": y_star, "Gamma": Gamma_star, "D": D, "rho": rho, "beta0": beta0,
        "x_rel": x_rel, "K_ring": K_ring_y, "Lambda_ring": Lambda_ring_y, "zeta_ring": zeta_ring,
        # spectrum & Coulomb
        "K_spec": K_spec, "L_list": L_list, "D_C": D_C, "D_phys_S4": D_phys_S4,
        # Λ-chain & ζ(α)
        "P_ir": LAMBDA["P_ir"], "Lambda_ind": LAMBDA["Lambda_ind"],
        "Delta_UVIR": LAMBDA["Delta_UVIR"], "Delta_out": LAMBDA["Delta_out"],
        "Delta_sync": LAMBDA["Delta_sync"], "Lambda_eff": LAMBDA["Lambda_eff"],
        "zeta_alpha": zeta_alpha, "zeta_soft": zeta_soft,
        # kernels (A/B)
        "K_EM_raw": K_EM_raw,
        "K_EM_mid_A": K_EM_mid_A, "K_EM_eff_A": K_EM_eff_A,
        "N_eff_raw_A": N_eff_raw_A, "N_eff_eff_A": N_eff_eff_A,
        "K_T": K_T, "K_T_eff": K_T_eff, "K_EM_raw_eff_B": K_EM_raw_eff_B,
        # RG & masses
        "sigma_A": sigma_A, "sigma_B": sigma_B, "mA": mA, "mB": mB,
        "rbarA": rbarA, "rbarB": rbarB,
        # equality alpha
        "alpha_hat": alpha_hat
    }

# ============================== Reporting layer ==============================
def reduced_compton_radius_from_MeV(m_MeV):
    HBAR  = PHYS["hbar"]["value"]; C = PHYS["c"]["value"]; E = PHYS["e"]["value"]
    return (HBAR*C)/(m_MeV*1e6*E)

def effect_line(label, m_ref, m_alt):
    dMeV = m_ref - m_alt
    d_eV = eV_from_MeV(dMeV)
    rel  = dMeV / m_ref
    return [label, nstr(d_eV,14), nstr(100*rel,12), nstr(ppm(rel),12)]

def _analytic_mass_exponents():
    """
    Exact scaling of m with constants (holding α fixed):
      m ∝ ( √ħ * c^(5/2) ) / ( √G * e ) × F(α, geometry)
    ⇒ log-derivatives:
      d ln m / d ln ħ = +1/2
      d ln m / d ln c = +5/2
      d ln m / d ln G = -1/2
      d ln m / d ln e = -1
    """
    return {"hbar": mp.mpf('0.5'), "c": mp.mpf('2.5'), "G": -mp.mpf('0.5'), "e": -mp.mpf('1.0')}

def _alpha_log_sensitivity(alpha0, rel_step=mp.mpf('1e-6')):
    a0 = alpha0; dr = mp.mpf(rel_step)
    Rp = compute_all(a0*(1+dr)); Rm = compute_all(a0*(1-dr))
    sA = (mp.log(Rp["mA"]) - mp.log(Rm["mA"])) / (2*dr)
    sB = (mp.log(Rp["mB"]) - mp.log(Rm["mB"])) / (2*dr)
    return sA, sB

def main():
    alpha_in = PHYS["alpha_in"]["value"]
    R = compute_all(alpha_in)

    blocks = []
    def show(block, add_blank=True):
        print(block)
        blocks.append(block)
        if add_blank:
            print()
            blocks.append("")

    # 0) Run provenance & configuration
    prov_rows = [
        ["Date (UTC)", datetime.datetime.utcnow().isoformat() + "Z"],
        ["Python", sys.version.split()[0]],
        ["Platform", platform.platform()],
        ["mpmath", getattr(mp, "__version__", "unknown")],
        ["mp.dps", str(mp.mp.dps)],
        ["OUT_MODE", OUT_MODE], ["OUT_LBASE", OUT_LBASE], ["GL_NODES", GL_NODES],
        ["SPEC_M_MAX", SPEC_M_MAX], ["CURV_SERIES_ON", CURV_SERIES_ON],
        ["CURV_SERIES_ORDER", CURV_SERIES_ORDER], ["CHI_LADDER_ON", CHI_LADDER_ON],
        ["SELF_LADDER_ON", SELF_LADDER_ON], ["USE_PARALLEL", USE_PARALLEL],
        ["MAX_WORKERS", MAX_WORKERS],
    ]
    block = _format_table(["Key", "Value"], prov_rows, title="=== RUN PROVENANCE & CONFIG ===")
    show(block)

    # 1) Shared geometry & constants
    shared_rows = [
        ["α (input)",                nstr(R["alpha_in"], 20), PHYS["alpha_in"]["source"]],
        ["ℓ_P [m] = sqrt(ħG/c^3)",  nstr(R["ELL_P"], 26), ""],
        ["t_P [s] = sqrt(ħG/c^5)",  nstr(R["t_P"], 26), ""],
        ["y* (shape root)",         nstr(R["y_star"], 24), "3 Γ'/Γ + (2/y − y)=0"],
        ["Γ_str(y*)",               nstr(R["Gamma"], 24), "2+ζ_ring+αK_ring/8π²"],
        ["D = 4Γ_str/3",            nstr(R["D"], 24), ""],
        ["ρ = 1/D",                 nstr(R["rho"], 24), ""],
        ["x = ρ/y*",                nstr(R["x_rel"], 24), ""],
        ["K_ring(y*)",              nstr(R["K_ring"], 24), ""],
        ["Λ_ring(y*)",              nstr(R["Lambda_ring"], 24), ""],
        ["ζ_ring(y*)",              nstr(R["zeta_ring"], 24), "[K_ring/(2π²)]Λ_ring"],
        ["β0 (from D/Γ_str)",       nstr(R["beta0"], 24), "≈ 8/π"],
    ]
    block = _format_table(["Variable / Formula", "Value", "Note"], shared_rows, title="=== SHARED GEOMETRY & CONSTANTS ===")
    show(block)

    # 2) Spectrum & 𝓓_C(α)
    spec_rows = [
        ["K (spectral)",  nstr(R["K_spec"], 26), ""],
        ["𝓓_C(α)",        nstr(R["D_C"], 26),    ""],
        ["𝓓_phys^S4",     nstr(R["D_phys_S4"], 26), "D_C − (α/π)ζ + (αζ/π)^2/(4 D_C)"],
    ]
    block = _format_table(["Quantity", "Value", "Comment"], spec_rows, title="=== SPECTRUM & 𝓓_C(α) ===")
    show(block)

    # 3) Λ-chain (exact + L-extrapolated)
    lam_rows = [
        ["P_IR^(χ)(ℓ0)",           nstr(R["P_ir"], 26), ""],
        ["Λ_ind",                   nstr(R["Lambda_ind"], 26), ""],
        ["ΔΛ_UV→IR",                nstr(R["Delta_UVIR"], 26), ""],
        ["ΔΛ_out (extrapolated)",   nstr(R["Delta_out"], 26), ""],
        ["ΔΛ_sync (γ_eff + ladders)", nstr(R["Delta_sync"], 26), ""],
        ["Λ_eff",                   nstr(R["Lambda_eff"], 26), ""],
        ["ζ(α) = (K/2π²)Λ_eff",    nstr(R["zeta_alpha"], 24), ""],
        ["ζ_soft = (K/2π²)Λ_ind",  nstr(R["zeta_soft"], 24), "Path A local only"],
    ]
    block = _format_table(["Λ-chain term", "Value", "Note"], lam_rows, title="=== Λ-CHAIN (Alpha-Lock) ===")
    show(block)

    # 4) Path A — scalar (mid→eff and raw→eff)
    pathA_rows = [
        ["K_EM^{raw}",         nstr(R["K_EM_raw"], 26), ""],
        ["K_EM^{mid,A}",       nstr(R["K_EM_mid_A"], 26), "raw × s_G(κ)"],
        ["K_EM^{eff,A}",       nstr(R["K_EM_eff_A"], 26), "mid × (1 − ½ρ² ζ_soft)"],
        ["N_eff^{raw,A}",      nstr(R["N_eff_raw_A"], 26), "2(1+κ)(1−½κ²)"],
        ["N_eff^{eff,A}",      nstr(R["N_eff_eff_A"], 26), "× (1 + ρ² 𝓓_phys^S4)"],
        ["σ_A",                nstr(R["sigma_A"], 24), ""],
        ["m_A [MeV]",          nstr(R["mA"], 18), ""],
    ]
    block = _format_table(["Path A Quantity", "Value", "Comment"], pathA_rows, title="=== PATH A — SCALAR CHANNEL ===")
    show(block)

    # 5) Path B — vector (Ward product)
    pathB_rows = [
        ["K_T",              nstr(R["K_T"], 24), "2κ"],
        ["K_T^{eff}",        nstr(R["K_T_eff"], 24), "·(1−½ρ²ζ_ring)(1−½ρ²ζ(α))"],
        ["K_EM^{raw,eff}",   nstr(R["K_EM_raw_eff_B"], 24), "raw·(1+¼ρ²ζ_ring)(1+¼ρ²ζ(α))"],
        ["σ_B",              nstr(R["sigma_B"], 24), ""],
        ["m_B [MeV]",        nstr(R["mB"], 18), ""],
    ]
    block = _format_table(["Path B Quantity", "Value", "Comment"], pathB_rows, title="=== PATH B — VECTOR CHANNEL ===")
    show(block)

    # 6) Errors vs experiment (eV & ppm)
    m_ref = PHYS["m_e_exp_MeV"]["value"]
    def err_line(tag, m):
        d = m - m_ref; rel = d/m_ref
        return [tag, nstr(eV_from_MeV(d), 14), nstr(100*rel, 12), nstr(ppm(rel), 12)]
    err_rows = [
        err_line("Path A — final", R["mA"]),
        err_line("Path B — final", R["mB"]),
    ]
    block = _format_table(["Case", "Δ vs exp [eV]", "Δ [%]", "Δ [ppm]"], err_rows, title="=== ERRORS vs EXPERIMENT ===")
    show(block)

    # 7) Path-to-path gap (energy & length)
    rA = R["rbarA"]; rB = R["rbarB"]
    gap_rows = [
        ["Mass gap (A−B) [eV]", nstr(eV_from_MeV(R["mA"] - R["mB"]), 14), ""],
        ["Radius gap (A−B) [m]", nstr(rA - rB, 26), ""],
        ["r̄_C(A), r̄_C(B) [m]", f"{nstr(rA,26)} ; {nstr(rB,26)}", ""],
    ]
    block = _format_table(["Quantity", "Value", "Note"], gap_rows, title="=== PATH-TO-PATH GAP ===")
    show(block)

    # 8) α̂ from equality (diagnostic)
    block = _format_table(
        ["Quantity", "Value", "Note"],
        [["α̂ (from equality)", nstr(R["alpha_hat"], 22), "diagnostic; compare vs input α"]],
        title="=== α FROM PATH EQUALITY (DIAGNOSTIC) ==="
    )
    show(block)

    # 9) Per-term effects (toggles) — reviewer audit
    # Path A:
    K_final_gauss_off = R["K_EM_raw"] * (1 - mp.mpf('0.5')*R["rho"]**2 * R["zeta_soft"])
    mA_noGaussian = mass_path_A(R["ELL_P"], R["D"], R["rho"], K_final_gauss_off, R["N_eff_eff_A"], pi, R["sigma_A"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"])
    mA_noScalar   = mass_path_A(R["ELL_P"], R["D"], R["rho"], R["K_EM_eff_A"], R["N_eff_raw_A"],     pi, R["sigma_A"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"])
    effA_rows = [
        effect_line("A: Gaussian softening off (mid→raw)", R["mA"], mA_noGaussian),
        effect_line("A: scalar DC off (N_eff^{raw,A})",    R["mA"], mA_noScalar),
        effect_line("A: RG locking off (σ_A=0)",           R["mA"], mass_path_A(R["ELL_P"], R["D"], R["rho"], R["K_EM_eff_A"], R["N_eff_eff_A"], pi, mp.mpf('0'), R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"])),
    ]
    block = _format_table(["Toggle (Path A term)", "Δ [eV]", "Δ [%]", "Δ [ppm]"], effA_rows, title="=== PER-TERM EFFECTS — PATH A ===")
    show(block)

    # Path B:
    mB_noF     = mass_path_B(R["ELL_P"], R["D"], R["rho"], R["K_EM_raw_eff_B"], R["K_T_eff"], R["sigma_B"], R["y_star"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"], use_fM1=False)
    mB_noVecT  = mass_path_B(R["ELL_P"], R["D"], R["rho"], R["K_EM_raw_eff_B"], R["K_T"],     R["sigma_B"], R["y_star"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"], use_fM1=True)
    mB_noVecE  = mass_path_B(R["ELL_P"], R["D"], R["rho"], R["K_EM_raw"],       R["K_T_eff"], R["sigma_B"], R["y_star"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"], use_fM1=True)
    mB_bothOff = mass_path_B(R["ELL_P"], R["D"], R["rho"], R["K_EM_raw"],       R["K_T"],     R["sigma_B"], R["y_star"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"], use_fM1=True)
    effB_rows = [
        effect_line("B: f_M1 off",                    R["mB"], mB_noF),
        effect_line("B: vector on K_T off",           R["mB"], mB_noVecT),
        effect_line("B: vector on K_EM^raw off",      R["mB"], mB_noVecE),
        effect_line("B: both vector completions off", R["mB"], mB_bothOff),
        effect_line("B: RG locking off (σ_B=0)",      R["mB"], mass_path_B(R["ELL_P"], R["D"], R["rho"], R["K_EM_raw_eff_B"], R["K_T_eff"], mp.mpf('0'), R["y_star"], R["alpha_in"], PHYS["c"]["value"], PHYS["hbar"]["value"], PHYS["e"]["value"], use_fM1=True)),
    ]
    block = _format_table(["Toggle (Path B term)", "Δ [eV]", "Δ [%]", "Δ [ppm]"], effB_rows, title="=== PER-TERM EFFECTS — PATH B ===")
    show(block)

    # 10) Uncertainty propagation (1σ)
    exps = _analytic_mass_exponents()
    def lines_from_analytic(const_key, urel, mA, mB):
        if urel == 0:
            dA_eV = mp.mpf('0'); dB_eV = mp.mpf('0'); pA = mp.mpf('0'); pB = mp.mpf('0')
        else:
            s = abs(exps[const_key])
            dA_eV = eV_from_MeV(s * urel * mA)
            dB_eV = eV_from_MeV(s * urel * mB)
            pA    = s * urel * mp.mpf('1e6')
            pB    = s * urel * mp.mpf('1e6')
        return [const_key, nstr(urel,10), nstr(dA_eV,12), nstr(pA,10), nstr(dB_eV,12), nstr(pB,10), "analytic"]

    sA, sB = _alpha_log_sensitivity(alpha_in, rel_step=mp.mpf('1e-6'))
    dA_alpha = abs(sA) * PHYS["alpha_in"]["u_rel"] * R["mA"]
    dB_alpha = abs(sB) * PHYS["alpha_in"]["u_rel"] * R["mB"]
    row_alpha = ["alpha_in", nstr(PHYS["alpha_in"]["u_rel"],10),
                 nstr(eV_from_MeV(dA_alpha),12), nstr(abs(sA)*PHYS["alpha_in"]["u_rel"]*mp.mpf('1e6'),10),
                 nstr(eV_from_MeV(dB_alpha),12), nstr(abs(sB)*PHYS["alpha_in"]["u_rel"]*mp.mpf('1e6'),10),
                 "sensitivity"]

    unc_rows = [
        lines_from_analytic("hbar",  PHYS["hbar"]["u_rel"],  R["mA"], R["mB"]),
        lines_from_analytic("c",     PHYS["c"]["u_rel"],     R["mA"], R["mB"]),
        lines_from_analytic("e",     PHYS["e"]["u_rel"],     R["mA"], R["mB"]),
        lines_from_analytic("G",     PHYS["G"]["u_rel"],     R["mA"], R["mB"]),
        row_alpha,
    ]
    block = _format_table(
        ["Constant", "rel σ", "Δm_A [eV]", "Δm_A [ppm]", "Δm_B [eV]", "Δm_B [ppm]", "Method"],
        unc_rows,
        title="=== UNCERTAINTY PROPAGATION (1σ EFFECT ON MASS) ==="
    )
    show(block)

    # 11) Constants summary (values, sources, uncertainties)
    const_rows = []
    for key in ["c","e","hbar","G","alpha_in","m_e_exp_MeV"]:
        d = PHYS[key]
        tag = key if key != "alpha_in" else "alpha_in (input)"
        if key == "m_e_exp_MeV": tag = "m_e_exp_MeV (reporting only)"
        const_rows.append([tag, nstr(d["value"], 20), d["source"], nstr(d["u_rel"], 12)])
    block = _format_table(["Constant", "Value", "Source", "rel 1σ"], const_rows, title="=== CONSTANTS (values, sources, 1σ) ===")
    show(block)

    # Optional: save consolidated report (+ SHA-256 checksum of content)
    if SAVE_REPORT:
        out_text = "\n".join(blocks) + "\n"
        content_bytes = out_text.encode("utf-8")
        sha256 = hashlib.sha256(content_bytes).hexdigest()
        footer = _format_table(
            ["Field", "Value"],
            [["SHA-256 (of above content)", sha256], ["Bytes (UTF-8)", len(content_bytes)], ["Lines", out_text.count("\n")]],
            title="=== CHECKSUM ==="
        )
        out_text_with_checksum = out_text + footer + "\n"
        with open(SAVE_REPORT, "w", encoding="utf-8") as f:
            f.write(out_text_with_checksum)
        print(f"\nSaved report to: {os.path.abspath(SAVE_REPORT)}")
        print(f"SHA-256: {sha256}")

if __name__ == "__main__":
    main()
