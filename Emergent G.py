#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2025 by Mehrdad Pajuhaan (pajuhaan@gmail.com)
# =============================================================================
#  Relator — Emergent Newton's Constant G from the electron mass (Path A & B)
#  ---------------------------------------------------------------------------
#  Goal (reviewer-facing):
#    Given the *experimental* electron mass m_e (MeV) and the closed, alpha-locked
#    dimensionless geometry/micro blocks from the Relator program, compute the
#    emergent Newton constant G via *inverting* the ring–Planck prefactor.
#
#  Key inversion (shared skeleton):
#    ℓ_P = sqrt(ħ G / c^3)  ⇒  1/ℓ_P = sqrt(c^3/(ħ G))
#    m_path(MeV) = [ D ħ /(c ℓ_P) ] · F_path(α,geometry,...) · (c^2 / e) / 1e6
#                 = ( D·√ħ·c^(5/2) / (e·1e6) ) · F_path · G^(−1/2)
#    ⇒  G_path = [ D·√ħ·c^(5/2) · F_path / ( m_e(MeV) · e · 1e6 ) ]^2
#
#  Paths:
#    Path A (scalar):
#      F_A = [ 3 N_eff^eff ζ_C^4 ρ^3 e^{−σ_A/α} / (32π α K_EM^eff) ]^(1/4),
#      with K_EM^eff = K_EM^raw · softening(κ), and a *local* softening-only
#      micro on K_EM (use ζ_soft from Λ_ind to avoid double counting bulk ζ(α)).
#
#    Path B (vector, Ward-preserving product):
#      F_B = sqrt{ ρ K_T^eff e^{−σ_B/α} (1+f_M1) / [ 2 α K_EM^{raw,eff} ] },
#      with K_T^eff = K_T · (1 − ½ ρ² ζ_geom)(1 − ½ ρ² ζ(α)),
#           K_EM^{raw,eff} = K_EM^raw · (1 + ¼ ρ² ζ_geom)(1 + ¼ ρ² ζ(α)),
#           f_M1 = (α/2) e^{−y*²/2} y*².
#
#  Outputs:
#    • G_A, G_B, and geometric mean G_geo = sqrt(G_A G_B)
#    • Comparison vs CODATA G (Δ, ppm)
#    • Minimal audit block (y*, S, D, ρ, ζ_geom, ζ(α), etc.)
#
#  Notes:
#    - α is taken as an input constant (no fitting).
#    - No use of experimental G inside the pipeline beyond the final comparison table.
# =============================================================================

from dataclasses import dataclass
import mpmath as mp
import os

# ============================= Global precision ==============================
mp.mp.dps = 140

# ========================= Physical constants (SI) ===========================
# Values and uncertainties (1σ) aligned with SI/CODATA 2022. Sources are noted.
PHYS = {
    "c": {  # exact by SI
        "value": mp.mpf('299792458.0'),
        "source": "SI definition (exact)",
        "u_rel": mp.mpf('0')
    },
    "e": {  # exact by SI
        "value": mp.mpf('1.602176634e-19'),
        "source": "SI definition (exact)",
        "u_rel": mp.mpf('0')
    },
    "hbar": {  # exact via fixed h; numerical rounding only
        "value": mp.mpf('1.054571817e-34'),
        "source": "SI (via fixed h); numerical rounding only",
        "u_rel": mp.mpf('0')
    },
    "G_CODATA": {  # for comparison only
        "value": mp.mpf('6.67430e-11'),
        "source": "CODATA 2022 (comparison only)",
        "u_rel": mp.mpf('2.2e-5')  # ~22 ppm
    },
    "alpha_in": {  # α is an input to the model (not fitted)
        "value": mp.mpf('7.2973525643e-3'),
        "source": "CODATA 2022 (input)",
        "u_rel": mp.mpf('1.5e-10')
    },
    "m_e_exp_MeV": {  # used to *invert* for G
        "value": mp.mpf('0.51099895069'),
        "source": "CODATA 2022 (electron mass, MeV)",
        "u_rel": mp.mpf('0')
    }
}

# ================================ Utilities ==================================
def nstr(x, n=28): return mp.nstr(x, n)
def ppm(x): return mp.mpf('1e6') * x

# ====================== Paper notation: shared constants =====================
pi      = mp.pi
C0_UNI  = (1/pi) * (mp.mpf('4')/3 + 1/(4*pi**2))  # C0^{uni}
C0_GAUSS= mp.mpf('0.5') * (mp.log(2) + mp.euler)  # ½(ln2+γ)
EPSILON = 1/mp.sqrt(pi)                           # ε = 1/√π
ETA0    = 1/pi                                    # η = 1/π
ELL0    = EPSILON*ETA0                            # ℓ0 = εη
ITOT    = mp.mpf('1')/6 - mp.mpf('1')/(4*pi**2)   # ∫_0^1 x^2 sin^2(πx)dx

# ============================= Ring integrals (S¹) ===========================
def erfcx(z): return mp.e**(z**2) * mp.erfc(z)  # scaled erfc

def a_of_delta(Delta): return 2*mp.sin(Delta/2)

def J_of(a, y):
    # J(a,y)=√π/(2y)·erfcx(a/(2y))
    return mp.sqrt(pi)/(2*y) * erfcx(a/(2*y))

def K0_of(y):
    # K0(y)=(1/π)∫_0^π J(a(Δ),y)dΔ
    return mp.quad(lambda D: J_of(a_of_delta(D), y), [0, pi]) / pi

def Lam_of(y):
    # Λ(y)=(1/π)∫_0^π cosΔ·J(a(Δ),y)dΔ
    return mp.quad(lambda D: mp.cos(D) * J_of(a_of_delta(D), y), [0, pi]) / pi

def S_of(y, alpha):
    # S(y;α)=2 + ζ_geom(y) + α K0(y)/(8π²),  ζ_geom=[K0/(2π²)]Λ
    K0 = K0_of(y); Lam = Lam_of(y)
    zeta_geom = (K0/(2*pi**2))*Lam
    S = (2 + zeta_geom) + alpha * K0 / (8*pi**2)
    return S, zeta_geom, K0, Lam

def dS_dy(y, alpha):
    h = mp.mpf('5e-7')*(1+abs(y))
    Sp,_,_,_ = S_of(y+h, alpha); Sm,_,_,_ = S_of(y-h, alpha)
    return (Sp - Sm)/(2*h)

def F_shape(y, alpha):
    # Stationarity: 3 S'/S + (2/y − y) = 0
    S,_,_,_ = S_of(y, alpha)
    return 3*dS_dy(y, alpha)/S + (2/y - y)

def solve_y_star(alpha, y_lo=mp.mpf('0.35'), y_hi=mp.mpf('3.0'), step=mp.mpf('0.02')):
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

def precompute_spectrum(M: int = 30, tol: mp.mpf = mp.mpf('1e-40')):
    K = series_K(tol=tol)
    L_list = []
    if M >= 2:
        for m in range(2, M+1):
            L_list.append((m, L_2m(m, tol=tol)))
    return K, L_list

def DC_of_alpha_fixedK(alpha: mp.mpf, K: mp.mpf, L_list, M: int = 30) -> mp.mpf:
    # 𝓓_C(α) = (α/π)[ √(1−ξ) − (ξ/2)K − Σ_{m≥2} (ξ/2)^m L_{2m} ], ξ=2 C0^{uni} α
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
    # Tangential field on the unit sphere (x=cosθ), evaluated at R=r*.
    s = mp.sqrt(1 - x**2); rho = rstar*s; z = rstar*x
    return _B_rho(rho, z) * x - _B_z(rho, z) * s

def Lambda_OUT_exact(eta: mp.mpf, lmax: int = 23, gl_nodes: int = 512) -> mp.mpf:
    # Exact OUT energy via spherical multipoles up to lmax; then Λ_out = −2 U_out.
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

def Lambda_OUT_extrapolated(eta: mp.mpf, lbase: int = 23) -> mp.mpf:
    def _eval(lm): return Lambda_OUT_exact(eta, lmax=lm, gl_nodes=512)
    S1 = _eval(lbase); S2 = _eval(lbase+2); S3 = _eval(lbase+4)
    denom = (S3 - 2*S2 + S1)
    if denom == 0 or abs(denom) < mp.mpf('1e-30')*max(1,abs(S3)): return S3
    Sout = S1 - (S2 - S1)**2 / denom
    if abs(Sout) > 10*max(abs(S1),abs(S2),abs(S3)): return S3
    return Sout

# ============================= Λ-chain & ζ(α) ================================
def P_IR_chi(ell: mp.mpf) -> mp.mpf:
    # IR acceptance P_IR^(χ)(ℓ) with exact S¹ weight (no small-κ expansions).
    w = lambda x: (x**2) * (mp.sin(pi*x)**2)
    num = mp.quad(lambda x: w(x) * (1 - mp.mpf('1')/3 * ((1 - x)**2/(((1 - x)**2) + ell**2))) *
                            mp.e**(-((1 - x)/ell)**2), [0, 1])
    return num / ITOT

def Lambda_ind(eps: mp.mpf) -> mp.mpf:
    return mp.log(8/eps) - 2

def curvature_series_eta(eta: mp.mpf, order: int = 20) -> mp.mpf:
    if order < 2: return mp.mpf('0')
    s = mp.mpf('0'); max_k = int(order // 2)
    for k in range(1, max_k + 1):
        s += (eta**(2*k)) / mp.factorial(2*k + 1)
    return s

def gamma_geom(eta: mp.mpf) -> mp.mpf:
    return mp.mpf('0.5') * (1 + curvature_series_eta(eta, order=20))

def gamma_eff(eta: mp.mpf, K: mp.mpf, DC_lock: mp.mpf, P_ir: mp.mpf) -> mp.mpf:
    return gamma_geom(eta) + (K / (2 * DC_lock)) * C0_GAUSS * P_ir

def build_Lambda_eff(alpha: mp.mpf, K: mp.mpf, L_list):
    # Full Λ_eff with exact OUT + L-extrap; χ/self ladders folded in compactly.
    P_ir  = P_IR_chi(ELL0)
    Lin   = Lambda_ind(EPSILON)
    dUV   = C0_GAUSS * P_ir
    dOUT  = Lambda_OUT_extrapolated(ETA0, lbase=23)
    base  = Lin + dUV + dOUT
    DC    = DC_of_alpha_fixedK(alpha, K, L_list, M=30)
    g_eff = gamma_eff(ETA0, K, DC, P_ir)
    dSYNC = g_eff * P_ir * dOUT
    Lam_eff = base + dSYNC
    return {"Lambda_eff": Lam_eff, "Lambda_ind": Lin, "P_ir": P_ir, "Delta_out": dOUT, "DC": DC}

# ========================= Ring-level kernels (Path A/B) =====================
def kappa_tension(x): return (x**2)/(1+x**2)

def softening_factor(kappa):
    if kappa <= mp.mpf('0'): return mp.mpf('1')
    if kappa < mp.mpf('1e-20'): return mp.mpf('1') - kappa/2 + (3*(kappa**2))/4
    return mp.sqrt(mp.pi)/mp.sqrt(kappa) * erfcx(1/mp.sqrt(kappa))

def K_EM_raw_closed(x):
    # K_EM^{raw}(x) = (2x/π) K(m=-4x²)
    m = -4*(x**2)
    return (2*x/pi) * mp.ellipk(m)

def K_EM_eff_closed(x): return K_EM_raw_closed(x) * softening_factor(kappa_tension(x))
def N_eff_closed(x):
    κ = kappa_tension(x)
    return 2*(1+κ)*(1 - mp.mpf('0.5')*κ**2)
def K_T_closed(x): return 2*kappa_tension(x)

# ============================= Geometric RG locks ============================
def beta0_from_DS(D_over_S): return (6/pi)*(D_over_S)       # 4/3 → 8/π
def sigmas_from_beta0(beta0): return 4/beta0, 2/beta0       # (σ_A, σ_B) = (π/2, π/4)

# =============================== F-path builders =============================
def build_blocks_and_F(alpha):
    # Shared geometry
    y_star = solve_y_star(alpha)
    S_star, zeta_geom, K0_y, Lam_y = S_of(y_star, alpha)
    D      = (mp.mpf('4')/3)*S_star
    rho    = 1/D
    x_rel  = rho / y_star

    # Kernels
    K_EM_raw = K_EM_raw_closed(x_rel)
    K_EM_eff = K_EM_eff_closed(x_rel)
    N_eff    = N_eff_closed(x_rel)
    K_T      = K_T_closed(x_rel)

    # RG locks
    beta0 = beta0_from_DS(D/S_star)     # = 8/π exactly here
    sigma_A, sigma_B = sigmas_from_beta0(beta0)

    # Spectrum & ζ(α)
    K_spec, L_list = precompute_spectrum(M=30, tol=mp.mpf('1e-40'))
    LAMBDA = build_Lambda_eff(alpha, K_spec, L_list)
    zeta_alpha = (K_spec/(2*pi**2)) * LAMBDA["Lambda_eff"]

    # Scalar S4 completion (Path A)
    D_C = LAMBDA["DC"]
    D_phys_S4 = D_C - (alpha/pi)*zeta_alpha + ((alpha*zeta_alpha)/pi)**2 / (4*D_C)
    N_eff_eff = N_eff * (1 + rho**2 * D_phys_S4)

    # Local ζ_soft for Gaussian EM-soft only (avoid double-count with bulk ζ(α))
    zeta_soft = (K_spec/(2*pi**2)) * LAMBDA["Lambda_ind"]
    K_EM_eff_soft = K_EM_eff * (1 - mp.mpf('0.5')*rho**2 * zeta_soft)

    # Vector (Path B) Ward-preserving product
    K_T_eff   = K_T * (1 - mp.mpf('0.5') * rho**2 * zeta_geom) * (1 - mp.mpf('0.5') * rho**2 * zeta_alpha)
    K_EMr_eff = K_EM_raw * (1 + mp.mpf('0.25') * rho**2 * zeta_geom) * (1 + mp.mpf('0.25') * rho**2 * zeta_alpha)

    # F_A and F_B (dimensionless)
    ZETA_C = pi
    F_A = ( 3*N_eff_eff*(ZETA_C**4)*(rho**3)*mp.e**(-sigma_A/alpha) / (32*pi*alpha*K_EM_eff_soft) )**(mp.mpf('1')/4)
    fM1 = (alpha/2) * mp.e**(-(y_star**2)/2) * (y_star**2)
    F_B = mp.sqrt( (rho*K_T_eff*mp.e**(-sigma_B/alpha)*(1+fM1)) / (2*alpha*K_EMr_eff) )

    return {
        "y_star": y_star, "S": S_star, "D": D, "rho": rho, "x_rel": x_rel,
        "zeta_geom": zeta_geom, "zeta_alpha": zeta_alpha,
        "K_EM_raw": K_EM_raw, "K_EM_eff": K_EM_eff, "K_EM_eff_soft": K_EM_eff_soft,
        "N_eff": N_eff, "N_eff_eff": N_eff_eff, "K_T": K_T, "K_T_eff": K_T_eff, "K_EMr_eff": K_EMr_eff,
        "sigma_A": sigma_A, "sigma_B": sigma_B,
        "F_A": F_A, "F_B": F_B
    }

def invert_for_G(F_path, D, m_e_MeV, HBAR, C, ECHG):
    # G_path = [ D·√ħ·c^(5/2) · F_path / ( m_e(MeV) · e · 1e6 ) ]^2
    num = D * mp.sqrt(HBAR) * (C**(mp.mpf('2.5'))) * F_path
    den = m_e_MeV * ECHG * mp.mpf('1e6')
    return (num/den)**2

# ============================== Reporting layer ==============================
def _format_table(headers, rows, title=None):
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
            for j in range(cols): widths[j] = max(widths[j], len(str(r[j])))
        def line(vals): return "  ".join(str(vals[j]).ljust(widths[j]) for j in range(cols))
        out = []
        if title: out.append(title)
        out.append(line(headers))
        out.append(line(["-"*w for w in widths]))
        for r in rows: out.append(line(r))
        return "\n".join(out)

def main():
    # Inputs
    alpha = PHYS["alpha_in"]["value"]
    m_e   = PHYS["m_e_exp_MeV"]["value"]
    HBAR  = PHYS["hbar"]["value"]
    C     = PHYS["c"]["value"]
    ECHG  = PHYS["e"]["value"]

    # Build dimensionless blocks and F_A, F_B
    B = build_blocks_and_F(alpha)

    # Invert for G (both paths)
    G_A = invert_for_G(B["F_A"], B["D"], m_e, HBAR, C, ECHG)
    G_B = invert_for_G(B["F_B"], B["D"], m_e, HBAR, C, ECHG)
    G_geo = mp.sqrt(G_A * G_B)

    # Comparison vs CODATA (optional, for context)
    G_cod = PHYS["G_CODATA"]["value"]
    rows_G = [
        ["G_A (Path A)", nstr(G_A, 16), nstr((G_A-G_cod), 12), nstr(ppm((G_A-G_cod)/G_cod), 10)],
        ["G_B (Path B)", nstr(G_B, 16), nstr((G_B-G_cod), 12), nstr(ppm((G_B-G_cod)/G_cod), 10)],
        ["G_geo = √(G_A G_B)", nstr(G_geo, 16), nstr((G_geo-G_cod), 12), nstr(ppm((G_geo-G_cod)/G_cod), 10)],
    ]
    print(_format_table(["Quantity", "Value [m^3 kg^-1 s^-2]", "Δ vs CODATA", "Δ [ppm]"], rows_G,
                        title="=== EMERGENT NEWTON'S CONSTANT FROM m_e (THIS RUN) ==="))
    print()

    # Minimal audit block (helps reproducibility without cluttering)
    rows_audit = [
        ["α (input)",           nstr(alpha, 20)],
        ["y* (shape root)",     nstr(B["y_star"], 22)],
        ["S(y*)",               nstr(B["S"], 22)],
        ["D = 4S/3",            nstr(B["D"], 22)],
        ["ρ = 1/D",             nstr(B["rho"], 22)],
        ["x = ρ/y*",            nstr(B["x_rel"], 22)],
        ["ζ_geom(y*)",          nstr(B["zeta_geom"], 22)],
        ["ζ(α) bulk",           nstr(B["zeta_alpha"], 22)],
        ["σ_A, σ_B",            f"{nstr(B['sigma_A'],16)} ; {nstr(B['sigma_B'],16)}"],
        ["F_A, F_B",            f"{nstr(B['F_A'],20)} ; {nstr(B['F_B'],20)}"],
    ]
    print(_format_table(["Quantity", "Value"], rows_audit, title="=== AUDIT — KEY DIMENSIONLESS BLOCKS ==="))
    print()

    # Constants table
    rows_const = []
    for key in ["c","e","hbar","alpha_in","m_e_exp_MeV","G_CODATA"]:
        d = PHYS[key]
        tag = key if key != "alpha_in" else "alpha_in (input)"
        if key == "m_e_exp_MeV": tag = "m_e_exp_MeV (input mass)"
        if key == "G_CODATA": tag = "G_CODATA (comparison)"
        rows_const.append([tag, nstr(d["value"], 20), d["source"], nstr(d["u_rel"], 12)])
    print(_format_table(["Constant", "Value", "Source", "rel 1σ"], rows_const, title="=== CONSTANTS (VALUES, SOURCES, 1σ) ==="))

if __name__ == "__main__":
    main()
