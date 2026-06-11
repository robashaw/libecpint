#!/usr/bin/env python3
"""
Independent, closed-form reference for on-centre (single-atom) ECP integrals.

For an ECP and two Cartesian Gaussian shells all sharing the same centre, the ECP
integral separates exactly into an angular part (monomials integrated over the unit
sphere) and a radial part (half-integer Gamma functions). This is computed here with
sympy in exact arithmetic, completely independently of libecpint's Bessel/recursion
machinery, and emitted to ref_atom.output for tests/lib/ref_atom/test.cpp to check.

Operator (matching what libecpint::ECPIntegral::compute_shell_pair evaluates):

    U = U_L(r) . 1  +  sum_{l=0}^{L-1} U_l(r) P_l

with U_l(r) = sum_k d_k r^(m_k - 2) exp(-zeta_k r^2), where m_k is the power passed to
ECP::addPrimitive (libecpint stores m_k - 2 and evaluates that power; verified empirically),
and P_l = sum_m |lm><lm| the projector onto angular momentum l about the centre.

Cartesian component ordering matches libecpint: for x = L..0, y = L-x..0, z = L-x-y.
"""
import sympy as sp
from sympy import Rational, gamma, factorial2, legendre, symbols, expand, Poly, pi

mp_dps = 25

# ---- angular: integral over the unit sphere of n_x^i n_y^j n_z^k ----
def sphere_monomial(i, j, k):
    if (i % 2) or (j % 2) or (k % 2):
        return sp.Integer(0)
    # 4*pi * (i-1)!! (j-1)!! (k-1)!! / (i+j+k+1)!!   (with (-1)!! = 1)
    return 4*pi * factorial2(i-1)*factorial2(j-1)*factorial2(k-1) / factorial2(i+j+k+1)

# ---- radial: int_0^inf r^P exp(-c r^2) dr = 1/2 Gamma((P+1)/2) / c^((P+1)/2) ----
def radial_int(P, c):
    return Rational(1, 2) * gamma(Rational(P+1, 2)) / c**Rational(P+1, 2)

nx, ny, nz, mx, my, mz = symbols('nx ny nz mx my mz')

def angular_double(l, A, B):
    """(2l+1)/(4pi) * int dOmega dOmega' P_l(n.n') n^A n'^B, with A=(a,b,c), B=(d,e,f)."""
    a, b, c = A
    d, e, f = B
    t = nx*mx + ny*my + nz*mz
    poly = expand(legendre(l, t) * nx**a*ny**b*nz**c * mx**d*my**e*mz**f)
    total = sp.Integer(0)
    p = Poly(poly, nx, ny, nz, mx, my, mz)
    for monom, coeff in p.terms():
        i1, j1, k1, i2, j2, k2 = monom
        total += coeff * sphere_monomial(i1, j1, k1) * sphere_monomial(i2, j2, k2)
    return (2*l + 1) / (4*pi) * total

def cart_list(L):
    out = []
    for x in range(L, -1, -1):
        for y in range(L - x, -1, -1):
            out.append((x, y, L - x - y))
    return out

def prim_integral(A, alpha, B, beta, ecp):
    """ECP integral for one primitive pair; A,B cart tuples; ecp: dict l -> list of (m,zeta,d)."""
    L = max(ecp.keys())
    lA, lB = sum(A), sum(B)
    a, b, c = A; d, e, f = B
    # type 1: local channel L applied with no projector
    ang1 = sphere_monomial(a+d, b+e, c+f)
    rad1 = sum(dk * radial_int(lA + lB + mk, alpha + beta + zk) for (mk, zk, dk) in ecp[L])
    val = ang1 * rad1
    # type 2: semilocal channels l < L
    for l in range(L):
        if l not in ecp:
            continue
        t2ang = angular_double(l, A, B)
        if t2ang == 0:
            continue
        rad = sum(dk * radial_int(lA + lB + mk, alpha + beta + zk) for (mk, zk, dk) in ecp[l])
        val += t2ang * rad
    return val

def shell_pair(LA, primsA, LB, primsB, ecp):
    """Return flat list (libecpint order) of contracted integrals for shells A,B."""
    cartsA, cartsB = cart_list(LA), cart_list(LB)
    out = []
    for A in cartsA:
        for B in cartsB:
            s = sp.Integer(0)
            for (alpha, cA) in primsA:
                for (beta, cB) in primsB:
                    s += cA * cB * prim_integral(A, alpha, B, beta, ecp)
            out.append(s)
    return out

if __name__ == "__main__":
    mp = sp  # exact, evaluated to float at the end
    a = lambda s: sp.nsimplify(s)
    # ECP at origin, L = 2 (local = l=2, semilocal = l=0,1)
    # --- Scenario A: moderate exponents, mixed channels (general sanity) ---
    ecp_A = {
        0: [(2, a('1.5'), a('2.0')),  (3, a('0.7'), a('0.4'))],
        1: [(2, a('1.2'), a('-1.0'))],
        2: [(2, a('0.9'), a('0.5')),  (4, a('0.6'), a('0.25'))],
    }
    s_A = (0, [(a('1.1'), a('1.0'))])
    p_A = (1, [(a('0.8'), a('1.0')), (a('2.5'), a('0.7'))])
    d_A = (2, [(a('0.7'), a('1.0'))])

    # --- Scenario B: STEEP local channel (g up to 1e5), regression for the type-1 grid-scale
    #     collapse (g >> p) and the l>0 screening-gate collapse at large p+g. L=0 so this is a
    #     pure local (type-1) potential exercised across s/p/d shells and diffuse AOs. ---
    ecp_B = {0: [(2, a('1000.0'), a('1.0')), (2, a('100000.0'), a('0.5'))]}
    s_B = (0, [(a('0.05'), a('1.0'))])
    p_B = (1, [(a('0.05'), a('1.0'))])
    d_B = (2, [(a('0.5'),  a('1.0'))])

    # --- Scenario C: diffuse AO pair vs a standard local primitive, g/p ~ 250 (the "ordinary"
    #     augmented-basis diffuse-diffuse regime the screening bug silently zeroed out). ---
    ecp_C = {0: [(2, a('10.0'), a('1.0'))]}
    s_C = (0, [(a('0.02'), a('1.0'))])
    p_C = (1, [(a('0.02'), a('1.0'))])

    scenarios = [
        (ecp_A, [("A_ss", s_A, s_A), ("A_ps", p_A, s_A), ("A_pp", p_A, p_A),
                 ("A_ds", d_A, s_A), ("A_dp", d_A, p_A), ("A_dd", d_A, d_A)]),
        (ecp_B, [("B_ss", s_B, s_B), ("B_ps", p_B, s_B), ("B_pp", p_B, p_B),
                 ("B_ds", d_B, s_B), ("B_dp", d_B, p_B), ("B_dd", d_B, d_B)]),
        (ecp_C, [("C_ss", s_C, s_C), ("C_ps", p_C, s_C), ("C_pp", p_C, p_C)]),
    ]

    flat = []
    for ecp, pairs in scenarios:
        for name, (LA, pa), (LB, pb) in pairs:
            flat.extend(shell_pair(LA, pa, LB, pb, ecp))

    with open("ref_atom.output", "w") as fh:
        fh.write("\n".join("{:.15e}".format(float(sp.N(v, mp_dps))) for v in flat))
    print(f"wrote {len(flat)} reference values")
