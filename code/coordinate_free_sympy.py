#!/usr/bin/env python3
"""Exact symbolic checks for the 4D single-scalar DBI stress-tensor flow.

The script verifies, with SymPy, the matrix identities used in the
coordinate-free derivation.  It deliberately does not use neural-network
output or floating-point fitting.

Conventions:
  g = diag(-,+,+,+)
  u_mu = sqrt(lambda) partial_mu phi
  K^mu_nu = sqrt(lambda) F^mu_nu  (g-skew-adjoint)
  P^mu_nu = u^mu u_nu
  C = I + K + P
  D^2 = det C
  M = I - lambda T

Temporary scalar contractions used only in the proof:
  s     = u.u
  tau   = (1/2) tr K^2
  c     = u.K^2.u
  delta = det K
"""

from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

OUT = Path('/mnt/data')
SUMMARY = OUT / 'dbi_4d_coordinate_free_sympy_summary.json'

# ------------------------------------------------------------------
# 1. Fully generic 4x4 Lorentz-covariant symbolic matrix check
# ------------------------------------------------------------------
a = sp.symbols('a')
u0,u1,u2,u3 = sp.symbols('u0 u1 u2 u3')
f01,f02,f03,f12,f13,f23 = sp.symbols('f01 f02 f03 f12 f13 f23')

g = sp.diag(-1,1,1,1)
Bcov = sp.Matrix([
    [0,   f01, f02, f03],
    [-f01, 0,  f12, f13],
    [-f02,-f12, 0,  f23],
    [-f03,-f13,-f23,0],
])
K = g * Bcov
u = sp.Matrix([u0,u1,u2,u3])
uflat = u.T * g
P = u * uflat
I = sp.eye(4)

s = sp.expand((uflat*u)[0])
tau = sp.expand(sp.trace(K*K)/2)
delta = sp.expand(K.det())
c = sp.expand((uflat*K*K*u)[0])

# Four-dimensional Cayley-Hamilton identity for a g-skew endomorphism.
CH = K**4 - tau*K**2 + delta*I
assert all(sp.simplify(entry) == 0 for entry in CH)

# Matrix determinant lemma + Cayley-Hamilton prediction.
p_direct = sp.expand((a*I + K + P).det())
p_formula = sp.expand(
    a**4 + s*a**3 - tau*a**2 + (c-s*tau)*a + delta
)
assert sp.factor(p_direct-p_formula) == 0

# Rank-one determinant identity for I+P.
assert sp.factor((I+P).det()-(1+s)) == 0

# ------------------------------------------------------------------
# 2. Abstract exact algebra for E1, E3, E4 and the DBI lambda-flow
# ------------------------------------------------------------------
S, TAU, C2, DELTA, D = sp.symbols('s tau c delta D')
D2 = 1 + S - TAU + C2 - S*TAU + DELTA

# p'(1) = D E1.
DE1 = 4 + 3*S - 2*TAU + C2 - S*TAU
# E3 = E4 tr(M^{-1}), with E4=1+s.
DE3 = (1+S)*(4+S-2*TAU) + C2
E4 = 1+S

# At fixed physical fields, S and TAU scale as lambda, whereas
# C2-S*TAU+DELTA scales as lambda^2. Therefore
# 2 D (y+1) = 2+s-tau, y=lambda^2 Q.
flow_factor = 2 + S - TAU

# Main pre-division identity: valid also at S=0.
assert sp.simplify((DE1-DE3) + S*flow_factor) == 0

# Hence E1-E3 = 2(1-E4)(y+1).
# The rational formula follows when E4 != 1.

# ------------------------------------------------------------------
# 3. Characteristic-polynomial factorization and the invariant identity
# ------------------------------------------------------------------
alpha,beta,r,z = sp.symbols('alpha beta r z')
chi = sp.expand((z**2-alpha*z+r)*(z**2-beta*z+1))
# chi = z^4 - E1 z^3 + E2 z^2 - E3 z + E4
Ee1 = sp.expand(alpha+beta)
Ee2 = sp.expand(1+r+alpha*beta)
Ee3 = sp.expand(alpha+beta*r)
Ee4 = r

identity = sp.expand(
    (Ee4-1)**2*Ee2
    -(-Ee1**2*Ee4 + Ee1*Ee3*(Ee4+1) - Ee3**2
      + Ee4**3 - Ee4**2 - Ee4 + 1)
)
assert sp.factor(identity) == 0

# beta is the coefficient of the quadratic with unit constant term.
# For r != 1: beta=(E1-E3)/(1-r)=2(y+1).
assert sp.factor((Ee1-Ee3)-beta*(1-r)) == 0

# Null locus r=1: the two factors both have unit constant term.
# E1=alpha+beta and E2=2+alpha beta, so the discriminant is
# (alpha-beta)^2.
null_disc = sp.expand((alpha+beta)**2 - 4*((2+alpha*beta)-2))
assert sp.factor(null_disc-(alpha-beta)**2) == 0

# ------------------------------------------------------------------
# 4. Convert the generic formula to invariants of T itself
# ------------------------------------------------------------------
lam = sp.symbols('lambda', nonzero=True)
s1,s2,s3,s4 = sp.symbols('e1T e2T e3T e4T')
ME1 = 4-lam*s1
ME3 = 4-3*lam*s1+2*lam**2*s2-lam**3*s3
ME4 = 1-lam*s1+lam**2*s2-lam**3*s3+lam**4*s4

y_from_M = sp.factor((ME1-ME3+2*ME4-2)/(2*(1-ME4)))
Q_from_T = sp.factor((2*lam*s4-s3)/(2*(s1-lam*s2+lam**2*s3-lam**3*s4)))
assert sp.simplify(y_from_M-lam**2*Q_from_T) == 0

# ------------------------------------------------------------------
# 5. Independent exact rational examples using the actual matrices
# ------------------------------------------------------------------
def exact_example(vals_u, vals_f):
    sub = {
        u0: sp.Rational(vals_u[0]), u1: sp.Rational(vals_u[1]),
        u2: sp.Rational(vals_u[2]), u3: sp.Rational(vals_u[3]),
        f01: sp.Rational(vals_f[0]), f02: sp.Rational(vals_f[1]),
        f03: sp.Rational(vals_f[2]), f12: sp.Rational(vals_f[3]),
        f13: sp.Rational(vals_f[4]), f23: sp.Rational(vals_f[5]),
    }
    Kn = K.subs(sub)
    Pn = P.subs(sub)
    Cn = I + Kn + Pn
    Cdn = I - Kn + Pn
    D2n = sp.factor(Cn.det())
    Dn = sp.sqrt(D2n)
    sn = sp.factor(s.subs(sub))
    taun = sp.factor(tau.subs(sub))
    cn = sp.factor(c.subs(sub))
    deltan = sp.factor(delta.subs(sub))
    predicted_D2 = sp.factor(1+sn-taun+cn-sn*taun+deltan)
    assert sp.factor(D2n-predicted_D2) == 0

    # Actual mixed matrix M=D/2(C^{-1}+C^{-dagger}) for this exact sample.
    Mn = sp.simplify(Dn * sp.Rational(1,2) * (Cn.inv()+Cdn.inv()))
    coeff = Mn.charpoly().all_coeffs()  # 1,-E1,E2,-E3,E4
    E1n = sp.factor(-coeff[1])
    E2n = sp.factor(coeff[2])
    E3n = sp.factor(-coeff[3])
    E4n = sp.factor(coeff[4])

    # y=lambda^2 Q; the exact scaling identity gives y+1=(2+s-tau)/(2D).
    yn = sp.factor((2+sn-taun)/(2*Dn)-1)
    assert sp.simplify(E1n-E3n-2*(1-E4n)*(yn+1)) == 0
    algn = sp.factor(
        (E4n-1)**2*E2n
        -(-E1n**2*E4n + E1n*E3n*(E4n+1) - E3n**2
          + E4n**3-E4n**2-E4n+1)
    )
    assert algn == 0

    null_check = None
    if sn == 0:
        disc = sp.factor(E1n**2-4*E2n+8)
        # On the standard real branch c<=0, the plus root is physical.
        null_check = sp.simplify((E1n+sp.sqrt(disc))/4-1-yn)
        assert null_check == 0

    return {
        's': str(sn), 'tau': str(taun), 'c': str(cn),
        'detK': str(deltan), 'D2': str(D2n),
        'matrix_predivision_identity': 'exact zero',
        'matrix_algebraic_constraint': 'exact zero',
        'null_flow_formula': ('exact zero' if sn == 0 else 'not applicable'),
    }

examples = [
    ((sp.Rational(1,4),0,0,0), (sp.Rational(1,5),sp.Rational(1,7),0,sp.Rational(1,6),0,sp.Rational(1,8))),
    ((sp.Rational(1,3),sp.Rational(1,5),sp.Rational(1,7),0),
     (sp.Rational(1,4),0,sp.Rational(1,9),sp.Rational(1,6),sp.Rational(1,10),sp.Rational(1,8))),
    ((sp.Rational(1,4),sp.Rational(1,4),0,0),
     (sp.Rational(1,5),sp.Rational(1,8),sp.Rational(1,9),sp.Rational(1,7),sp.Rational(1,11),sp.Rational(1,6))),
]
example_results = [exact_example(*ex) for ex in examples]

summary = {
    'generic_4x4_cayley_hamilton': 'exact zero',
    'generic_4x4_determinant_polynomial': 'exact zero',
    'det_I_plus_rank1': 'exact zero',
    'main_predivision_identity': 'exact zero',
    'characteristic_factorization_identity': 'exact zero',
    'null_discriminant_identity': 'exact zero',
    'conversion_to_T_invariants': 'exact zero',
    'exact_matrix_level_checks': 'exact zero on all rational examples',
    'exact_rational_examples': example_results,
}
SUMMARY.write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
print(f'Wrote {SUMMARY}')
