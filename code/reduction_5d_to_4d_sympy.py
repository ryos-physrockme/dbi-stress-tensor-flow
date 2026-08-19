#!/usr/bin/env python3
"""Exact and numerical checks for 5D BI -> 4D single-scalar DBI reduction.

The purpose is to verify independently that the known five-dimensional
Born-Infeld stress-tensor flow reduces to the four-dimensional single-scalar
DBI stress-tensor flow discovered in the preceding calculation.

Conventions:
  g4 = diag(-,+,+,+), g5 = diag(-,+,+,+,+)
  y is a compact spacelike coordinate, partial_y A_M = 0, g_yy = 1
  phi = A_y
  u^mu = sqrt(lambda) partial^mu phi
  K^mu_nu = sqrt(lambda) F^mu_nu
  C = I + K + u u_flat
  D^2 = det C
  M = I - lambda T_4

The five-dimensional mixed BI matrix before taking the symmetric inverse is
  C5 = [[I+K, u],[-u_flat, 1]].
Its Schur complement is C, so the reduced action is exactly 4D single-scalar
DBI.  All symbolic checks below are exact SymPy identities.
"""

from __future__ import annotations
import json
from pathlib import Path
import math
import numpy as np
import pandas as pd
import sympy as sp

OUT = Path('/mnt/data')
JSON = OUT / 'dbi_5d_to_4d_reduction_sympy_summary.json'
CSV = OUT / 'dbi_5d_to_4d_reduction_validation.csv'

u0,u1,u2,u3 = sp.symbols('u0 u1 u2 u3')
f01,f02,f03,f12,f13,f23 = sp.symbols('f01 f02 f03 f12 f13 f23')

g4 = sp.diag(-1,1,1,1)
g5 = sp.diag(-1,1,1,1,1)
Fcov = sp.Matrix([
    [0, f01, f02, f03],[-f01,0,f12,f13],[-f02,-f12,0,f23],[-f03,-f13,-f23,0]
])
K = g4 * Fcov
u = sp.Matrix([u0,u1,u2,u3])
uf = u.T * g4
P = u * uf
I4 = sp.eye(4)
I5 = sp.eye(5)
C = I4 + K + P
Cdag = I4 - K + P

C5 = sp.zeros(5)
C5[:4,:4] = I4 + K
C5[:4,4] = u
C5[4,:4] = -uf
C5[4,4] = 1

# 5D Born-Infeld determinant -> 4D DBI determinant.
assert sp.factor(C5.det() - C.det()) == 0

s = sp.expand((uf*u)[0])
tau = sp.expand(sp.trace(K*K)/2)
c = sp.expand((uf*K*K*u)[0])
delta = sp.expand(K.det())
D2_formula = sp.expand(1+s-tau+c-s*tau+delta)
assert sp.factor(C.det() - D2_formula) == 0

# R = u_flat C^{-1} u; avoid a generic inverse using adj(C).
R_num = sp.expand((uf * C.adjugate() * u)[0])
assert sp.factor(R_num - (s+c-s*tau)) == 0

# Pure algebra: eliminate the scalar contraction from the reduced 5D flow.
S,TAU,C2,D = sp.symbols('s tau c D', nonzero=True)
DE1 = 4 + 3*S - 2*TAU + C2 - S*TAU
DE3 = (1+S)*(4+S-2*TAU) + C2
E4 = 1+S
D2R = S + C2 - S*TAU
assert sp.factor(DE1*(1+E4) - 2*DE3 - (E4-1)*D2R) == 0

E1,E2,E3,E4s = sp.symbols('E1 E2 E3 E4')
DR_from_E = sp.factor((E1*(1+E4s)-2*E3)/(E4s-1))
Y_reduced = sp.factor((E1-4-DR_from_E)/4)
Y_expected = sp.factor((E1-E3+2*E4s-2)/(2*(1-E4s)))
assert sp.factor(Y_reduced-Y_expected) == 0

# Convert to elementary symmetric polynomials of T itself.
lam = sp.symbols('lambda', nonzero=True)
t1,t2,t3,t4 = sp.symbols('e1T e2T e3T e4T')
ME1 = 4-lam*t1
ME3 = 4-3*lam*t1+2*lam**2*t2-lam**3*t3
ME4 = 1-lam*t1+lam**2*t2-lam**3*t3+lam**4*t4
Y_T = sp.factor(Y_expected.subs({E1:ME1,E3:ME3,E4s:ME4}))
Q_T = sp.factor((2*lam*t4-t3)/(2*(t1-lam*t2+lam**2*t3-lam**3*t4)))
assert sp.simplify(Y_T-lam**2*Q_T) == 0

# Null branch discriminant.
alpha,beta = sp.symbols('alpha beta', real=True)
null_E1 = alpha+beta
null_E2 = 2+alpha*beta
null_disc = sp.factor(null_E1**2 - 4*null_E2 + 8)
assert sp.factor(null_disc-(alpha-beta)**2) == 0


def exact_example(vals_u, vals_f, label):
    sub = {
        u0: sp.Rational(vals_u[0]), u1: sp.Rational(vals_u[1]),
        u2: sp.Rational(vals_u[2]), u3: sp.Rational(vals_u[3]),
        f01: sp.Rational(vals_f[0]), f02: sp.Rational(vals_f[1]),
        f03: sp.Rational(vals_f[2]), f12: sp.Rational(vals_f[3]),
        f13: sp.Rational(vals_f[4]), f23: sp.Rational(vals_f[5]),
    }
    Cn = C.subs(sub); Cdn = Cdag.subs(sub); C5n = C5.subs(sub)
    C5dn = g5 * C5n.T * g5
    D2n = sp.factor(Cn.det()); Dn = sp.sqrt(D2n)
    M4 = sp.simplify(Dn*sp.Rational(1,2)*(Cn.inv()+Cdn.inv()))
    M5 = sp.simplify(Dn*sp.Rational(1,2)*(C5n.inv()+C5dn.inv()))
    coeff = M4.charpoly().all_coeffs()
    E1n,E2n,E3n,E4n = map(sp.factor,(-coeff[1],coeff[2],-coeff[3],coeff[4]))
    sn = sp.factor(s.subs(sub)); taun = sp.factor(tau.subs(sub)); cn = sp.factor(c.subs(sub))
    Rn = sp.factor((uf*C.adjugate()*u)[0].subs(sub)/D2n)
    Ydirect = sp.factor((2+sn-taun)/(2*Dn)-1)
    Y5 = sp.factor((sp.trace(M5)-4-M5.det())/4)
    assert sp.simplify(Y5-Ydirect) == 0
    assert sp.simplify(M5.det()-Dn) == 0
    assert sp.simplify(sp.trace(M5) - (E1n + Dn*(1-Rn))) == 0
    if sn != 0:
        Y4 = sp.factor((E1n-E3n+2*E4n-2)/(2*(1-E4n)))
    elif any(sp.Rational(x) != 0 for x in vals_u):
        disc = sp.factor(E1n**2-4*E2n+8)
        Y4 = sp.factor((E1n+sp.sqrt(disc))/4-1)
    else:
        Y4 = sp.factor(E1n/4-1)
    assert sp.simplify(Y4-Ydirect) == 0
    return {
        'sector': label, 'u_squared': str(sn), 'D_squared': str(D2n),
        '5d_flow_minus_direct_4d_flow': 'exact zero',
        'det_I_minus_lambda_T5_minus_D': 'exact zero',
        'block_trace_relation': 'exact zero',
        '4d_stress_only_minus_direct_flow': 'exact zero',
    }

examples = [
    exact_example((sp.Rational(1,4),0,0,0),(sp.Rational(1,5),sp.Rational(1,7),0,sp.Rational(1,6),0,sp.Rational(1,8)),'timelike'),
    exact_example((0,sp.Rational(1,4),sp.Rational(1,5),0),(sp.Rational(1,6),sp.Rational(1,9),sp.Rational(1,10),sp.Rational(1,7),sp.Rational(1,11),sp.Rational(1,8)),'spacelike'),
    exact_example((sp.Rational(1,4),sp.Rational(1,4),0,0),(sp.Rational(1,5),sp.Rational(1,8),sp.Rational(1,9),sp.Rational(1,7),sp.Rational(1,11),sp.Rational(1,6)),'null'),
    exact_example((0,0,0,0),(sp.Rational(1,5),sp.Rational(1,8),sp.Rational(1,9),sp.Rational(1,7),sp.Rational(1,11),sp.Rational(1,6)),'pure Born-Infeld'),
]

# Independent float64 random validation.
g4n=np.diag([-1.,1.,1.,1.]); g5n=np.diag([-1.,1.,1.,1.,1.]); I4n=np.eye(4); I5n=np.eye(5)

def f_matrix(vals):
    f01,f02,f03,f12,f13,f23=vals
    return np.array([[0,f01,f02,f03],[-f01,0,f12,f13],[-f02,-f12,0,f23],[-f03,-f13,-f23,0]],float)

def elementary_4(M):
    tr1=np.trace(M); tr2=np.trace(M@M); tr3=np.trace(M@M@M)
    return tr1,0.5*(tr1*tr1-tr2),(tr1**3-3*tr1*tr2+2*tr3)/6.0,np.linalg.det(M)

def one_numeric(p_contra,fvals,lam):
    p=np.asarray(p_contra,float); pf=p@g4n; F=f_matrix(fvals); K0=g4n@F
    K=np.sqrt(lam)*K0; u=np.sqrt(lam)*p; uf=u@g4n
    C=I4n+K+np.outer(u,uf); Cdag=I4n-K+np.outer(u,uf)
    detC=np.linalg.det(C)
    if detC<=1e-8: return None
    D=np.sqrt(detC); M4=0.5*D*(np.linalg.inv(C)+np.linalg.inv(Cdag))
    C5=np.zeros((5,5)); C5[:4,:4]=I4n+K; C5[:4,4]=u; C5[4,:4]=-uf; C5[4,4]=1.0
    C5dag=g5n@C5.T@g5n; M5=0.5*D*(np.linalg.inv(C5)+np.linalg.inv(C5dag))
    T5=(I5n-M5)/lam
    Q5=1/(4*lam**2)-np.trace(T5)/(4*lam)-np.linalg.det(M5)/(4*lam**2)
    Cprime=(0.5/np.sqrt(lam))*K0 + np.outer(p,p@g4n)
    Dprime=0.5*D*np.trace(np.linalg.solve(C,Cprime))
    Qdirect=(D-1)/lam**2-Dprime/lam
    E1n,E2n,E3n,E4n=elementary_4(M4); sphys=p@g4n@p
    if abs(sphys)>1e-9: Y=(E1n-E3n+2*E4n-2)/(2*(1-E4n))
    elif np.linalg.norm(p)>1e-12: Y=(E1n+np.sqrt(max(E1n*E1n-4*E2n+8,0.0)))/4-1
    else: Y=E1n/4-1
    return Q5,Qdirect,Y/lam**2,np.linalg.det(M5)-D

rng=np.random.default_rng(20260820); rows=[]
for sector in ['timelike','spacelike','null','pure Born-Infeld']:
  for lamv in [0.1,0.3,0.7]:
    errs5=[]; errs4=[]; deterrs=[]; n=0
    while n<1500:
      fvals=rng.normal(0,0.25,6)
      if sector=='timelike':
        spatial=rng.normal(0,0.12,3); p0=0.25+abs(rng.normal(0,0.08)); p=np.r_[p0,spatial]
        if p@g4n@p>=-1e-4: continue
      elif sector=='spacelike':
        p=np.r_[rng.normal(0,0.08),rng.normal(0,0.22,3)]
        if p@g4n@p<=1e-4: continue
      elif sector=='null':
        nv=rng.normal(size=3); nv/=np.linalg.norm(nv); a=0.15+0.12*rng.random(); p=np.r_[a,a*nv]
      else: p=np.zeros(4)
      out=one_numeric(p,fvals,lamv)
      if out is None: continue
      q5,qd,q4,de=out
      if not (np.isfinite(q5) and np.isfinite(qd) and np.isfinite(q4)): continue
      errs5.append(q5-qd); errs4.append(q4-qd); deterrs.append(de); n+=1
    rows.append({'sector':sector,'lambda':lamv,'samples':n,
      'rms_5d_reduction_vs_direct':float(np.sqrt(np.mean(np.square(errs5)))),
      'max_5d_reduction_vs_direct':float(np.max(np.abs(errs5))),
      'rms_4d_stress_formula_vs_direct':float(np.sqrt(np.mean(np.square(errs4)))),
      'max_4d_stress_formula_vs_direct':float(np.max(np.abs(errs4))),
      'max_abs_detM5_minus_D':float(np.max(np.abs(deterrs)))})

pd.DataFrame(rows).to_csv(CSV,index=False)
summary={'generic_5d_determinant_reduction':'exact zero','generic_D2_R_identity':'exact zero',
'elimination_of_R_without_lambda_derivative':'exact zero','reduced_5d_flow_equals_generic_4d_stress_only_formula':'exact zero',
'conversion_to_T_invariants':'exact zero','null_branch_discriminant':'exact square',
'exact_rational_matrix_examples':examples,'numerical_validation_csv':str(CSV)}
JSON.write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda x:f'{x:.3e}'))
