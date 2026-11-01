#!/usr/bin/env python3
"""Exact two-field null DBI check; run python check.py. SymPy only, no network."""
from pathlib import Path
import json, platform, time
import sympy as s
start=time.perf_counter(); checks={}
def zero(name,expression):
    residual=s.factor(s.simplify(expression))
    assert residual==0,(name,residual)
    checks[name]=str(residual)
l,a=s.symbols('lambda a',positive=True)
e,b,t=s.symbols('e b t',real=True)
g=s.diag(-1,1,1,1);v=s.Matrix([a,a,0,0]);F=s.zeros(4)
F[0,2]=e;F[2,0]=-e;F[1,2]=b;F[2,1]=-b
h=g+l*v*v.T;G=h+s.sqrt(l)*F
P=s.factor(G.det()/g.det());D=s.sqrt(P)
k=a**2*(e-b)**2
zero('determinant',P-(1+l*(b*b-e*e)-l*l*k))
Gi=G.inv();Sym=(Gi+Gi.T)/2
T=(g-D*Sym)/l;Theta=(h.inv()-D*Sym)/l
L=(1-D)/l;Q=s.diff(L,l)
Tvv=s.factor((v.T*T*v)[0])
zero('directional_contraction',Tvv-k/D)
C=s.factor((v.T*(Theta-h.inv()*L)*v)[0]/2)
zero('metric_chain_correction',C-k/(2*D))
H=Theta*h;Qh=(s.trace(H*H)-s.trace(H)**2/2)/8
zero('CH_chain_rule',Q-Qh-C)
zero('trace_directional_identity',Q+s.trace(T*g)/(4*l)-Tvv/4)
M=s.eye(4)-l*T*g;E1=s.trace(M);E2=(E1**2-s.trace(M*M))/2
Zdisc=s.factor(E1**2-4*E2+8)
zero('discriminant_square',Zdisc-l**4*k**2/P)
# Positive branch P>0 and k>=0 fixes sqrt(discriminant)=lambda**2*k/D.
root=l*l*k/D
zero('null_invariant_flow',Q-(E1+root-4)/(4*l*l))
Fm=g*F;X=(v.T*g*v)[0];Z=(v.T*Fm*Fm*g*v)[0]
w=F*g*v
zero('scalar_norm',X)
zero('Z_equals_minus_w_squared',Z+(w.T*g*w)[0])
zero('twofield_Z',Z+k)
Gt=g+l*t*v*v.T+s.sqrt(l)*F
Pt=s.factor(Gt.det()/g.det())
zero('amplitude_determinant',Pt-(1+l*(b*b-e*e)-l*l*t*k))
zero('amplitude_Jacobi',s.diff(Pt,t).subs(t,1)-P*l*(v.T*Gi*v)[0])
zero('general_null_contraction_specialized',Tvv+Z/D)
# Test the general-null bridge at the level of known invariant identities.
Ds,SS,ZZ=s.symbols('Delta S Z',real=True,nonzero=True)
alpha=(2-2*l*SS+l*l*ZZ)/Ds; beta=(2-2*l*SS)/Ds
E1_known=alpha+beta;E2_known=2+alpha*beta
zero('known_factorization_discriminant',E1_known**2-4*E2_known+8-l**4*ZZ**2/Ds**2)
Qchain_known=(E1_known-4)/(4*l*l)-ZZ/(4*Ds)
zero('general_null_CH_bridge',Qchain_known-(E1_known-l*l*ZZ/Ds-4)/(4*l*l))
zero('transverse_limit',P.subs(b,0)-(1-l*e**2-l*l*a*a*e*e))
zero('coincident_e_b_determinant',P.subs(b,e)-1)
zero('coincident_e_b_Q',s.simplify(Q.subs(b,e)))
zero('coincident_e_b_correction',s.simplify(C.subs(b,e)))
zero('opposite_e_b_determinant',P.subs(b,-e)-(1-4*l*l*a*a*e*e))
zero('pure_BI_limit',s.simplify((Q-Qh).subs(a,0)))
zero('weak_coupling_limit',s.limit(Q,l,0,dir='+')-((b*b-e*e)**2/8+k/2))
Fa=s.zeros(4);Fa[0,1]=e;Fa[1,0]=-e
Ga=h+s.sqrt(l)*Fa;Pa=s.factor(Ga.det()/g.det());Da=s.sqrt(Pa)
Gai=Ga.inv();Ta=(g-Da*(Gai+Gai.T)/2)/l
zero('aligned_determinant',Pa-(1-l*e*e))
zero('aligned_contraction',(v.T*Ta*v)[0])
# Canonical nonzero null vector with all six independent field components.
f01,f02,f03,f12,f13,f23,r=s.symbols('f01 f02 f03 f12 f13 f23 r',real=True)
F6=s.Matrix([[0,f01,f02,f03],[-f01,0,f12,f13],[-f02,-f12,0,f23],[-f03,-f13,-f23,0]])
w6=F6*g*v
Z6=-(w6.T*g*w6)[0]
zero('sixfield_null_Z',Z6+a*a*((f02-f12)**2+(f03-f13)**2))
P6=(g+r*r*t*v*v.T+r*F6).det(method='domain-ge')/g.det()
zero('sixfield_amplitude_derivative',s.diff(P6,t)-r**4*Z6)
zero('sixfield_amplitude_linearity',s.diff(P6,t,2))
result={'status':'all_passed','source_commit':'c304d28526e7ed1b8deff80890eb1f7ac060f712','manuscript_blob_sha1':'4e78787d298956360ec2971bc06a1c5a973adc9a','python':platform.python_version(),'sympy':s.__version__,'platform':platform.platform(),'branch_assumptions':'lambda>0,a>0,e,b real; P=1+lambda*(b^2-e^2)-lambda^2*a^2*(e-b)^2>0; D=positive sqrt(P)','exact_residuals':checks,'formulas':{'P':str(P),'k':str(k),'Tvv':str(k/D),'metric_correction':str(k/(2*D)),'discriminant':str(l**4*k*k/P),'weak_limit':str((b*b-e*e)**2/8+k/2)},'scope':'Two-field symbolic verification and six-field polynomial amplitude checks; general-null CH bridge reuses known manuscript determinant and characteristic-polynomial identities, not an independent new matrix theorem.','elapsed_seconds':round(time.perf_counter()-start,3)}
Path(__file__).with_name('results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
