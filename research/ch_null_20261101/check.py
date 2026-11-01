#!/usr/bin/env python3
"""CPU-only exact follow-up to PR5; python check.py writes results.json."""
import contextlib, hashlib, io, json, platform, runpy, time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
start=time.perf_counter()
source=HERE/'source_ch_dbi_sympy.py'
b=source.read_bytes()
blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
assert blob=='8f5056aa25bc51d10f03e1aea3d64df2ae435d37'
buf=io.StringIO()
with contextlib.redirect_stdout(buf): runpy.run_path(str(source),run_name='__main__')
checks={}
def zero(name,expr):
    residual=s.factor(s.simplify(expr))
    assert residual==0,(name,residual)
    checks[name]=str(residual)
l,a,e=s.symbols('lambda a e',positive=True)
g=s.diag(-1,1,1,1);v=s.Matrix([a,a,0,0]);F=s.zeros(4)
F[0,2]=e;F[2,0]=-e
h=g+l*v*v.T;G=h+s.sqrt(l)*F
D2=s.factor(G.det()/g.det());D=s.sqrt(D2)
zero('null_norm',(v.T*g*v)[0]);zero('null_induced_determinant',h.det()/g.det()-1)
zero('determinant_polynomial',D2-(1-l*e**2-l**2*a**2*e**2))
Gi=G.inv();S=(Gi+Gi.T)/2
T=(g-D*S)/l;Theta=(h.inv()-D*S)/l
L=(1-D)/l;Q=s.diff(L,l)
H=Theta*h
Qbi=(s.trace(H*H)-s.trace(H)**2/2)/8
C=(v.T*(Theta-h.inv()*L)*v)[0]/2
zero('null_chain_rule',Q-Qbi-C)
zero('nonzero_metric_correction',C-a**2*e**2/(2*D))
M=s.eye(4)-l*T*g
E1=s.trace(M);E2=(E1**2-s.trace(M*M))/2
Tvv=(v.T*T*v)[0]
zero('directional_contraction',Tvv-a**2*e**2/D)
zero('trace_plus_directional_flow',Q+s.trace(T*g)/(4*l)-Tvv/4)
Z=s.factor(E1**2-4*E2+8)
zero('null_discriminant',Z-l**4*a**4*e**4/D2)
# On lambda>0,a>0,e>0,D2>0, the positive square root is l**2*a**2*e**2/D.
root=l**2*a**2*e**2/D
zero('null_invariant_flow_positive_branch',Q-(E1+root-4)/(4*l**2))
zero('pure_BI_limit',s.simplify((Q-Qbi).subs(a,0)))
# Independent pure-scalar limit with X real and q=sqrt(1+lambda X)>0.
X=s.symbols('X',real=True);q=s.sqrt(1+l*X)
Qs=s.diff((1-q)/l,l)
zero('pure_scalar_chain',Qs-((q-1)/l**2-X/(2*l*q)))
zero('pure_scalar_rationalized',Qs-X**2/(2*q*(q+1)**2))
# Generic chain-rule cancellation expressed through independent contractions A,B,C.
Delta,q0,A,B,C0=s.symbols('Delta q A B C',nonzero=True)
Lb=(1-Delta/q0)/l
Qfixed=(Delta/q0-1)/l**2-Delta/q0*C0/(4*l*s.sqrt(l))
Qchain=(q0-1)/l**2+(q0*A/2)*(Lb-1/l)+q0*Qfixed+Delta*(A-B)/(2*l)
Qdirect=(Delta-1)/l**2-Delta*B/(2*l)-Delta*C0/(4*l*s.sqrt(l))
zero('generic_chain_cancellation',Qchain-Qdirect)
traceT=(4-Delta*(4-l*B-s.sqrt(l)*C0))/l
zero('generic_trace_directional_identity',Qdirect+traceT/(4*l)+Delta*B/(4*l))
zero('weak_coupling_limit',s.limit(Q,l,0,dir='+')-(a**2*e**2/2+e**4/8))
# Aligned electric field is a control: a nonzero null gradient need not give a correction.
Fa=s.zeros(4);Fa[0,1]=e;Fa[1,0]=-e
Ga=h+s.sqrt(l)*Fa
Da2=s.factor(Ga.det()/g.det());Da=s.sqrt(Da2)
Ga_inv=Ga.inv();Ta=(g-Da*(Ga_inv+Ga_inv.T)/2)/l
zero('aligned_determinant',Da2-(1-l*e**2))
zero('aligned_directional_contraction',(v.T*Ta*v)[0])
Theta_a=(h.inv()-Da*(Ga_inv+Ga_inv.T)/2)/l
La=(1-Da)/l
zero('aligned_metric_correction',(v.T*(Theta_a-h.inv()*La)*v)[0]/2)
point={l:s.Rational(1,2),a:s.Rational(1,3),e:s.Rational(1,5)}
result={'status':'all_passed','source_commit':'62e0d609e6d5109675834b0d5f14e5afdfc47550','source_blob_sha1':blob,'python':platform.python_version(),'sympy':s.__version__,'platform':platform.platform(),'baseline_stdout':buf.getvalue(),'exact_residuals':checks,'null_family':{'det_ratio':str(D2),'metric_correction':str(a**2*e**2/(2*D)),'Tvv':str(a**2*e**2/D),'discriminant':str(Z),'branch':'lambda>0, a>0, e>0, 1-lambda*e^2-lambda^2*a^2*e^2>0','sample_metric_correction':str(s.simplify(C.subs(point)))},'limitations':['No new general original-EMT invariant closure proved','Null invariant formula already exists in source; this independently verifies and derives its root in a restricted family','Targeted primary-literature comparison is recorded in LITERATURE.md; no literature novelty claim'],'elapsed_seconds':round(time.perf_counter()-start,3)}
(HERE/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
