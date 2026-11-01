#!/usr/bin/env python3
"""Bounded exact algebra for lambda elimination; no external data needed."""
import sympy as s
import json,time,platform
from pathlib import Path
start=time.perf_counter();checks={}
def zero(name,x):
 r=s.factor(s.simplify(x)); assert r==0,(name,r);checks[name]=str(r)
l=s.symbols('lambda');e1,e2,e3,e4=s.symbols('e1 e2 e3 e4');es=(e1,e2,e3,e4)
E1=4-l*e1;E2=6-3*l*e1+l*l*e2;E3=4-3*l*e1+2*l*l*e2-l**3*e3;E4=1-l*e1+l*l*e2-l**3*e3+l**4*e4
C=(E4-1)**2*E2+E1**2*E4-E1*E3*(E4+1)+E3**2-E4**3+E4**2+E4-1
H=s.cancel(C/l**6)
assert s.Poly(H,l).degree()==6
xs=s.symbols('t1:5');sub={es[i]:sum(s.prod(c) for c in __import__('itertools').combinations(xs,i+1)) for i in range(4)}
pairs=[(i,j) for i in range(4) for j in range(i+1,4)]
product=s.prod(l*xs[i]*xs[j]-xs[i]-xs[j] for i,j in pairs)
zero('constraint_six_pair_factorization',H.subs(sub)+product)
u,v=s.symbols('u v',nonzero=True)
lam_pair=(u+v)/(u*v)
zero('reciprocal_pair',(1-lam_pair*u)*(1-lam_pair*v)-1)
zero('flow_pair',(2-lam_pair*(u+v)-2)/(2*lam_pair**2)+u*v/2)
Q=s.symbols('Q')
PQ=s.prod(Q+xs[i]*xs[j]/2 for i,j in pairs)
formal,remainder,mapping=s.symmetrize(s.expand(PQ),xs,formal=True)
zero('flow_polynomial_symmetry',remainder)
K=s.expand(formal.subs({mapping[i][0]:es[i] for i in range(4)}))
zero('flow_polynomial_pair_product',K.subs(sub)-PQ)
# Degenerate BI spectrum: cross-pair lambda roots coincide, yet the same Q remains.
a,b=s.symbols('a b',nonzero=True)
zero('BI_repeated_constraint',product.subs(dict(zip(xs,[a,a,b,b])))-(l*a*a-2*a)*(l*b*b-2*b)*(l*a*b-a-b)**4)
zero('BI_repeated_flow',PQ.subs(dict(zip(xs,[a,a,b,b])))-(Q+a*a/2)*(Q+b*b/2)*(Q+a*b/2)**4)
# A generic exact physical DBI point establishes a nonempty simple-root region.
g=s.diag(-1,1,1,1);vv=s.Matrix([s.Rational(1,3),s.Rational(1,4),0,s.Rational(1,5)])
F=s.zeros(4)
for (i,j),value in {(0,1):s.Rational(1,5),(0,2):s.Rational(1,7),(1,2):s.Rational(1,6),(2,3):s.Rational(1,4)}.items():F[i,j]=value;F[j,i]=-value
h=g+vv*vv.T;G=h+F;P=s.factor(G.det()/g.det());r=s.factor(h.det()/g.det());assert P>0 and r>0
D=s.sqrt(P);Gi=G.inv();T=(g-D*(Gi+Gi.T)/2)*g
tr1=s.trace(T);tr2=s.trace(T*T);tr3=s.trace(T*T*T)
vals=[tr1,(tr1**2-tr2)/2,(tr1**3-3*tr1*tr2+2*tr3)/6,T.det()]
vals=[s.simplify(z) for z in vals];point=dict(zip(es,vals))
zero('physical_point_constraint',H.subs(point).subs(l,1))
Hp=s.factor(s.diff(H,l).subs(point).subs(l,1));assert Hp!=0
Qdirect=(D-1)-D*s.trace(Gi*(vv*vv.T+F/2))/2
zero('physical_point_autonomous_polynomial',K.subs(point).subs(Q,Qdirect))
Qp=s.factor(s.diff(K,Q).subs(point).subs(Q,Qdirect));assert Qp!=0
# Two real, regular DBI configurations with the same full mixed stress tensor.
target=s.diag(-2,-1,s.Rational(1,2),s.Rational(1,2))
cp={e1:s.trace(target),e2:sum(target[i,i]*target[j,j] for i,j in pairs),e3:sum(s.prod(target[i,i] for i in c) for c in __import__('itertools').combinations(range(4),3)),e4:target.det()}
counterexample={}
for name,lstar,vec,ef,expected_D,expected_Q in [
 ('timelike',s.Rational(3,2),[s.Rational(1,2),0,0,0],s.sqrt(s.Rational(3,8)),s.Rational(1,4),s.Rational(1,2)),
 ('spacelike',s.Integer(1),[0,1/s.sqrt(2),0,0],s.sqrt(5)/2,s.Rational(1,2),s.Rational(1,4))]:
 cv=s.Matrix(vec);cf=s.zeros(4);cf[0,1]=ef;cf[1,0]=-ef
 ch=g+l*cv*cv.T;cg=ch+s.sqrt(l)*cf
 cr=s.factor(ch.det()/g.det());cP=s.factor(cg.det()/g.det());cD=s.sqrt(cP)
 cGi=cg.inv();cT=s.simplify(((g-cD*(cGi+cGi.T)/2)/l*g).subs(l,lstar))
 cL=(1-cD)/l;cQ=s.simplify(s.diff(cL,l).subs(l,lstar))
 diff=cT-target
 zero(name+'_same_full_T',s.trace(diff.T*diff))
 zero(name+'_direct_Q',cQ-expected_Q)
 zero(name+'_same_L',s.simplify(cL.subs(l,lstar))-s.Rational(1,2))
 zero(name+'_Delta',cD.subs(l,lstar)-expected_D)
 zero(name+'_autonomous_polynomial',K.subs(cp).subs(Q,cQ))
 # r and P are affine in lambda; positive at both interval endpoints proves
 # positivity on the entire fixed-field interval 0 <= lambda <= lstar.
 assert s.degree(cr,l)<=1 and s.degree(cP,l)<=1
 assert cr.subs(l,0)>0 and cr.subs(l,lstar)>0
 assert cP.subs(l,0)>0 and cP.subs(l,lstar)>0
 counterexample[name]={'lambda':str(lstar),'v_cov':str(cv.T),'F01':str(ef),'r_lambda':str(cr),'Delta_squared_lambda':str(cP),'T_mixed':str(cT),'L':str(s.simplify(cL.subs(l,lstar))),'Q':str(cQ),'physical_interval_verified':True}
assert s.Rational(counterexample['timelike']['Q'])!=s.Rational(counterexample['spacelike']['Q'])
result={'status':'all_passed','source_commit':'c304d28526e7ed1b8deff80890eb1f7ac060f712','python':platform.python_version(),'sympy':s.__version__,'exact_residuals':checks,'H_lambda':str(s.collect(s.expand(H),l)),'K_Q':str(s.collect(K,Q)),'physical_point':{'lambda':'1','h_det_ratio':str(r),'G_det_ratio':str(P),'H_derivative_nonzero':str(Hp),'H_derivative_numeric':str(s.N(Hp,15)),'K_derivative_nonzero':str(Qp),'K_derivative_numeric':str(s.N(Qp,15)),'Q_direct':str(s.simplify(Qdirect))},'counterexample':counterexample,'conclusion':'Branch-local autonomous Q exists at simple roots. Two physical configurations with identical full mixed T and L but different Q rule out a single global T-only autonomous function covering both timelike and spacelike scalar sectors. No no-go for a restricted sector is claimed.','elapsed_seconds':round(time.perf_counter()-start,3)}
Path(__file__).with_name('results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,indent=2))
