"""Rational reconstruction identities proving persistence near the exact point."""
import sympy as S,json,time,platform
from pathlib import Path
start=time.perf_counter();r,s,p,q=S.symbols('rho s p q');t=S.symbols('t');zero={};rows={}
def check(name,x):
 v=S.factor(x);assert v==0,(name,v);zero[name]=str(v)
for name,c,o,sign in [('A',r,s,-1),('B',s,r,1)]:
 l=(c-q)/(c*q);D=1-l*p;R=D*(1+l*o);X=(R-1)/l
 E=(1+R-D*(2+l*(c-q)))/l
 B=(R-l*E-D*D)/(l*l*E);A=sign*(X-B)
 P=1+t*X-t*E-t*t*B*E
 check(name+'_X',X-(B+sign*A))
 check(name+'_determinant',P.subs(t,l)-D*D)
 check(name+'_block_trace',(1+R-l*E)/D-(2+l*(c-q)))
 check(name+'_block_determinant',P.subs(t,l)/D**2-1)
 check(name+'_reciprocal',(1+l*c)*(1-l*q)-1)
 check(name+'_remaining_axis',R/D-(1+l*o))
 check(name+'_L',(1-D)/l-p)
 Q=(D-1)/l**2-S.diff(P,t).subs(t,l)/(2*l*D)
 check(name+'_direct_Q',Q-c*q/2)
 vals=[l,D,R,X,E,A,B];base={r:3,s:2,p:S.Rational(1,2),q:S.Rational(7,10)}
 at=[S.factor(x.subs(base)) for x in vals];assert all(v>0 for v in at)
 rows[name]=dict(zip(['lambda','Delta','R','X','electric_squared','a_squared','b_squared'],map(str,at)))
# Six pair-products distinct implies both selected roots of K are simple.
eigs=[-r,-s,p,q];K=S.prod(S.Symbol('Q')+eigs[i]*eigs[j]/2 for i in range(4) for j in range(i+1,4));derivatives={}
for name,root in [('A',r*q/2),('B',s*q/2)]:
 der=S.factor(S.diff(K,S.Symbol('Q')).subs(S.Symbol('Q'),root).subs(base));assert der!=0;derivatives[name]=str(der)
result={'source_commit':'1b546752fb96dc1ee65a3b7bb1827110b1d2eb20','python':platform.python_version(),'sympy':S.__version__,'status':'all_passed','identities':zero,'exact_zero_count':len(zero),'positive_base_values':rows,'K_derivatives':derivatives,'elapsed_seconds':round(time.perf_counter()-start,3),'conclusion':'Rational reconstruction and all strict physical conditions persist on an open neighborhood of (rho,s,p,q)=(3,2,1/2,7/10). Both X>0. Distinct eigenvalues and nonzero pair-product-root derivatives persist; local branches differ.'}
Path(__file__).with_name('reconstruction_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
