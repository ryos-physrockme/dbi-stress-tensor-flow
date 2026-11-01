#!/usr/bin/env python3
"""Exact axial DBI two-sector family; CPU-small, no search or external data."""
import json, platform, time
from pathlib import Path
import sympy as S
start=time.perf_counter(); checks={}
def zero(name, expression):
    r=S.factor(S.cancel(expression)); assert r==0,(name,r); checks[name]=str(r)
# Positive coordinates enforce rho > s > q > 0 without symbolic inequalities.
q,u,w=S.symbols('q u w',positive=True)
rho=q+u+w; s=q+u
lam=S.symbols('lambda',positive=True)
g=S.diag(-1,1,1,1)
target=S.diag(-rho,-s,q,q)
cases={
 'timelike':((rho-q)/(rho*q),q*(rho-s)/rho,q*(q+s)/rho,0,rho*q/2,q/rho),
 'spacelike':((s-q)/(s*q),q*(rho-s)/s,q*(rho+q)/s,1,s*q/2,q/s)}
results={}
for name,(ls,v2,e2,index,expectedQ,expectedD) in cases.items():
 v=S.zeros(4,1);v[index]=S.sqrt(v2)
 F=S.zeros(4);F[0,1]=S.sqrt(e2);F[1,0]=-F[0,1]
 h=g+lam*v*v.T; G=h+S.sqrt(lam)*F
 R=S.factor(h.det()/g.det()); P=S.factor(G.det()/g.det())
 zero(name+'_Delta_squared',P.subs(lam,ls)-expectedD**2)
 # expectedD>0, so this is the physical square-root branch at the target.
 Gi=G.subs(lam,ls).inv(); T=(g-expectedD*(Gi+Gi.T)/2)*g/ls
 for i in range(4):
  for j in range(4):zero(name+f'_T{i}{j}',T[i,j]-target[i,j])
 zero(name+'_L',(1-expectedD)/ls-q)
 Q=(expectedD-1)/ls**2-S.diff(P,lam).subs(lam,ls)/(2*ls*expectedD)
 zero(name+'_Q',Q-expectedQ)
 # P and R affine: positivity at endpoints proves full interval positivity.
 assert S.degree(P,lam)<=1 and S.degree(R,lam)<=1
 for expr in [R.subs(lam,0),R.subs(lam,ls),P.subs(lam,0),P.subs(lam,ls)]:
  assert S.factor(expr).is_positive
 results[name]={'lambda':str(S.factor(ls)),'v_squared':str(S.factor(v2)),
  'F01_squared':str(S.factor(e2)),'Delta':str(expectedD),'L':str(q),
  'Q':str(expectedQ),'h_det_ratio_path':str(R),'Delta_squared_path':str(P),
  'whole_interval_positive':True}
# Sector-specific reciprocal-pair equations select the unique nonzero lambda.
a,b=S.symbols('rho s',positive=True)
zero('timelike_reciprocal',(1+((a-q)/(a*q))*a)*(1-((a-q)/(a*q))*q)-1)
zero('spacelike_reciprocal',(1+((b-q)/(b*q))*b)*(1-((b-q)/(b*q))*q)-1)
zero('Q_difference',cases['timelike'][4]-cases['spacelike'][4]-q*w/2)
# u=w=epsilon gives a shared Maxwell limit, with both lambda -> 0+.
eps=S.symbols('epsilon',positive=True)
for name,(ls,v2,e2,_,Q,D) in cases.items():
 for label,expr,expected in [('lambda',ls,0),('scalar_amplitude_squared',v2,0),('electric_squared',e2,2*q),('Q',Q,q*q/2)]:
  zero(name+'_weak_'+label,S.limit(expr.subs({u:eps,w:eps}),eps,0,dir='+')-expected)
zero('weak_Q_difference',(cases['timelike'][4]-cases['spacelike'][4]).subs({u:eps,w:eps})-q*eps/2)
result={'source_commit':'45c11bc3f2396bcb523a6e4906bba89316e16ce9',
 'python':platform.python_version(),'sympy':S.__version__,'status':'all_passed',
 'domain':'q,u,w>0; s=q+u; rho=q+u+w; same reference units',
 'mixed_stress':'diag(-rho,-s,q,q)','cases':results,'exact_residuals':checks,
 'zero_residual_count':len(checks),'elapsed_seconds':round(time.perf_counter()-start,3),
 'scope':'Axial purely electric configurations only. Both causal sectors required; no fixed-sector or generic-simple-spectrum no-go.'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','zero_residual_count','elapsed_seconds']},indent=2))
