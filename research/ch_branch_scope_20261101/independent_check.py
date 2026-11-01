"""Independent axial DBI review: differentiate metric and lambda BEFORE substitution."""
import sympy as S
from pathlib import Path
q,u,w=S.symbols('q u w',positive=True);rho=q+u+w;s=q+u
l,a,b=S.symbols('lambda a b',positive=True);d=S.symbols('d',positive=True)
gx,gy=S.symbols('g0 g1',nonzero=True)
g=S.diag(gx,gy,1,1);flat={gx:-1,gy:1}
checks=0;out=[]
def zero(x):
 global checks
 assert S.factor(S.cancel(x))==0,S.factor(x)
 checks+=1
for name,index,ls,aa,bb,DD,Q in [('A',0,(rho-q)/(rho*q),q*(rho-s)/rho,q*(q+s)/rho,q/rho,rho*q/2),('B',1,(s-q)/(s*q),q*(rho-s)/s,q*(rho+q)/s,q/s,s*q/2)]:
 v=S.zeros(4,1);v[index]=S.sqrt(a);F=S.zeros(4);F[0,1]=S.sqrt(b);F[1,0]=-S.sqrt(b)
 G=g+l*v*v.T+S.sqrt(l)*F;P=S.factor(G.det()/g.det());Pf=S.factor(P.subs(flat))
 # Derivatives evaluated before field/target substitution. Square root replaced by d at endpoint.
 L=(1-S.sqrt(P))/l
 rawT=[S.simplify((2*g[i,i]*S.diff(L,z)+L).subs(flat)) for i,z in [(0,gx),(1,gy)]]
 rawQ=S.diff((1-S.sqrt(Pf))/l,l)
 sub={a:aa,b:bb,l:ls}
 # square roots simplify automatically in positive q,u,w variables
 for expr,expected in zip(rawT,[-rho,-s]):zero(S.simplify(expr.subs(sub))-expected)
 zero((1-DD)/ls-q);zero(S.simplify(rawQ.subs(sub))-Q)
 zero(Pf.subs(sub)-DD**2)
 # Full tensor independently from matrix inverse, hence off-diagonals checked too.
 eta=g.subs(flat);Gi=G.subs(flat).subs(sub).inv()
 TT=(eta-DD*(Gi+Gi.T)/2)*eta/ls
 for x in TT-S.diag(-rho,-s,q,q):zero(x)
 R=S.factor(((g+l*v*v.T).det()/g.det()).subs(flat).subs({a:aa}))
 pp=S.factor(Pf.subs({a:aa,b:bb}))
 assert S.degree(R,l)<=1 and S.degree(pp,l)<=1
 for z in [R.subs(l,0),R.subs(l,ls),pp.subs(l,0),pp.subs(l,ls)]:assert S.factor(z).is_positive
 eps=S.symbols('epsilon',positive=True)
 limits=[S.limit(expr.subs({u:eps,w:eps}),eps,0,dir='+') for expr in [ls,aa,bb,Q,DD]]
 assert limits==[0,0,2*q,q*q/2,1]
 out.append(f'{name}: lambda={S.factor(ls)}, v_squared={S.factor(aa)}, electric_squared={S.factor(bb)}, Q={S.factor(Q)}, limits={limits}')
zero(rho*q/2-s*q/2-q*w/2)
out.append(f'PASS: {checks} exact residuals; direct diagonal metric variation, full T, fixed-field lambda derivative, interval positivity, and limits.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(out)+'\n');print('\n'.join(out))
