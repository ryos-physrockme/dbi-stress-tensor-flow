import sympy as s,json,time,platform
from pathlib import Path
start=time.perf_counter(); checks={};l=s.symbols('lambda',positive=True);g=s.diag(-1,1,1,1)
target=s.diag(-3,-2,s.Rational(1,2),s.Rational(7,10));z=s.symbols('z')
def zero(k,x):
 r=s.simplify(s.radsimp(x));assert r==0,(k,r);checks[k]=str(r)
cases=[('A',s.Rational(23,21),s.Matrix([s.sqrt(s.Rational(361,238)),0,s.sqrt(s.Rational(98,51)),0]),s.sqrt(s.Rational(51,140)),s.Rational(19,42),s.Rational(21,20)),('B',s.Rational(13,14),s.Matrix([0,s.sqrt(s.Rational(1665,2212)),s.sqrt(s.Rational(28,79)),0]),s.sqrt(s.Rational(79,56)),s.Rational(15,28),s.Rational(7,10))]
output={}
for name,lv,v,e,D,Qexp in cases:
 F=s.zeros(4);F[0,1]=e;F[1,0]=-e
 h=g+l*v*v.T;G=h+s.sqrt(l)*F;P=s.factor(-G.det());R=s.factor(-h.det());X=(v.T*g*v)[0]
 zero(name+'_Delta',P.subs(l,lv)-D*D)
 Gi=G.subs(l,lv).inv();T=s.simplify((g-D*(Gi+Gi.T)/2)*g/lv)
 zero(name+'_characteristic',T.charpoly(z).as_expr()-target.charpoly(z).as_expr())
 Q=s.simplify((D-1)/lv**2-s.diff(P,l).subs(l,lv)/(2*lv*D))
 zero(name+'_Q',Q-Qexp);zero(name+'_L',(1-D)/lv-s.Rational(1,2))
 assert s.simplify(X)>0 and R.subs(l,lv)>0 and s.diff(P,l,2)<0 and P.subs(l,0)>0 and P.subs(l,lv)>0
 # Distinct real eigenvectors; normalize in Minkowski metric and put in target order.
 columns=[]
 for i,t in enumerate(target.diagonal()):
  vs=(T-t*s.eye(4)).nullspace();assert len(vs)==1
  w=vs[0];norm=s.simplify((w.T*g*w)[0]);assert s.sign(norm)==g[i,i]
  columns.append(s.simplify(w/s.sqrt(abs(norm))))
 E=s.Matrix.hstack(*columns)
 for i in range(4):
  for j in range(4):
   zero(name+f'_Lorentz{i}{j}',(E.T*g*E-g)[i,j])
   zero(name+f'_diagonal{i}{j}',(T*E-E*target)[i,j])
 vp=s.simplify(E.T*v);Fp=s.simplify(E.T*F*E)
 output[name]={'lambda':str(lv),'v_cov':str(v.T),'F01':str(e),'X':str(s.simplify(X)),'P_path':str(P),'R_path':str(R),'Delta':str(D),'Q':str(Q),'T_original':str(T),'frame_columns':str(E),'v_diagonal_frame':str(vp.T),'F_diagonal_frame':str(Fp),'interval_positive_by_concavity':True}
result={'source_commit':'1b546752fb96dc1ee65a3b7bb1827110b1d2eb20','paper_blob':'4e78787d298956360ec2971bc06a1c5a973adc9a','status':'all_passed','python':platform.python_version(),'sympy':s.__version__,'cases':output,'target_T_mixed':str(target),'exact_residual_count':len(checks),'residuals':checks,'elapsed_seconds':round(time.perf_counter()-start,3),'scope':'Both X>0, distinct stress eigenvalues, same full T after verified Lorentz frame alignment, same L. Not same lambda.'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
