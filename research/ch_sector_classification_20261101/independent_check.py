import sympy as s
r=s.Rational
t,z=s.symbols('t z',real=True)
g=s.diag(-1,1,1,1)
I=s.eye(4)
data=[('A',r(23,21),s.Matrix([s.sqrt(r(361,238)),0,s.sqrt(r(98,51)),0]),s.sqrt(r(51,140))),('B',r(13,14),s.Matrix([0,s.sqrt(r(1665,2212)),s.sqrt(r(28,79)),0]),s.sqrt(r(79,56)))]
for name,lam,v,e in data:
 F=s.zeros(4);F[0,1]=e;F[1,0]=-e
 h=g+t*v*v.T;G=h+s.sqrt(t)*F
 P=s.factor(G.det()/g.det());P0=P.subs(t,lam)
 Delta=s.sqrt(P0);L=(1-Delta)/lam
 Q=s.simplify(-(1-Delta)/lam**2-s.diff(P,t).subs(t,lam)/(2*lam*Delta))
 Gi=s.simplify(G.subs(t,lam).inv());symGi=(Gi+Gi.T)/2
 Tcontra=s.simplify((g.inv()-Delta*symGi)/lam);T=s.simplify(Tcontra*g)
 spec=T.eigenvects()
 print('\nCANDIDATE',name,'lambda',lam)
 print('v²=',s.simplify((v.T*g.inv()*v)[0]),'X=1+lambda v²=',s.simplify(1+lam*(v.T*g.inv()*v)[0]))
 print('P(t)=',P,'; endpoints',P.subs(t,0),P0,'; Psecond=',s.diff(P,t,2))
 print('Delta,L,Q=',Delta,L,Q,'\nT=',T)
 print('spectrum and eigenvector g-norms',[(val,mult,[s.simplify((w.T*g*w)[0]) for w in ws]) for val,mult,ws in spec])
 # Direct variation of action density sqrt(-detg)*L, independently differentiating determinants.
 checks=[]
 for a in range(4):
  for b in range(a,4):
   E=s.zeros(4);E[a,b]=1;E[b,a]=1
   gz=g+z*E
   Gz=gz+lam*v*v.T+s.sqrt(lam)*F
   Pz=s.factor(Gz.det()/gz.det())
   # stress = 2/sqrt(-g) derivative density; off-diagonal z varies both entries
   deriv=s.diff(s.sqrt(-gz.det())*(1-s.sqrt(Pz))/lam,z).subs(z,0)
   expected=Tcontra[a,b]*(1 if a!=b else r(1,2))
   checks.append(s.simplify(deriv-expected)==0)
 print('All 10 independent metric variations PASS=',all(checks))
 assert all(checks)
 assert set(T.eigenvals())=={s.Integer(-3),s.Integer(-2),r(1,2),r(7,10)}
 assert L==r(1,2)
 assert Q==({'A':r(21,20),'B':r(7,10)}[name])
