"""Direct Hilbert and fixed-field derivative check (SymPy, CPU-small).

All quantities use the same fixed reference units as note.tex. g varies while
covariant F and v are held fixed. For diagonal metrics,
T^i_i = 2 g_ii dL/dg_ii + L (no sum). The independent full-inverse formula
also checks all off-diagonal components. Q is dL/dlambda at fixed (g,F,v).
Run: python metric_variation_check.py. No network or external data required.
"""
import sympy as s
l=s.symbols('lambda',positive=True);a=s.symbols('g:4',nonzero=True)
g=s.diag(*a);flat={a[0]:-1,a[1]:1,a[2]:1,a[3]:1}
target=s.diag(-2,-1,s.Rational(1,2),s.Rational(1,2))
for name,lv,v,e,expected in [('A',s.Rational(3,2),s.Matrix([s.Rational(1,2),0,0,0]),s.sqrt(s.Rational(3,8)),s.Rational(1,2)),('B',s.Integer(1),s.Matrix([0,1/s.sqrt(2),0,0]),s.sqrt(5)/2,s.Rational(1,4))]:
 F=s.zeros(4);F[0,1]=e;F[1,0]=-e
 h=g+l*v*v.T;G=h+s.sqrt(l)*F;ratio=s.factor(G.det()/g.det());D=s.sqrt(ratio);L=(1-D)/l
 td=[s.simplify((2*a[i]*s.diff(L,a[i])+L).subs(flat).subs(l,lv)) for i in range(4)]
 Q=s.simplify(s.diff(L,l).subs(flat).subs(l,lv))
 gi=G.subs(flat).subs(l,lv).inv();eta=g.subs(flat)
 full=s.simplify((eta-D.subs(flat).subs(l,lv)*(gi+gi.T)/2)*eta/lv)
 assert full==target and td==list(target.diagonal()) and Q==expected
 print(name,{'lambda':str(lv),'Delta_squared_path':str(s.factor(ratio.subs(flat))),'h_det_ratio_path':str(s.factor((h.det()/g.det()).subs(flat))),'T':str(full),'Q':str(Q),'L':str(s.simplify(L.subs(flat).subs(l,lv)))})
print('PASS: same full mixed T, different directly differentiated Q.')
