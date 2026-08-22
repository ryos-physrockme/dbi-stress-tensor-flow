#!/usr/bin/env python3
"""Independent numerical validation of the frame-independent 4D single-scalar DBI flow.

The script uses float64 matrix algebra in Minkowski signature (-,+,+,+). It tests:
  * det(I-lambda T) = 1 + lambda (partial phi)^2;
  * the generic stress-tensor-only flow for non-null scalar gradient;
  * the regular null-gradient branch;
  * the pure 4D Born-Infeld limit;
  * the algebraic relation among elementary symmetric polynomials of I-lambda T.

Dimensionless variables u=sqrt(lambda) d phi and B=sqrt(lambda) F are sampled directly.
"""
from __future__ import annotations
import os
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('DBI_OUTPUT_DIR', REPO_ROOT / 'data')).resolve()
OUT.mkdir(parents=True, exist_ok=True)
CSV=OUT/'coordinate_free_validation.csv'
FIG=OUT/'coordinate_free_validation_errors.png' 
ETA=np.diag([-1.0,1.0,1.0,1.0])
I4=np.eye(4)


def anti6(x):
    B=np.zeros((4,4),dtype=np.float64)
    pairs=((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
    for val,(i,j) in zip(x,pairs):
        B[i,j]=val; B[j,i]=-val
    return B


def esym4(M):
    p1=np.trace(M)
    M2=M@M
    p2=np.trace(M2)
    p3=np.trace(M2@M)
    e1=p1
    e2=(p1*p1-p2)/2
    e3=(p1**3-3*p1*p2+2*p3)/6
    e4=np.linalg.det(M)
    return np.array([e1,e2,e3,e4],dtype=np.float64)


def dbi(u,B):
    G=ETA+np.outer(u,u)+B
    detG=np.linalg.det(G)
    if not np.isfinite(detG) or detG>=-1e-12:
        return None
    D=math.sqrt(-detG)
    if D<1e-10:
        return None
    Gi=np.linalg.inv(G)
    M=D*0.5*(Gi+Gi.T)@ETA
    E=esym4(M)
    # At fixed physical dphi,F: lambda D' = D/2 Tr[G^{-1}(u u^T + B/2)].
    lDp=0.5*D*np.trace(Gi@(np.outer(u,u)+0.5*B))
    y=D-1.0-lDp
    s=float(u@ETA@u)
    return D,M,E,y,s


def y_generic(E):
    e1,_,e3,e4=E
    return (e1-e3+2*e4-2)/(2*(1-e4))


def y_null(E):
    e1,e2,_,_=E
    disc=max(e1*e1-4*e2+8,0.0)
    return 0.25*(e1+math.sqrt(disc))-1.0


def y_pure_bi(E):
    return E[0]/4.0-1.0


def algebraic_constraint(E):
    e1,e2,e3,e4=E
    return ((e4-1)**2*e2
        -(-e1*e1*e4+e1*e3*(e4+1)-e3*e3+e4**3-e4**2-e4+1))


def sample_u(rng,sector):
    if sector=='timelike':
        spatial=rng.normal(scale=0.18,size=3)
        t=math.sqrt(float(spatial@spatial)+rng.uniform(0.05,0.35))
        return np.r_[t,spatial]
    if sector=='spacelike':
        t=rng.normal(scale=0.12)
        vec=rng.normal(size=3); vec/=np.linalg.norm(vec)
        r=math.sqrt(t*t+rng.uniform(0.05,0.45))
        return np.r_[t,r*vec]
    if sector=='null nonzero':
        vec=rng.normal(size=3); vec/=np.linalg.norm(vec)
        a=rng.uniform(0.08,0.55)
        return np.r_[a,a*vec]
    if sector=='pure Born-Infeld':
        return np.zeros(4)
    raise ValueError(sector)


def collect(sector,n,seed):
    rng=np.random.default_rng(seed)
    generic_err=[]; null_err=[]; bi_err=[]; det_err=[]; alg_err=[]; cvals=[]
    count=0; attempts=0
    while count<n and attempts<100*n:
        attempts+=1
        u=sample_u(rng,sector)
        B=anti6(rng.normal(scale=0.28,size=6))
        out=dbi(u,B)
        if out is None: continue
        D,M,E,y,s=out
        if D<0.18 or np.max(np.abs(E))>30 or abs(y)>10: continue
        # Enforce the intended causal sector numerically.
        if sector=='timelike' and not s<-1e-7: continue
        if sector=='spacelike' and not s>1e-7: continue
        if sector=='null nonzero' and abs(s)>2e-14: continue

        det_err.append(abs(E[3]-(1+s)))
        alg_err.append(abs(algebraic_constraint(E)))

        if sector in ('timelike','spacelike'):
            generic_err.append(abs(y_generic(E)-y))
        elif sector=='null nonzero':
            null_err.append(abs(y_null(E)-y))
            # c=u_mu (K^2)^mu_nu u^nu, K=g^{-1}B.
            K=ETA@B
            c=float((u@ETA)@(K@K)@u)
            cvals.append(c)
        else:
            # The direct pure-BI limit is numerically much better conditioned
            # than the square-root null formula at its double-root point.
            bi_err.append(abs(y_pure_bi(E)-y))
        count+=1
    if count<n:
        raise RuntimeError((sector,count,attempts))

    flow = generic_err or null_err or bi_err
    row={
        'sector':sector,
        'n':count,
        'flow_rms_abs_error':float(np.sqrt(np.mean(np.square(flow)))),
        'flow_max_abs_error':float(np.max(flow)),
        'flow_p99_abs_error':float(np.quantile(flow,0.99)),
        'det_identity_max_abs_error':float(np.max(det_err)),
        'algebraic_constraint_max_abs_error':float(np.max(alg_err)),
    }
    if cvals:
        row['max_c=uK2u']=float(np.max(cvals))
        row['fraction_c_positive']=float(np.mean(np.asarray(cvals)>1e-12))
    else:
        row['max_c=uK2u']=np.nan
        row['fraction_c_positive']=np.nan
    return row, np.asarray(flow)


if __name__=='__main__':
    specs=[
        ('timelike',10000,202608201),
        ('spacelike',10000,202608202),
        ('null nonzero',15000,202608203),
        ('pure Born-Infeld',10000,202608204),
    ]
    rows=[]; errors={}
    for sector,n,seed in specs:
        row,err=collect(sector,n,seed)
        rows.append(row); errors[sector]=err
    df=pd.DataFrame(rows)
    df.to_csv(CSV,index=False)
    print(df.to_string(index=False,float_format=lambda x:f'{x:.6e}'))

    labels=[x[0] for x in specs]
    med=[np.median(errors[s]) for s in labels]
    p99=[np.quantile(errors[s],0.99) for s in labels]
    mx=[np.max(errors[s]) for s in labels]
    x=np.arange(len(labels))
    plt.figure(figsize=(7.4,4.8))
    plt.semilogy(x,med,marker='o',label='median')
    plt.semilogy(x,p99,marker='s',label='99th percentile')
    plt.semilogy(x,mx,marker='^',label='maximum')
    plt.xticks(x,labels,rotation=12)
    plt.ylabel('absolute error in dimensionless flow')
    plt.title('Independent validation of the 4D DBI stress-tensor flow')
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG,dpi=180)
    print(f'Wrote {CSV}')
    print(f'Wrote {FIG}')
