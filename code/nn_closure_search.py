#!/usr/bin/env python3
"""4D single-scalar DBI stress-tensor flow discovery/validation.

Conventions
-----------
Flat metric eta = diag(-1,+1,+1,+1).  With
    A_mu = sqrt(lambda) d_mu phi,
    B_mn = sqrt(lambda) F_mn,
write
    G_mn = eta_mn + A_m A_n + B_mn,
    D = sqrt(det G / det eta) = sqrt(-det G).
The mixed stress tensor satisfies
    M = I - lambda T = D G^{-1 (mn)} eta_{n r}.
The dimensionless flow target is
    y = lambda^2 * (1/sqrt(-g)) d L_lambda / d lambda
at fixed physical fields.

The program:
  * generates generic random field configurations;
  * trains a neural network O_theta(E1,E2,E3,E4) as a scout;
  * performs rational SVD readout of the neural target;
  * validates an exact closed-form candidate on independent samples;
  * checks an algebraic identity among EMT invariants;
  * checks the Jacobian rank in a canonical Lorentz frame;
  * compares the lambda->0 limit with the p=3 perturbative result of
    Blair et al., JHEP 02 (2025) 024, arXiv:2410.03591.

The exact candidate found by the readout is
    y = (E1 - E3 + 2 E4 - 2) / (2 (1-E4)),
where Ek = e_k(M).
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import time
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

np.set_printoptions(precision=8, suppress=True)
torch.set_default_dtype(torch.float32)

OUT = Path('/mnt/data')
CSV_SUMMARY = OUT / 'dbi_4d_single_scalar_validation.csv'
CSV_LAMBDA = OUT / 'dbi_4d_single_scalar_error_by_lambda.csv'
FIG_SCATTER = OUT / 'dbi_4d_nn_scout_scatter.png'
FIG_LAMBDA = OUT / 'dbi_4d_exact_formula_error_vs_lambda.png'
FIG_SV = OUT / 'dbi_4d_jacobian_singular_values.png'
JSON_READOUT = OUT / 'dbi_4d_rational_readout.json'

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
I4 = np.eye(4)


def elementary_symmetric_4(M: np.ndarray) -> np.ndarray:
    """Return e1,e2,e3,e4 from traces/determinant for a 4x4 matrix."""
    tr1 = np.trace(M)
    M2 = M @ M
    tr2 = np.trace(M2)
    tr3 = np.trace(M2 @ M)
    e1 = tr1
    e2 = 0.5 * (tr1 * tr1 - tr2)
    e3 = (tr1**3 - 3.0 * tr1 * tr2 + 2.0 * tr3) / 6.0
    e4 = np.linalg.det(M)
    return np.array([e1, e2, e3, e4], dtype=float)


def antisymmetric_from_six(x: np.ndarray) -> np.ndarray:
    """Six independent components -> covariant antisymmetric 4x4 matrix."""
    B = np.zeros((4, 4), dtype=float)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for val, (i, j) in zip(x, pairs):
        B[i, j] = val
        B[j, i] = -val
    return B


def dbi_dimensionless(A: np.ndarray, B: np.ndarray):
    """Compute D, mixed M=I-lambda T, invariants E_k, and y=lambda^2 Q.

    A=sqrt(lambda) q and B=sqrt(lambda)F are dimensionless.
    Returns None outside the real DBI branch.
    """
    G = ETA + np.outer(A, A) + B
    detG = np.linalg.det(G)
    if not np.isfinite(detG) or detG >= -1e-10:
        return None
    D = math.sqrt(-detG)
    if D < 1e-8:
        return None
    try:
        Ginv = np.linalg.inv(G)
    except np.linalg.LinAlgError:
        return None
    Ginv_sym = 0.5 * (Ginv + Ginv.T)
    M = D * Ginv_sym @ ETA
    E = elementary_symmetric_4(M)

    # lambda D' = D/2 Tr[G^{-1}(A A^T + B/2)] at fixed q,F.
    lambda_Dprime = 0.5 * D * np.trace(Ginv @ (np.outer(A, A) + 0.5 * B))
    y = D - 1.0 - lambda_Dprime
    return D, M, E, y


def exact_y_from_E(E: np.ndarray) -> float:
    e1, e2, e3, e4 = E
    return (e1 - e3 + 2.0 * e4 - 2.0) / (2.0 * (1.0 - e4))


def invariant_identity_residual(E: np.ndarray) -> float:
    e1, e2, e3, e4 = E
    lhs = (e4 - 1.0)**2 * e2
    rhs = (-e1*e1*e4 + e1*e3*(e4+1.0) - e3*e3
           + e4**3 - e4**2 - e4 + 1.0)
    return lhs - rhs


def sample_dimensionless(n: int, seed: int, a_scale=0.42, b_scale=0.38,
                         d_min=0.18, e4_gap=2e-3):
    rng = np.random.default_rng(seed)
    rows = []
    attempts = 0
    while len(rows) < n and attempts < 80*n:
        attempts += 1
        A = rng.normal(scale=a_scale, size=4)
        six = rng.normal(scale=b_scale, size=6)
        B = antisymmetric_from_six(six)
        out = dbi_dimensionless(A, B)
        if out is None:
            continue
        D, M, E, y = out
        if D < d_min or abs(1.0-E[3]) < e4_gap:
            continue
        if not np.all(np.isfinite(E)) or not np.isfinite(y):
            continue
        # Avoid extreme outliers that make the scout NN needlessly difficult.
        if np.max(np.abs(E)) > 20 or abs(y) > 5:
            continue
        rows.append((A, B, D, E, y))
    if len(rows) < n:
        raise RuntimeError(f'Could only generate {len(rows)} / {n} samples')
    return rows


def canonical_B(e, p, b):
    B = np.zeros((4,4), dtype=float)
    B[0,1], B[1,0] = e, -e
    B[2,3], B[3,2] = p, -p
    B[3,1], B[1,3] = b, -b
    return B


def sample_canonical(n: int, seed: int):
    rng = np.random.default_rng(seed)
    X, y = [], []
    configs = []
    while len(X) < n:
        a = rng.uniform(0.18, 0.72)
        e = rng.uniform(-0.42, 0.42)
        p = rng.uniform(-0.52, 0.52)
        b = rng.uniform(-0.52, 0.52)
        A = np.array([a,0,0,0.0])
        B = canonical_B(e,p,b)
        out = dbi_dimensionless(A,B)
        if out is None: continue
        D,M,E,target = out
        if D < 0.30 or abs(1-E[3]) < 5e-3 or abs(target)>2:
            continue
        X.append(E)
        y.append(target)
        configs.append([a,e,p,b])
    return np.asarray(X), np.asarray(y), np.asarray(configs)


class ScoutNet(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(4, 64), torch.nn.Tanh(),
            torch.nn.Linear(64, 64), torch.nn.Tanh(),
            torch.nn.Linear(64, 64), torch.nn.Tanh(),
            torch.nn.Linear(64, 1),
        )
    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_neural_scout(seed=42, n=18000, epochs=260):
    E, y, _ = sample_canonical(n, seed=seed)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    ntr = int(0.8*n)
    tr, te = idx[:ntr], idx[ntr:]

    mean = E[tr].mean(axis=0)
    std = E[tr].std(axis=0) + 1e-8
    ymean = y[tr].mean()
    ystd = y[tr].std() + 1e-8

    Xtr = torch.tensor((E[tr]-mean)/std, dtype=torch.float32)
    Ytr = torch.tensor((y[tr]-ymean)/ystd, dtype=torch.float32)
    Xte = torch.tensor((E[te]-mean)/std, dtype=torch.float32)
    Yte_raw = y[te]

    torch.manual_seed(seed)
    model = ScoutNet()
    opt = torch.optim.Adam(model.parameters(), lr=2e-3)
    batch = 512
    for ep in range(epochs):
        perm = torch.randperm(ntr)
        for s in range(0, ntr, batch):
            bidx = perm[s:s+batch]
            pred = model(Xtr[bidx])
            loss = torch.mean((pred-Ytr[bidx])**2)
            opt.zero_grad(); loss.backward(); opt.step()

    with torch.no_grad():
        pred = (model(Xte).numpy()*ystd + ymean)
    err = pred - Yte_raw
    metrics = {
        'nn_test_rms': float(np.sqrt(np.mean(err**2))),
        'nn_test_max_abs': float(np.max(np.abs(err))),
        'nn_test_relative_l2': float(np.linalg.norm(err)/np.linalg.norm(Yte_raw)),
        'nn_test_n': int(len(te)),
    }
    return model, E[te], Yte_raw, pred, metrics


def rational_svd_readout(E: np.ndarray, y: np.ndarray):
    """Fit y=(a.f)/(b.f), f=(1,E1,E2,E3,E4), with homogeneous SVD."""
    f = np.column_stack([np.ones(len(E)), E])
    A = np.column_stack([f, -y[:,None]*f])
    U, S, Vh = np.linalg.svd(A, full_matrices=False)
    vec = Vh[-1]
    a = vec[:5]
    b = vec[5:]
    # Normalize denominator constant to +2 when possible.
    if abs(b[0]) > 1e-12:
        scale = 2.0/b[0]
    else:
        scale = 1.0/np.max(np.abs(vec))
    a *= scale; b *= scale
    pred = (f@a)/(f@b)
    return a,b,S, float(np.sqrt(np.mean((pred-y)**2)))


def arbitrary_validation(n=30000, seed=20260820):
    rows = sample_dimensionless(n, seed)
    abs_errors=[]; id_errors=[]; ys=[]
    for A,B,D,E,y in rows:
        yp = exact_y_from_E(E)
        abs_errors.append(abs(yp-y))
        id_errors.append(abs(invariant_identity_residual(E)))
        ys.append(y)
    ae=np.asarray(abs_errors); ie=np.asarray(id_errors); ys=np.asarray(ys)
    return {
        'n': n,
        'formula_rms_abs': float(np.sqrt(np.mean(ae**2))),
        'formula_max_abs': float(ae.max()),
        'formula_median_abs': float(np.median(ae)),
        'formula_p99_abs': float(np.quantile(ae,0.99)),
        'formula_max_normalized': float(np.max(ae/(1+np.abs(ys)))),
        'identity_rms_abs': float(np.sqrt(np.mean(ie**2))),
        'identity_max_abs': float(ie.max()),
    }


def physical_lambda_validation(n_per_lambda=2500, seed=17):
    rng=np.random.default_rng(seed)
    lambdas=[0.02,0.05,0.1,0.3,0.7,1.2]
    records=[]
    for lam in lambdas:
        errors=[]; norm_errors=[]; kept=0; attempts=0
        # Physical q,F held as randomly sampled and then A,B scale with sqrt(lambda).
        # Scales are chosen so all lambda bins have a healthy real-branch sample.
        while kept<n_per_lambda and attempts<100*n_per_lambda:
            attempts+=1
            q=rng.normal(scale=0.36,size=4)
            F=antisymmetric_from_six(rng.normal(scale=0.31,size=6))
            A=math.sqrt(lam)*q; B=math.sqrt(lam)*F
            out=dbi_dimensionless(A,B)
            if out is None: continue
            D,M,E,y=out
            if D<0.18 or abs(1-E[3])<2e-4 or abs(y)>5: continue
            yp=exact_y_from_E(E)
            # Q error = |y-y_pred|/lambda^2. Report dimensionless y error to avoid
            # artificial lambda^-2 amplification of floating point noise.
            e=abs(yp-y)
            errors.append(e); norm_errors.append(e/(1+abs(y))); kept+=1
        e=np.asarray(errors); ne=np.asarray(norm_errors)
        records.append({
            'lambda':lam,'n':kept,
            'rms_abs_y_error':float(np.sqrt(np.mean(e**2))),
            'median_abs_y_error':float(np.median(e)),
            'p99_abs_y_error':float(np.quantile(e,0.99)),
            'max_abs_y_error':float(np.max(e)),
            'max_normalized_y_error':float(np.max(ne)),
        })
    return pd.DataFrame(records)


def canonical_map(params, include_y=False):
    a,e,p,b=params
    out=dbi_dimensionless(np.array([a,0,0,0.]),canonical_B(e,p,b))
    if out is None: return None
    D,M,E,y=out
    return np.r_[E,y] if include_y else E


def finite_difference_jacobian(params, include_y=False, h=2e-6):
    params=np.asarray(params,dtype=float)
    f0=canonical_map(params,include_y)
    J=np.zeros((len(f0),4))
    for j in range(4):
        pp=params.copy(); pm=params.copy(); pp[j]+=h; pm[j]-=h
        fp=canonical_map(pp,include_y); fm=canonical_map(pm,include_y)
        J[:,j]=(fp-fm)/(2*h)
    return J


def rank_check():
    points=np.array([
        [0.40,0.25,0.35,0.20],
        [0.60,0.30,0.15,0.40],
        [0.30,-0.20,0.50,-0.25],
        [0.52,-0.31,-0.22,0.36],
        [0.45,0.18,-0.42,-0.28],
    ])
    rows=[]
    sv_E=[]; sv_Ey=[]
    for i,p in enumerate(points):
        s=np.linalg.svd(finite_difference_jacobian(p,False),compute_uv=False)
        sy=np.linalg.svd(finite_difference_jacobian(p,True),compute_uv=False)
        sv_E.append(s); sv_Ey.append(sy)
        rec={'point':i}
        for k,val in enumerate(s,1): rec[f'sv_E_{k}']=val
        for k,val in enumerate(sy,1): rec[f'sv_Ey_{k}']=val
        rows.append(rec)
    return pd.DataFrame(rows),np.asarray(sv_E),np.asarray(sv_Ey)


def sigma_from_T(T):
    return elementary_symmetric_4(T)


def free_stress_mixed(q: np.ndarray, Fcov: np.ndarray):
    q2=float(q@ETA@q)
    Tsc_cov=np.outer(q,q)-0.5*ETA*q2
    Fup=ETA@Fcov@ETA
    F2scalar=float(np.sum(Fcov*Fup))
    F_nu_up=Fcov@ETA
    first=Fcov@F_nu_up.T
    Tem_cov=first-0.25*ETA*F2scalar
    Tcov=Tsc_cov+Tem_cov
    Tmix=ETA@Tcov
    return Tmix


def perturbative_check(n=10000,seed=91):
    rng=np.random.default_rng(seed)
    diffs=[]
    for _ in range(n):
        q=rng.normal(scale=0.65,size=4)
        F=antisymmetric_from_six(rng.normal(scale=0.55,size=6))
        T=free_stress_mixed(q,F)
        s1,s2,s3,s4=sigma_from_T(T)
        if abs(s1)<2e-5: continue
        candidate=-s3/(2*s1)
        trT2=np.trace(T@T); trT=np.trace(T)
        Fmix=ETA@F
        trF2=np.trace(Fmix@Fmix)
        trF4=np.trace(Fmix@Fmix@Fmix@Fmix)
        blair=0.25*trT2-0.125*trT**2-0.125*(trF4-0.25*trF2**2)
        diffs.append(candidate-blair)
    d=np.asarray(diffs)
    return {
        'n':int(len(d)),
        'rms_abs_difference':float(np.sqrt(np.mean(d**2))),
        'max_abs_difference':float(np.max(np.abs(d))),
        'median_abs_difference':float(np.median(np.abs(d))),
    }


def main():
    t0=time.time()
    print('Training neural scout...')
    model,Etest,ytest,ypred,nn_metrics=train_neural_scout()
    print(nn_metrics)

    # Rational readout uses exact targets on an independent canonical set.
    Eread,yread,_=sample_canonical(12000,seed=707)
    a,b,S,read_rms=rational_svd_readout(Eread,yread)
    readout={
        'numerator_coefficients_[1,E1,E2,E3,E4]':a.tolist(),
        'denominator_coefficients_[1,E1,E2,E3,E4]':b.tolist(),
        'smallest_singular_value':float(S[-1]),
        'second_smallest_singular_value':float(S[-2]),
        'readout_rms':read_rms,
    }
    JSON_READOUT.write_text(json.dumps(readout,indent=2),encoding='utf-8')
    print('Rational readout:',readout)

    print('Arbitrary covariant validation...')
    validation=arbitrary_validation()
    print(validation)

    print('Lambda-binned validation...')
    lambda_df=physical_lambda_validation()
    lambda_df.to_csv(CSV_LAMBDA,index=False)
    print(lambda_df)

    print('Jacobian rank check...')
    rank_df,svE,svEy=rank_check()
    print(rank_df)

    print('Perturbative p=3 check...')
    pert=perturbative_check()
    print(pert)

    # Summary CSV as key-value rows plus rank info.
    summary_rows=[]
    for k,v in nn_metrics.items(): summary_rows.append({'category':'neural_scout','quantity':k,'value':v})
    for k,v in validation.items(): summary_rows.append({'category':'exact_validation','quantity':k,'value':v})
    for k,v in pert.items(): summary_rows.append({'category':'perturbative_check','quantity':k,'value':v})
    summary_rows.append({'category':'rational_readout','quantity':'rms','value':read_rms})
    summary_rows.append({'category':'rational_readout','quantity':'smallest_singular_value','value':float(S[-1])})
    summary_rows.append({'category':'rational_readout','quantity':'second_smallest_singular_value','value':float(S[-2])})
    for _,r in rank_df.iterrows():
        for k in range(1,5):
            summary_rows.append({'category':f'jacobian_point_{int(r.point)}','quantity':f'sv_E_{k}','value':r[f'sv_E_{k}']})
            summary_rows.append({'category':f'jacobian_point_{int(r.point)}','quantity':f'sv_Ey_{k}','value':r[f'sv_Ey_{k}']})
    pd.DataFrame(summary_rows).to_csv(CSV_SUMMARY,index=False)

    # Figure 1: NN scout scatter.
    plt.figure(figsize=(6.4,5.2))
    take=np.linspace(0,len(ytest)-1,min(2500,len(ytest))).astype(int)
    plt.scatter(ytest[take],ypred[take],s=8,alpha=0.45)
    lo=min(ytest[take].min(),ypred[take].min()); hi=max(ytest[take].max(),ypred[take].max())
    plt.plot([lo,hi],[lo,hi])
    plt.xlabel(r'direct DBI target $y=\lambda^2 Q$')
    plt.ylabel(r'neural prediction $O_\theta(E_1,E_2,E_3,E_4)$')
    plt.title('Neural scout on held-out single-scalar DBI configurations')
    plt.tight_layout(); plt.savefig(FIG_SCATTER,dpi=180); plt.close()

    # Figure 2: exact candidate error by lambda.
    plt.figure(figsize=(6.6,4.8))
    plt.semilogy(lambda_df['lambda'],lambda_df['median_abs_y_error'],marker='o',label='median')
    plt.semilogy(lambda_df['lambda'],lambda_df['p99_abs_y_error'],marker='s',label='99th percentile')
    plt.semilogy(lambda_df['lambda'],lambda_df['max_abs_y_error'],marker='^',label='maximum')
    plt.xlabel(r'deformation parameter $\lambda$')
    plt.ylabel(r'absolute error in $y=\lambda^2 Q$')
    plt.title('Exact stress-tensor formula: independent numerical validation')
    plt.legend(); plt.tight_layout(); plt.savefig(FIG_LAMBDA,dpi=180); plt.close()

    # Figure 3: Jacobian singular values.
    plt.figure(figsize=(6.6,4.8))
    x=np.arange(1,5)
    for i,row in enumerate(svE):
        plt.semilogy(x,row,marker='o',alpha=0.65,label=f'point {i+1}')
    plt.xlabel('singular-value index')
    plt.ylabel('singular value')
    plt.xticks(x)
    plt.title(r'Jacobian of $(a,e,p,b)\mapsto(E_1,E_2,E_3,E_4)$')
    plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(FIG_SV,dpi=180); plt.close()

    print(f'Runtime {time.time()-t0:.1f}s')
    print('Saved files:')
    for p in [CSV_SUMMARY,CSV_LAMBDA,FIG_SCATTER,FIG_LAMBDA,FIG_SV,JSON_READOUT]: print(p)

if __name__=='__main__':
    main()
