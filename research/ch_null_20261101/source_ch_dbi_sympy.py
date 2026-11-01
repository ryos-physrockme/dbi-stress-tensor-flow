#!/usr/bin/env python3
"""Exact checks of the Courant--Hilbert construction for four-dimensional DBI.

Conventions:
  metric signature (-,+,+,+); lambda > 0;
  h_{mu nu} = g_{mu nu} + lambda v_mu v_nu, v_mu = partial_mu phi;
  G_{mu nu} = h_{mu nu} + sqrt(lambda) F_{mu nu};
  L_DBI = (1 - sqrt(det(G)/det(g))) / lambda.

The script tests the ordinary Born--Infeld CH solution and the metric-chain-rule
identity relating its h-metric stress tensor to the full DBI stress tensor.
All examples use exact rational input and SymPy algebraic numbers.
The physical branch requires det(h)/det(g) > 0 and det(G)/det(g) > 0.

Reference:
  Pezhman, Babaei-Aghbolagh and Babaei Velni, arXiv:2610.23281v1.
"""

import sympy as sp


def assert_zero(expression, description):
    residual = sp.simplify(sp.factor(expression.rewrite(sp.exp)))
    if residual != 0:
        raise AssertionError(f"{description}: {residual}")


def check_courant_hilbert():
    """CH differential equation, BI solution, flow and auxiliary dictionary."""
    U, V, lam = sp.symbols("U V lam", real=True, positive=True)
    R = sp.sqrt((1 + 2 * lam * U) * (1 - 2 * lam * V))
    L = (1 - R) / lam
    assert_zero(sp.diff(L, U) * sp.diff(L, V) + 1, "CH PDE")

    s = sp.symbols("s", positive=True)
    a = (1 - s**2) / (2 * s)
    b = -(1 - s)**2 / (2 * s)
    flow_coefficient = -b / (4 * (a**2 + b**2))
    assert_zero(
        flow_coefficient - s / (4 * (1 + s**2)),
        "BI CH parametric flow coefficient",
    )
    r = (1 - s) / (1 + s)
    assert_zero(
        flow_coefficient - (1 - r**2) / (8 * (1 + r**2)),
        "BI CH autonomous flow coefficient",
    )

    # Russo--Townsend auxiliary-field dictionary for the same BI family.
    varphi = sp.symbols("varphi", positive=True, real=True)
    tau = (1 - sp.exp(-2 * varphi)) / (2 * lam)
    ell = (1 - sp.exp(-varphi)) / lam
    ell_prime = sp.exp(varphi)
    W = (sp.cosh(varphi) - 1) / lam
    assert_zero(W - (tau * ell_prime - ell), "auxiliary potential sign")
    assert_zero(sp.diff(W, varphi) - tau * ell_prime, "auxiliary W prime")
    L_aux = -sp.exp(-varphi) * U + sp.exp(varphi) * V - W
    assert_zero(
        sp.diff(L_aux, varphi).subs(V, tau - U * sp.exp(-2 * varphi)),
        "auxiliary equation of motion",
    )
    assert_zero(
        L_aux.subs(V, tau - U * sp.exp(-2 * varphi))
        - (ell - 2 * U / ell_prime),
        "auxiliary on-shell action",
    )
    print("PASS: CH equation, BI flow, and auxiliary-field dictionary")


def check_pure_bi_symbolically():
    """The BI CH quadratic flow and trace flow for parallel E and B."""
    lam, electric, magnetic = sp.symbols(
        "lam electric magnetic", positive=True, real=True
    )
    g = sp.diag(-1, 1, 1, 1)
    F = sp.zeros(4)
    F[0, 1], F[1, 0] = electric, -electric
    F[2, 3], F[3, 2] = magnetic, -magnetic
    G = g + sp.sqrt(lam) * F
    D = sp.sqrt(sp.factor(G.det() / g.det()))
    L = (1 - D) / lam
    stress_contra = (g.inv() - D * (G.inv() + G.inv().T) / 2) / lam
    stress_mixed = stress_contra * g
    trace_t = sp.trace(stress_mixed)
    trace_t_squared = sp.trace(stress_mixed * stress_mixed)
    flow_direct = sp.diff(L, lam)
    flow_ch = (trace_t_squared - trace_t**2 / 2) / 8
    assert_zero(flow_direct - flow_ch, "BI quadratic stress flow")
    assert_zero(flow_direct + trace_t / (4 * lam), "BI trace flow")
    print("PASS: symbolic pure BI quadratic and trace flows")


def make_field_strength(entries):
    F = sp.zeros(4)
    for (row, col), value in entries.items():
        F[row, col] = value
        F[col, row] = -value
    return F


def elementary_symmetric_coefficients(matrix):
    t1 = sp.trace(matrix)
    t2 = sp.trace(matrix * matrix)
    t3 = sp.trace(matrix * matrix * matrix)
    return (
        t1,
        (t1**2 - t2) / 2,
        (t1**3 - 3 * t1 * t2 + 2 * t3) / 6,
        matrix.det(),
    )


def check_dbi_configuration(name, vector, entries, lam=sp.Rational(1, 2)):
    """Compare independent determinant, CH-chain-rule, and old DBI formulas."""
    g = sp.diag(-1, 1, 1, 1)
    v = sp.Matrix(vector)
    F = make_field_strength(entries)
    scalar_norm = (v.T * g.inv() * v)[0]
    h = g + lam * v * v.T
    G = h + sp.sqrt(lam) * F
    det_h = sp.factor(h.det() / g.det())
    det_G = sp.factor(G.det() / g.det())
    if not (det_h > 0 and det_G > 0):
        raise ValueError(f"{name}: outside positive determinant branch")
    q = sp.sqrt(det_h)
    B = sp.sqrt(sp.factor(G.det() / h.det()))
    Delta = sp.sqrt(det_G)
    assert_zero(Delta - q * B, f"{name}: determinant factorization")

    L_bi = (1 - B) / lam
    stress_dbi = (
        g.inv() - Delta * (G.inv() + G.inv().T) / 2
    ) / lam
    stress_bi = (
        h.inv() - B * (G.inv() + G.inv().T) / 2
    ) / lam
    stress_from_induced_metric = g.inv() / lam + q * (
        stress_bi - h.inv() / lam
    )
    for component in stress_dbi - stress_from_induced_metric:
        assert_zero(component, f"{name}: metric stress relation")

    # Jacobi's determinant formula at fixed (g, F, v).
    Gprime = v * v.T + F / (2 * sp.sqrt(lam))
    Q_direct = (
        (Delta - 1) / lam**2
        - Delta * sp.trace(G.inv() * Gprime) / (2 * lam)
    )
    Theta = stress_bi * h
    Q_bi_ch = (
        sp.trace(Theta * Theta) - sp.trace(Theta)**2 / 2
    ) / 8
    # Independent Jacobi formula at fixed induced metric h.
    Q_bi_fixed_h = (
        (B - 1) / lam**2
        - B * sp.trace(G.inv() * F) / (4 * lam * sp.sqrt(lam))
    )
    assert_zero(Q_bi_fixed_h - Q_bi_ch, f"{name}: BI CH identity")

    qprime = scalar_norm / (2 * q)
    Q_chain = (
        (q - 1) / lam**2
        + qprime * (L_bi - 1 / lam)
        + q * Q_bi_ch
        + q * (v.T * (stress_bi - h.inv() * L_bi) * v)[0] / 2
    )
    assert_zero(Q_direct - Q_chain, f"{name}: DBI chain-rule flow")

    Tmixed = stress_dbi * g
    e1, e2, e3, e4 = elementary_symmetric_coefficients(Tmixed)
    denominator = e1 - lam * e2 + lam**2 * e3 - lam**3 * e4
    assert_zero(denominator + scalar_norm, f"{name}: denominator identity")
    if scalar_norm != 0:
        Q_original = (2 * lam * e4 - e3) / (2 * denominator)
        assert_zero(Q_direct - Q_original, f"{name}: original DBI flow")
    else:
        M = sp.eye(4) - lam * Tmixed
        E1, E2, _, _ = elementary_symmetric_coefficients(M)
        Q_null = (
            E1 + sp.sqrt(E1**2 - 4 * E2 + 8) - 4
        ) / (4 * lam**2)
        assert_zero(Q_direct - Q_null, f"{name}: null-gradient flow")
    print(f"PASS: {name}")


def main():
    q = sp.Rational
    check_courant_hilbert()
    check_pure_bi_symbolically()
    cases = {
        "timelike scalar gradient with gauge field": (
            [q(1, 3), q(1, 4), 0, q(1, 5)],
            {(0, 1): q(1, 5), (0, 2): q(1, 7),
             (1, 2): q(1, 6), (2, 3): q(1, 4)},
        ),
        "spacelike scalar gradient with gauge field": (
            [q(1, 7), q(1, 3), q(1, 4), 0],
            {(0, 1): q(1, 5), (0, 2): q(1, 7),
             (1, 2): q(1, 6), (2, 3): q(1, 4)},
        ),
        "nonzero null scalar gradient with gauge field": (
            [q(1, 3), q(1, 3), 0, 0],
            {(0, 2): q(1, 5), (1, 2): q(1, 7), (2, 3): q(1, 4)},
        ),
        "pure scalar DBI": ([q(1, 3), q(1, 4), 0, q(1, 5)], {}),
        "pure Born--Infeld": (
            [0, 0, 0, 0],
            {(0, 1): q(1, 5), (2, 3): q(1, 4)},
        ),
    }
    for name, (vector, entries) in cases.items():
        check_dbi_configuration(name, vector, entries)
    print("All Courant--Hilbert / DBI exact checks passed.")


if __name__ == "__main__":
    main()
