# logic/berlekamp_welch.py
import numpy as np
import galois

def berlekamp_welch_decode(xs, ys, k, e, p):
    """
    xs, ys: 1D arrays in GF(p)
    k: message degree bound (deg f < k)
    e: max errors to correct
    Returns: message polynomial f (galois.Poly), or raises if unsolvable.
    """
    GF = galois.GF(p)
    xs = GF(xs); ys = GF(ys)

    n = len(xs)
    if n < k + 2*e:
        raise ValueError("Need n >= k + 2e for unique decoding")

    # Unknowns: Q has degree < (k + e), E has degree <= e with leading coeff = 1 (to fix scale).
    q_deg = k + e - 1
    e_deg = e

    # Number of unknowns: (q_deg+1) + e_deg  (since leading coeff of E is set to 1)
    num_q = q_deg + 1
    num_e = e_deg

    # Build linear system A * u = b over GF
    # Unknown vector u = [q_0,...,q_qdeg, e_0,...,e_{e-1}]  (E's top coeff is 1)
    A = GF.Zeros((n, num_q + num_e))
    b = GF.Zeros(n)

    for i, (x, y) in enumerate(zip(xs, ys)):
        # Row encodes: sum_{j=0..q_deg} q_j x^j  - y * ( sum_{j=0..e_deg-1} e_j x^j  + x^{e_deg} )
        # Move known term y*x^{e_deg} to RHS.
        # Left side has unknown q_j and e_j.
        # Q(x) - y * (E(x)) = 0  with E leading coeff 1
        # => sum q_j x^j  - y * sum_{j=0..e_deg-1} e_j x^j  =  y * x^{e_deg}
        # Fill Q part
        pow_x = GF.Range(0, q_deg + 1)  # exponents
        A[i, :num_q] = x ** pow_x
        # Fill E part (note the minus y)
        if num_e > 0:
            pow_e = GF.Range(0, e_deg)   # 0..e-1
            A[i, num_q:] = -y * (x ** pow_e)
        # RHS
        b[i] = y * (x ** e_deg)

    # Solve
    u = galois.linalg.solve(A, b)  # vector in GF

    q_coeffs = u[:num_q]
    e_coeffs = np.concatenate([u[num_q:], GF([1])])  # add leading 1
    Q = galois.Poly(q_coeffs, field=GF)
    E = galois.Poly(e_coeffs, field=GF)

    # f = Q / E (polynomial division if divisible)
    f, r = Q // E, Q % E
    if r.degree >= 0 and any(r.coeffs != 0):
        raise ValueError("Decoding failed: Q not divisible by E")
    return f
