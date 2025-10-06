import numpy as np
import galois

def gf(p: int):
    return galois.GF(p)

def random_poly(secret, t, GF):
    #degree t-1 polynomial with a0 = secret
    coeffs = GF.Random(t) # length t random in GF
    coeffs[0] = GF(secret)
    return galois.Poly(coeffs, field=GF) #a0 + a1x + ...

def make_shares(secret: int, n: int, t:int , p:int):
    """
    Returns x,y in GF(p) with x=1..n (nonzero, distinct)
    """

    if not (2 <= t <= n):
        raise ValueError("Need 2 <= trheshold <= n")
    
    GF = gf(p)
    f = random_poly(secret, t, GF)
    xs = GF.Range(1, n+1)
    ys = f(xs)
    return GF, f, xs, ys

def lagrange_interpolate_at_zero(xs, ys, GF):
    """
    Reconstruct f(0) via Lagrange interpolation at x=0.
    Words for any number of points; correct only if >= threshold points from same polynomial.
    """

    x0 = GF(0)
    s = GF(0)

    for i in range(len(xs)):
        num = GF(1)
        den = GF(1)
        for j in range(len(xs)):
            if i ==j:
                continue
            num *= (x0 - xs[j])
            den *= (xs[i]- xs[j])
        Li0 = num / den
        s += ys[i] * Li0
    return int(s)