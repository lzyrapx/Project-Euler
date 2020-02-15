'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-08 17:03:15
'''
# coding:utf-8

"""
Functions using to solve various equations.
Functions list:
        linear_modulo_equation(a, b, n): Solve linear modular equation: (a * x) % n = b 
                                         or written in modular arithmatic: a * x = b (mod n)
                                         Return (smallest non-negative solution, step, n)
        square_modulo_prime_equation(n, p): Solve linear modular equation: (x^2) % p = n
                                            or written in modular arithmatic: x^2 = n (mod p)
                                            where n is a quadratic residue of p
                                            Return smallest solution a, noticing that (p-a) is also a solution.
                                            Using Tonelli-Shanks algorithm.
                                            Details see: http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
        square_modulo_prime_power_equation(n, p, k): Solve linear modular equation: (x^2) % p^k = n
                                                     or written in modular arithmatic: x^2 = n (mod p^k)
                                                     where n is a quadratic residue of p
                                                     Return all solutions betweem 0 and p^k-1, since there may exist more than 2 solutions.
                                                     Using Tonelli-Shanks algorithm and Hensel's Lift.
                                                     Details see:
                                                     - http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
                                                     - https://en.wikipedia.org/wiki/Hensel%27s_lemma
        chinese_remainder_theorem(equation_sets): return the only solution 0 <= x < (m_1 * m_2 * ... * m_r) for the set of simultaneous congruences:
                                                  x = a_i (mod m_i),  1 <= i <= r
                                                  where m_i are pairwise relatively prime
                                                  equation_sets must in the form [(a_1, m_1), (a2, m_2), ...]
                                                  details see: https://en.wikipedia.org/wiki/Chinese_remainder_theorem
        generalized_pell_equation_base(d, n=1):  Solve generalized Pell equation: x**2 - d * y**2 = n
                                                 Return smallest positive basic solution set, or enough solutions according to nsol
        generalized_pell_equation_generator(d, n=1): Solve generalized Pell equation: x**2 - d * y**2 = n
                                                     generate many solutions
        berlekamp_massey(sequence): Find minimum linear recurrence equation of sequence using Berlekamp-Massey algorithm.
                                    Return C = [c0, c1, ..., cL(=1)] with c0*si + c1*s(i+1) + ... + cL*s(i+L) = 0 for i+L >= n/2.
                                    @sequence (np.array)
        berlekamp_massey_with_bound(sequence, n): Find minimum linear recurrence equation of sequence using Berlekamp-Massey algorithm.
                                                  Return A = [a0, a1, ..., aL(=1)] with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0 for i+L >= n.
                                                  @sequence (np.array): sequence with at least 2*n terms
                                                  @n (int): dimension upper bound of the recurrence equation
        berlekamp_massey_mod_p(sequence, p): Find minimum linear recurrence equation of sequence in Z/Zp using Berlekamp-Massey algorithm.
                                             Return C = [c0, c1, ..., cL(=1)] with c0*si + c1*s(i+1) + ... + cL*s(i+L) = 0 for i+L >= len(seq)/2.
                                             @sequence (np.array)
                                             @p (int): prime number of field Z/Zp
        berlekamp_massey_with_bound_mod_p(sequence, n, p):  Find minimum linear recurrence equation of sequence in Z/Zp using Berlekamp-Massey algorithm.
                                                            Return A = [a0, a1, ..., aL(=1)] with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0 for i+L >= n.
                                                            @sequence (np.array): sequence with at least 2*n terms
                                                            @n (int): dimension upper bound of the recurrence equation
                                                            @p (int): prime number of field Z/Zp
"""

import numpy as np
from math import gcd, sqrt
from gmpy2 import invert
from base import cprod, is_square, legendre_symbol
from polynomial import poly_truncate, poly_add, poly_divmod, poly_mul

def linear_modulo_equation(a, b, n):
    """
    Solve linear modular equation
        (a * x) % n = b
    or written in modular arithmatic
         a * x = b (mod n)
    Return (smallest non-negative solution, step, n)
    """

    d = gcd(a, n)
    if b % d:
        raise ValueError('No Solution for ({} * x) % {} = {}!'.format(a, n, b))

    aa = a // d
    bb = b // d
    nn = n // d

    x0, x1, p, q = 0, 1, aa, nn
    while q != 1:
        p, q = q, p % q
        x0, x1 = x1, x0 - p // q * x1
    if (aa * x0) % nn != 1:
        x0 = -x0
    if x0 < 0:
        x0 = nn + x0
    sol = (x0 * bb) % n
    while sol > nn:
        sol = (sol + nn) % n
    return sol, nn, n

def square_modulo_prime_equation(n, p):
    """
    Solve linear modular equation
        (x^2) % p = n
    or written in modular arithmatic
         x^2 = n (mod p)
    where n is a quadratic residue of p
    Return smallest solution a, noticing that (p-a) is also a solution.
    Using Tonelli-Shanks algorithm.
    Details see: http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
    """

    n %= p

    if p == 2:
        return n
    elif is_square(n):
        r = int(sqrt(n))
    elif legendre_symbol(n, p) != 1:
        raise ValueError("n is not a quadratic residue of p!")
    elif p % 4 == 3:
        r = pow(n, (p + 1) >> 2, p)
    else:
        z = 1
        while legendre_symbol(z, p) != -1:
            z += 1

        q, s = p - 1, 0
        while q & 1 == 0:
            q >>= 1
            s += 1

        c = pow(z, q, p)
        r = pow(n, (q + 1) >> 1, p)
        t = pow(n, q, p)
        while t > 1:
            t2 = t
            for i in range(1, s):
                t2 = (t2 * t2) % p
                if t2 == 1:
                    break

            b = pow(c, 1 << (s-i-1), p)
            r = (r * b) % p
            c = (b * b) % p
            t = (t * c) % p
            s = i

    if 2 * r + 1 > p:
        r = p - r
    return r

def square_modulo_prime_power_equation(n, p, k):
    """
    Solve linear modular equation
        (x^2) % p^k = n
    or written in modular arithmatic
         x^2 = n (mod p^k)
    where n is a quadratic residue of p
    Return all solutions betweem 0 and p^k-1, since there may exist more than 2 solutions.
    Using Tonelli-Shanks algorithm and Hensel's Lift.
    Details see:
    - http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
    - https://en.wikipedia.org/wiki/Hensel%27s_lemma
    """

    r = square_modulo_prime_equation(n, p)
    sols = {r, (-r) % p}
    p_power = p
    n_power = 1
    while n_power < k:
        sols_new = set()
        p_power_new = p_power * p
        n_power += 1
        for r in sols:
            f = (r*r - n) % p_power_new
            df = (2 * r) % p
            if df:
                r = (r - f * pow(df, p_power_new-p_power-1, p_power_new)) % p_power_new
                sols_new.add(r)
            elif f == 0:
                while r < p_power_new:
                    sols_new.add(r)
                    r += p_power
        p_power = p_power_new
        sols = sols_new
    return sorted(sols)

def chinese_remainder_theorem(equation_sets):
    """
    return the only solution 0 <= x < (m_1 * m_2 * ... * m_r) for the set of simultaneous congruences:
        x = a_i (mod m_i),  1 <= i <= r
    where m_i are pairwise relatively prime
    equation_sets must in the form [(a_1, m_1), (a2, m_2), ...]
    details see: https://en.wikipedia.org/wiki/Chinese_remainder_theorem
    """

    _, ms = list(zip(*equation_sets))
    M = cprod(ms)

    x = 0
    for ai, mi in equation_sets:
        u = M // mi
        x += (ai * u * int(invert(u, mi))) % M
        x %= M
    return x

def _pqa(d, p, q):
    """
    PQa algorithm for solving generalized Pell equation, which can be regarded as a
    generalized continued fraction expansion with convergent computation
    """

    sd = sqrt(d)
    a, i = int((p+sd)/q), 0
    X, Y, PQ = [q, a*q - p], [0, 1], [(p, q)]
    while 1:
        p = a*q - p
        q = (d - p*p) // q
        a = int((p + sd) / q)

        if i == 0 and (p, q) in PQ:
            l = len(PQ) - PQ.index((p, q))
            i = 2 * len(PQ) - PQ.index((p, q))
        if len(PQ) == i:
            return l, X[1:-1], Y[1:-1], list(zip(*PQ))[1][1:]

        X.append(a * X[-1] + X[-2])
        Y.append(a * Y[-1] + Y[-2])
        PQ.append((p, q))

def generalized_pell_equation_base(d, n=1):
    """
    Solve generalized Pell equation
        x**2 - d * y**2 = n
    Return smallest positive basic solution set, or enough solutions according to nsol
    """

    if d <= 0:
        raise ValueError("D must be positive non-perfect-square integer!")
    sd = int(sqrt(d))
    if sd * sd == d:
        raise ValueError("D must be positive non-perfect-square integer!")

    # Classical Pell Equation: x**2 - d * y**2 = 1 (or -1)
    # Using continued fraction expansion of d to solve

    if abs(n) == 1:
        l, X, Y, Q = _pqa(d, 0, 1)
        if l & 1:
            if n == 1:
                x, y = X[2*l-1], Y[2*l-1]
            else:
                x, y = X[l-1], Y[l-1]
        else:
            if n == 1:
                x, y = X[l-1], Y[l-1]
            else:
                x, y = 0, 0
        return [(x, y)]

    # Generalized Pell Equation: x**2 - d * y**2 = n (n != 0)
    # Using Lagrange-Matthews-Mollin (LMM) algorithm

    # 1. Find all int f > 0 satisfying:
    #     n % (f*f) == 0
    # Note: need high efficient prime divisor decompsition algorithm when n is large.

    # 2. For each f, set m = abs(n / (f*f)).

    # 3. For each m, find all int z satisfying:
    #     (z*z) % m = d % m
    #     -m/2 < z <= m/2
    # Note: need high efficient quadratic residue algorithm to solve the first modular equation when m is large.

    # Here we just use brute-search to find all value of z.

    f, zdict = 1, {}
    while 2*f*f <= abs(n):
        if n % (f*f) == 0:
            m = n // (f*f)
            ma = abs(m)
            mb = int(ma/2)
            for z in range(-mb-(ma&1)+1, mb+1):
                if z*z % abs(m) == d % abs(m):
                    if (f, m) in zdict:
                        zdict[(f, m)].append(z)
                    else:
                        zdict[(f, m)] = [z]
        f += 1

    f = int(sqrt(abs(n)))
    if f*f == abs(n):
        zdict[(f, n//abs(n))] = [0]

    # 4. For each z according to each (f, m), run pqa(d, z, abs(m)).

    # 5. Search for first Qi = 1 or -1.

    # 6. If Xi**2 - d * Y**2 = m, (f * Xi, f * Yi) is a fundamental solution.

    # 7. Let (t, u) be the minimal positive solution of x**2 - d * y**2 = -1.
    #    If Xi**2 - d * Y**2 = -m, then (f*(t*X[i] + d*u*Y[i]), f*(u*X[i] + t*Y[i])) is a fundamental solution.

    # 8. Let (r, s) be the minimal positive solution of x**2 - d * y**2 = 1.
    #    If some fundmental solution (x, y) are not positive, we can use transformation:
    #    x' + y'*sqrt(d) = (x + y*sqrt(d)) * (1 or -1) * (r + s*sqrt(d))**k
    #    to find the minimal positive solution (x', y') in the equivalent class of (x, y)

    # 9. When all z are done, then we have a set of fundmental solutions (or minimal positive solutions) which are all in different equivalent classes.

    r, s = generalized_pell_equation_base(d, 1)[0]
    t, u = generalized_pell_equation_base(d, -1)[0]

    sols = []
    for (f, m), zlist in zdict.items():
        for z in zlist:
            l, X, Y, Q = _pqa(d, z, abs(m))
            for i, q in enumerate(Q):
                if q in (-1, 1):
                    x, y = X[i], Y[i]
                    diff = x**2 - d * y**2
                    if diff == m:
                        xg, yg = f*x, f*y
                        while xg < 0 or yg < 0:
                            if xg < 0 and yg < 0:
                                xg, yg = -xg, -yg
                            else:
                                xg, yg = r*xg + d*s*yg, s*xg + r*yg
                        sols.append((xg, yg))
                    elif t and u and diff == -m:
                        xg, yg = f*(t*x + d*u*y), f*(u*x + t*y)
                        while xg < 0 or yg < 0:
                            if xg < 0 and yg < 0:
                                xg, yg = -xg, -yg
                            else:
                                xg, yg = r*xg + d*s*yg, s*xg + r*yg
                        sols.append((xg, yg))
                    break

    # 10. Let (r, s) be the minimal positive solution of x**2 - d * y**2 = 1.
    #     We can expand the solution set.

    sols = sorted(sols)
    return sols

def generalized_pell_equation_generator(d, n=1):
    """
    Solve generalized Pell equation
        x**2 - d * y**2 = n
    generate many solutions
    """
    r, s = generalized_pell_equation_base(d, 1)[0]
    sols = generalized_pell_equation_base(d, n)

    if not sols or sols == [(0, 0)]:
        raise ValueError("No solution for x^2 - {} y^2 = {}.".format(d, n))

    while True:
        for sol in sols:
            yield sol

        sols = [(r*x + d*s*y, s*x + r*y) for x, y in sols]

def berlekamp_massey(sequence):
    """
    Find minimum linear recurrence equation of sequence using Berlekamp-Massey algorithm.
    Return C = [c0, c1, ..., cL(=1)] with c0*si + c1*s(i+1) + ... + cL*s(i+L) = 0 for i+L >= n/2.
    
    @sequence (np.array)
    """

    dtype = sequence.dtype

    # C = [c0, c1, c2, ...] letting c0*sn + c1*s(n-1) + ... = 0 for any n
    C = np.ones(1, dtype=dtype)
    B = None
    b = 0

    for n in range(len(sequence)):
        d = 0
        for i in range(len(C)):
            d += C[i] * sequence[n-i]

        B = np.append(0, B)  # align to n
        if d == 0:
            continue
        elif b == 0:
            B = C.copy()  # initialize
            b = d
        else:
            B2 = C.copy()
            C = poly_add(C, (-1 * d / b * B))
            if len(B) > len(B2):  # obtain shortest B
                B = B2
                b = d

    # rearrange C as from lower terms to higher terms 
    return C[::-1]

def berlekamp_massey_with_bound(sequence, n):
    """
    Find minimum linear recurrence equation of sequence using Berlekamp-Massey algorithm.
    Return A = [a0, a1, ..., aL(=1)] with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0 for i+L >= n.
    
    @sequence (np.array): sequence with at least 2*n terms
    @n (int): dimension upper bound of the recurrence equation
    """
    if len(sequence) < 2 * n:
        raise ValueError("sequence must with at least 2*n terms")
    dtype = sequence.dtype
    R0 = np.zeros(2*n+1, dtype=dtype)
    R0[-1] = 1
    R1 = poly_truncate(sequence[:2*n], direction="right")
    V0 = np.zeros(1, dtype=dtype)
    V1 = np.ones(1, dtype=dtype)

    while n <= len(R1) - 1:
        Q, R = poly_divmod(R0, R1)
        V = poly_add(V0, (-1) * poly_mul(Q, V1))
        V = poly_truncate(V, direction="right")
        V0, V1, R0, R1 = V1, V, R1, R

    # Let V1 = [v0, v1, ..., vL], then for i > max(deg(V1), deg(R1)+1)
    # we have v0*si + v1*s(i-1) + ... + vL*s(i-L) = 0 holds.
    # Now just rearrange and normalize it to A = [a0, a1, ..., aL(=1)]
    # with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0.

    A = V1[::-1]
    A = A / A[-1]
    return A

def berlekamp_massey_mod_p(sequence, p):
    """
    Find minimum linear recurrence equation of sequence in Z/Zp using Berlekamp-Massey algorithm.
    Return C = [c0, c1, ..., cL(=1)] with c0*si + c1*s(i+1) + ... + cL*s(i+L) = 0 for i+L >= len(seq)/2.
    @sequence (np.array)
    @p (int): prime number of field Z/Zp
    """

    dtype = sequence.dtype

    # C = [c0, c1, c2, ...] letting c0*sn + c1*s(n-1) + ... = 0 for any n
    C = np.ones(1, dtype=dtype)
    B = None
    b = 0

    for n in range(len(sequence)):
        d = 0
        for i in range(len(C)):
            d += C[i] * sequence[n-i]
            d %= p

        B = np.append(0, B)  # align to n
        if d == 0:
            continue
        elif b == 0:
            B = C.copy()  # initialize
            b = d
        else:
            B2 = C.copy()
            coeff = d * pow(int(b), p-2, p) % p
            C = poly_add(C, (-1 * coeff * B) % p) % p
            if len(B) > len(B2):  # obtain shortest B
                B = B2
                b = d

            # d2 = 0
            # for i in range(L):
            #     d2 += C[i] * sequence[n-i]
            #     d2 %= p
            # assert d2 == 0
                
    # rearrange C as from lower terms to higher terms 
    return C[::-1]


def berlekamp_massey_with_bound_mod_p(sequence, n, p):
    """
    Find minimum linear recurrence equation of sequence in Z/Zp using Berlekamp-Massey algorithm.
    Return A = [a0, a1, ..., aL(=1)] with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0 for i+L >= n.
    
    @sequence (np.array): sequence with at least 2*n terms
    @n (int): dimension upper bound of the recurrence equation
    @p (int): prime number of field Z/Zp
    """
    if len(sequence) < 2 * n:
        raise ValueError("sequence must with at least 2*n terms")
    dtype = sequence.dtype
    R0 = np.zeros(2*n+1, dtype=dtype)
    R0[-1] = 1
    R1 = poly_truncate(sequence[:2*n], direction="right")
    V0 = np.zeros(1, dtype=dtype)
    V1 = np.ones(1, dtype=dtype)

    while n <= len(R1) - 1:
        Q, R = poly_divmod_mod_p(R0, R1, p)
        V = poly_add(V0, (-1) * poly_mul_mod_p(Q, V1, p)) % p
        V = poly_truncate(V, direction="right")
        V0, V1, R0, R1 = V1, V, R1, R

    # Let V1 = [v0, v1, ..., vL], then for i > max(deg(V1), deg(R1)+1)
    # we have v0*si + v1*s(i-1) + ... + vL*s(i-L) = 0 holds.
    # Now just rearrange and normalize it to A = [a0, a1, ..., aL(=1)]
    # with a0*si + a1*s(i+1) + ... + aL*s(i+L) = 0.

    A = V1[::-1]
    inv = pow(int(A[-1]), p-2, p)
    A = (A * inv) % p
    return A