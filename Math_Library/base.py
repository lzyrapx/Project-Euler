'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:23:16
'''
# coding: utf-8

"""
normal functions:
   _is_square(n): return whther n is a perfect square
   _factorial(n): return n!
   _isqrt(n) : return integer square root of n
   _iroot(n, m): return integer m-th root of n, and whether n is a perfect power
   cprod(seq):  return seq[0] * seq[1] * ... * seq[-1]
   ggcd(seq): return the greatest common divisor (gcd) for n integers, where n can larger than 2
   extended_gcd(a, b): return gcd(a, b), x, y that a*x + b*y = gcd(a, b) using Extended Euclid Algorithm, where a, b > 0
   lcm(a,b): return the least common multiple of a and b
   llcm(seq): return the greatest common divisor (gcd) for n integers, where n can larger than 2
   padic_base_p(n, p): change integer n from base 10 to base p
   max_subarray(array): return max sum of any continous subarray of an array

modulo functions:
   _sum_mod(n): return n % 2 + n % 3 + ... + n % (n-1)
   inv_mod(n, m): return n^(-1) mod m using Extended Euclid Algorith
   fac_mod(n, m): return return n! % m
   iter_associate(f, x, n): f is a bivariable function following associate law, namely f(a, f(b, c)) = f(f(a, b), c)
                   return doing f iteratively for n times with identical input x, namely f(f(...f(x, x)..., x), x)
   sum_power_series_mod(i, n, m): sum of (x^i mod m) from x=1 to n, for i = 0, 1, 2, 3
   sum_floor(n,xmin, xmax): sum up (n//x) from x = xmin to xmax
   legendre_symbol(a, p): return value only contains {-1,0,1}, details see: http://en.wikipedia.org/wiki/Legendre_symbol

continuous fraction functions:
   rational_continous_frac(p, q): return continued fraction expansion of rational number p / q
   continous_frac_convergent(cfrac): given continued fraction, return series of convergents
   irrational_continous_frac(d, p=0, q=1, limit=100): return continued fraction expansion of quadratic irrational number (p + sqrt(d)) / q
                                                      details see: https://en.wikipedia.org/wiki/Quadratic_irrational
"""
import numpy as np
from math import gcd, sqrt
from collections import deque

try:
    from gmpy2 import is_square, fac, isqrt, iroot
    from gmpy2 import powmod as pow_mod
except:
    from fractions import gcd
    pow_mod = pow
    is_square = None
    fac = None
    isqrt = None
    iroot = None

try:
    from . ext.c_formula_int64 import c_sum_mod_int64
    sum_mod = c_sum_mod_int64
except:
    sum_mod = None

# Supplementry implementations
def _is_square(n):
    """return whether n is a perfect square"""
    """perfect square: 1, 4, 9, 25, 36, 49, 64..."""

    s = int(sqrt(n))
    return s * s == n

if is_square is None:
    is_square = _is_square

def _factorial(n):
    """return n!"""

    if n < 0:
        raise ValueError("n in n! must be positive!")
    if n == 0 or n == 1:
        return 1

    output = 1
    for i in range(2, n+1):
        output *= i
    return output

if fac is None:
    fac = _factorial

def _isqrt(n):
    """return integer square root of n"""

    return int(n**0.5)

if isqrt is None:
    isqrt = _isqrt


def _iroot(n, m):
    """return integer m-th root of n, and whether n is a perfect power"""

    r = int(n**(1./m))
    return r, r**m == n

if iroot is None:
    iroot = _iroot

def cprod(seq):
    """return seq[0] * seq[1] * ... * seq[-1]"""

    output = 1
    for i in iter(seq):
        output *= i
    return output

def ggcd(seq):
    """
    return the greatest common divisor (gcd) for n integers,
    where n can larger than 2
    """

    if len(seq) < 2:
        raise ValueError("There should be at least 2 integers!")
    elif len(seq) == 2:
        return gcd(seq[0], seq[1])
    else:
        g = gcd(seq[-2], seq[-1])
        if g == 1:
            return g
        for n in seq[:-2]:
            g = gcd(g, n)
            if g == 1:
                return 1
        return g

def extended_gcd(a, b):
    """
    return gcd(a, b), x, y that a*x + b*y = gcd(a, b)
    using Extended Euclid Algorithm, where a, b > 0
    """

    assert a >= 0 and b >= 0

    x = v = 0
    y = u = 1
    while a:
        q = b // a
        r = b - q * a
        m = x - u * q
        n = y - v * q

        b = a
        a = r
        x = u
        y = v
        u = m
        v = n

    gcd = b
    return gcd, x, y

def lcm(a, b):
    """return the least common multiple of a and b"""

    return a // gcd(a, b) * b

def llcm(seq):
    """
    return the greatest common divisor (gcd) for n integers,
    where n can larger than 2
    """

    if len(seq) < 2:
        raise ValueError("There should be at least 2 integers!")
    elif len(seq) == 2:
        return lcm(seq[0], seq[1])
    else:
        l = lcm(seq[-2], seq[-1])
        for n in seq[:-2]:
            if l % n:
                l = lcm(l, n)
        return l

def padic_base_p(n, p, ntype='s'):
    """change integer n from base 10 to base p"""

    base = '0123456789' + ''.join(map(chr, range(65, 92)))
    snp = ''
    while True:
        if n == 0:
            break
        n, r = divmod(n, p)
        snp += base[r]

    snp = snp[::-1]
    if ntype == 's':
        return snp
    elif ntype == 'n':
        return int(snp)
    elif ntype == 'l':
        return list(snp)

def max_subarray(array):
    """return max sum of any continous subarray of an array"""

    max_so_far = max_ending_here = 0
    for x in array:
        max_ending_here = max(0, max_ending_here + x)
        max_so_far = max(max_so_far, max_ending_here)
    return max_so_far

def _sum_mod(n):
    """return n % 2 + n % 3 + ... + n % (n-1)"""

    from itertools import takewhile, count

    sm = i = 0
    for i in takewhile(lambda x: n//x - n//(x+1) > 4, count(1)):
        a = n % (n//(i+1) + 1)
        b = n % (n//i) if i > 1 else 1
        c = (a-b) // i + 1
        sm += b*c + i*(c - 1)*c // 2
    sm += sum(n % j for j in range(2, n//(i+1) + 1))
    return sm

if sum_mod is None:
    sum_mod = _sum_mod

def inv_mod(n, m):
    """return n^(-1) mod m using Extended Euclid Algorithm"""

    n %= m
    if n == 0 or m <= 0:
        return 0
    
    m0, x0, x1 = m, 0, 1
    while n != 0:
        x0, x1 = x1, x0 - m // n * x1
        m, n = n, m % n

    if m == 1:
        return x0 % m0
    else:
        return 0

def fac_mod(n, m):
    """return n! % m"""

    if n < 0:
        raise ValueError("n in n! must be positive!")
    if n == 0 or n == 1:
        return 1

    output = 1
    for i in range(2, n+1):
        output *= i
        output %= m
    return output

def iter_associate(f, x, n):
    """
    f is a bivariable function following associate law, namely f(a, f(b, c)) = f(f(a, b), c)
    return doing f iteratively for n times with identical input x, namely f(f(...f(x, x)..., x), x)
    """

    n, r = n - 1, x
    while n:
        if n % 2:
            r = f(x, r)
        n >>= 1
        x = f(x, x)
    return r

def sum_power_series_mod(i, n, m):
    """sum of (x^i mod m) from x=1 to n, for i = 0, 1, 2, 3"""

    if i == 0:
        return n
    elif i == 1:
        n %= 2 * m
        return ((n*(n+1)) >> 1) % m
    elif i == 2:
        n %= 6 * m
        m3 = 3 * m
        res = ((n*(n+1)) >> 1) % m3
        return ((res * ((2*n+1) % m3)) % m3) // 3
    elif i == 3:
        n %= 2 * m
        res = ((n*(n+1)) >> 1) % m
        return (res * res) % m
    else:
        return 0

def sum_floor(n, xmin, xmax):
    """sum up n//x from x = xmin to xmax"""

    nrt = isqrt(n)
    res = 0
    if xmin > n:
        return 0
    if xmax > n:
        xmax = n

    if xmax <= nrt:
        for x in range(xmin, xmax+1):
            res += n // x
    elif xmin >= nrt:
        real_xmin = n // xmax
        real_xmax = n // xmin
        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, real_xmax+1):
            a0 = a1
            a1 = n // (x+1)
            ub = a0 if a0 < xmax else xmax
            lb = a1 if a1 >= xmin-1 else xmin-1
            res += (ub - lb) * x
    else:
        real_xmin = n // xmax
        if real_xmin > xmin:
            real_xmin = xmin

        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, nrt+1):
            a0 = a1
            a1 = n // (x+1)

            if x >= xmin:
                res += a0

            if a1 < xmax:
                ub = a0 if a0 < xmax else xmax
                res += (ub - a1) * x

        if x == n // x:
            res -= x
    return res

def legendre_symbol(a, p):
    """
    Legendre symbol
    Define if a is a quadratic residue of prime p
    Details see: http://en.wikipedia.org/wiki/Legendre_symbol
    """

    if p == 2 or a % p == 0:
        return 1

    ls = pow(a % p, (p - 1) >> 1, p)
    if ls == p - 1:
        return -1
    return ls

# Continuous Fraction Functions
def rational_continous_frac(p, q=1, limit=100):
    """
    return continued fraction expansion of rational number p / q
    for example: 415/93 -> (4,2,6,7)
    """

    p, q, cfrac = int(p), int(q), []
    while p % q:
        cfrac.append(p//q)
        p, q = q, p % q
    if q:
        cfrac.append(p)
    return tuple(cfrac)

def continous_frac_convergent(cfrac):
    """
    given continued fraction, return series of convergents
    for example: [4,2,6,7] -> [Fraction(4, 1), Fraction(9, 2), Fraction(58, 13), Fraction(415, 93)]
    """

    from fractions import Fraction

    if len(cfrac) < 2:
        raise ValueError("Continued fraction must longer than 2!")

    a0, a1 = cfrac[:2]
    p0, p1, q0, q1 = a0, a0*a1+1, 1, a1
    cvg = [Fraction(p0, q0), Fraction(p1, q1)]
    for a in cfrac[2:]:
        p0, p1 = p1, a*p1 + p0
        q0, q1 = q1, a*q1 + q0
        cvg.append(Fraction(p1, q1))
    return cvg

def irrational_continous_frac(d, p=0, q=1, limit=100):
    """
    return continued fraction expansion of quadratic irrational number (p + sqrt(d)) / q
    about quadratic irrational number: https://en.wikipedia.org/wiki/Quadratic_irrational
    """

    from fractions import Fraction

    d, p, q = int(d), int(p), int(q)

    sd = sqrt(d)
    if int(sd) * int(sd) == d:
        raise ValueError("D is perfect square!")

    a = int((p + sd) / q)
    repetend, pairspq = [a], []
    for _ in range(limit):
        p = Fraction(a*q - p, 1)
        q = Fraction(d - p*p, q)
        a = int((p + sd) / q)

        if (p, q) in pairspq:
            i = pairspq.index((p, q))
            return tuple(repetend[:i+1] + [tuple(repetend[i+1:])])

        pairspq.append((p, q))
        repetend.append(a)
    raise ValueError("Repetend is longer than {0:d}, please try higher limit!".format(limit))
