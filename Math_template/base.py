'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:23:16
'''
# coding: utf-8

"""
function:
   _is_square(n): return whther n is a perfect square
   _factorial(n): return n!
   _isqrt(n) : return integer square root of n
   _iroot(n, m): return integer m-th root of n, and whether n is a perfect power
   cprod(seq):  return seq[0] * seq[1] * ... * seq[-1]
   ggcd(seq): return the greatest common divisor (gcd) for n integers, where n can larger than 2
   extended_gcd(a, b): return gcd(a, b), x, y that a*x + b*y = gcd(a, b) using Extended Euclid Algorithm, where a, b > 0
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



