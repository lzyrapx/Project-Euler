'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:23:16
'''
# coding: utf-8

"""
normal function:
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

modulo Functions:
   _sum_mod: return n % 2 + n % 3 + ... + n % (n-1)
   inv_mod: return n^(-1) mod m using Extended Euclid Algorith
   fac_mod: return return n! % m
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







