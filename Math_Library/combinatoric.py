'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-06 22:08:13
'''
# coding: utf-8

from base import gcd, fac_mod
from prime import primes_list, euler_phi

"""
combinatoric related functions:
    C(n, k): return C(n,k) = n!/k!(n-k)!, example: C(3, 1) = 3, C(4, 2) = 12
    C_mod(n, k, m): return C(n, k) % m.
"""

def C(n, k):
    """
    return C(n,k)=n!/k!(n-k)!
    """
    if k > n:
        return 0
    if k > n/2:
        k = n - k
    if k == 0:
        return 1
    if k == 1:
        return n

    output = n
    for i in range(n-1, n-k, -1):
        output *= i
    for i in range(2, k+1):
        output //= i
    return output

def C_mod(n, k, m):
    """return C(n, k) % m"""

    if k > n:
        return 0
    if k > n//2:
        k = n - k
    if k == 0:
        return 1
    if k == 1:
        return n % m

    output = fac_mod(n, m)

    x = (fac_mod(k, m) * fac_mod(n - k, m)) % m
    if x and gcd(x, m) == 1:
        t = euler_phi(m)
        output *= pow(x, t-1, m)
        output %= m
    else:
        plist = primes_list(n+1)
        output = 1
        for p in plist:
            x = 0
            q = p
            while q <= n:
                x += n // q
                q *= p

            if p <= k:
                q = p
                while q <= k:
                    x -= k // q
                    q *= p

            if p <= n - k:
                q = p
                while q <= n - k:
                    x -= (n - k) // q
                    q *= p

            output *= pow(p, x, m)
            output %= m
    return output