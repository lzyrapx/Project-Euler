'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-04 00:16:33
'''
#coding: utf-8

"""
prime-related functions:
    _primes_list(n): Returns a array of primes, 2 <= p < n
    _is_prime(p, accuracy=100): Miller-Rabin primality test
                                          https://en.wikipedia.org/wiki/Miller-Rabin_primality_test
    _mobius_list(n): return mobius function mu(k) for 0 <= k <= 
"""
import random
import numpy as np

from base import cprod, gcd, sqrt, isqrt

try:
    from gmpy2 import is_prime
except:
    is_prime = None

primes_list = None
mobius_list = None
factor_sieve = None

def _primes_list(n): # n >= 6
    """input n>=6, Returns a array of primes, 2 <= p < n"""

    sieve = np.ones(n//3 + (n % 6 == 2), dtype=np.bool)
    for i in range(1, int(sqrt(n))//3+1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[k*k//3::2*k] = False
            sieve[k*(k - 2*(i & 1) + 4)//3::2*k] = False
    plist = np.r_[2, 3, ((3 * np.nonzero(sieve)[0][1:] + 1) | 1)]
    return [int(x) for x in plist]

if primes_list is None:
    primes_list = _primes_list

def _mr_decompose(n):
    exponentOfTwo = 0
    while n % 2 == 0:
        n //= 2
        exponentOfTwo += 1
    return exponentOfTwo, n

def _mr_isWitness(possibleWitness, p, exponent, remainder):
    possibleWitness = pow(possibleWitness, remainder, p)
    if possibleWitness == 1 or possibleWitness == p - 1:
        return False
    for _ in range(exponent):
        possibleWitness = pow(possibleWitness, 2, p)
        if possibleWitness == p - 1:
            return False
    return True

def _is_prime(p, accuracy=100):
    """
    Miller-Rabin primality test
    https://en.wikipedia.org/wiki/Miller-Rabin_primality_test
    """

    if p < 2:
        return False
    if p == 2 or p == 3:
        return True
    
    exponent, remainder = _mr_decompose(p - 1)

    for _ in range(accuracy):
        possibleWitness = random.randint(2, p - 2)
        if _mr_isWitness(possibleWitness, p, exponent, remainder):
            return False
    return True

if is_prime is None:
    is_prime = _is_prime

def _mobius_list(n):
    """return mobius function mu(k) for 0 <= k <= n"""

    plist = primes_list(isqrt(n)+1)
    mlist = np.ones(n+1, dtype=np.int64)

    for p in plist:
        mlist[::p] *= -p
        mlist[::p*p] = 0

    for i in range(1, n+1):
        if mlist[i]:
            if abs(mlist[i]) < i:
                mlist[i] *= -1

            if mlist[i] > 0:
                mlist[i] = 1
            else:
                mlist[i] = -1
    return mlist

if mobius_list is None:
    mobius_list = _mobius_list