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
                                details see: https://en.wikipedia.org/wiki/Miller-Rabin_primality_test
    _is_coprime(a, b): return whether a and b are coprime
    _mobius_list(n): return mobius function mu(k) for 0 <= k <= n
    _pollard_rho(n, rand=True): return a non-trivial(not one or n) factor of n
                                Pollard rho prime factorization algorithm
                                details see: https://en.wikipedia.org/wiki/Pollard's_rho_algorithm
    prime_divisor_decomposition(n, rand=True):
                                Prime factor decomposition
                                writing n as a product of prime factors. To factorise a number, 
                                divide it by the first possible prime number.
    all_divisors(n, rand=False): return all divisors of n as a sorted list
    euler_phi(n, rand=False): return Euler's totient value of n
                              http://oeis.org/A000010
                              Details see: https://en.wikipedia.org/wiki/Euler's_totient_function
    mobius(n):  return mobius function mu(n)
                details see: http://oeis.org/A008683
    _largest_prime_factor_sieve(n): return largest prime divisor of n <= sqrt(n), 1 if n is prime or 1. for 0 <= k <= n
                                    details: http://oeis.org/A217581
    prime_counting(n): a simple implementation of extended Meissel-Lehmer algorithm
                       return number of prime numbers <= x, with both time and space complexity O(x^2/3)
                       details see: http://oeis.org/A006880
    atkin_prime_sieve(limit=1000000): finding all prime numbers up to limit, O(n) time complexity, O(n) memory.
                                      Sieve of Atkin, which is a much stronger version than sieve of Eratosthenes
                                      details see: https://en.wikipedia.org/wiki/Sieve_of_Atkin
"""
import random
import numpy as np

from base import cprod, gcd, sqrt, isqrt

try:
    from gmpy2 import is_prime
    from gmpy2 import sqrt, gcd
except:
    from math import sqrt
    from fractions import gcd
    is_prime = None

is_coprime = None
primes_list = None
mobius_list = None
largest_prime_factor_sieve = None

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

def _is_coprime(a, b):
    """return whether a and b are coprime"""

    return gcd(a, b) == 1

if is_coprime is None:
    is_coprime = _is_coprime

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

def _pollard_rho(n, rand=True):
    """
    return a non-trivial(not one or n) factor of n.
    Pollard rho prime factorization algorithm
    https://en.wikipedia.org/wiki/Pollard's_rho_algorithm
    """

    f = lambda x, c: x*x + c
    if not rand:
        x, c = 1, 1
    else:
        x, c = random.randrange(2, 1e6), random.randrange(2, 1e6)

    y, d = x, 1
    while d == 1 and d != n:
        x = f(x, c) % n
        y = f(y, c) % n
        y = f(y, c) % n
        d = gcd(y-x, n)
    return int(d)

P10K = primes_list(10000)
P10Kset = set(P10K)
def prime_divisor_decomposition(n, rand=True):
    dlist, clist = [], []

    # 奇偶性判断
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    if c:
        dlist.append(2)
        clist.append(c)

    # 首先用10000以内的小素数试除
    for p in iter(P10K):
        c = 0
        while n % p == 0:
            n //= p
            c += 1
        if c:
            dlist.append(p)
            clist.append(c)

        if n == 1:
            return list(zip(dlist, clist))

        if n in P10Kset:  # set的in操作复杂度<=O(log(n))
            dlist.append(n)
            clist.append(1)
            return list(zip(dlist, clist))

        n = int(n)

    # 然后用Pollard rho方法生成素因子
    while 1:
        if n == 1:
            return list(zip(dlist, clist))

        if is_prime(n):
            dlist.append(n)
            clist.append(1)
            return list(zip(dlist, clist))

        p = _pollard_rho(n, rand)
        c = 0
        while n % p == 0:
            n //= p
            c += 1
        dlist.append(p)
        clist.append(c)

        n = int(n)

def all_divisors(n, rand=False):
    if n == 1:
        return [1]

    primefactors = prime_divisor_decomposition(n, rand)
    d = len(primefactors)
    clist = [0] * d
    output = []
    while 1:
        output.append(cprod([primefactors[i][0]**clist[i] for i in range(d)]))
        k = 0
        while 1:
            clist[k] += 1
            if clist[k] <= primefactors[k][1]:
                break
            clist[k] = 0
            k += 1
            if k >= d:
                return sorted(output)

def euler_phi(n, rand=False):
    """
    Return Euler's totient value of n
    http://oeis.org/A000010
    Details see: https://en.wikipedia.org/wiki/Euler's_totient_function
    """

    if n == 1:
        return 1

    primefactors = prime_divisor_decomposition(n, rand)
    phi = 1
    for p, a in primefactors:
        phi *= pow(p, a-1) * (p - 1)
    return phi

def mobius(n):
    """
    return mobius function mu(n)
    details see: http://oeis.org/A008683
    """

    flist = prime_divisor_decomposition(n)
    for p, a in flist:
        if a > 1:
            return 0

    if len(flist) & 1:
        return -1
    else:
        return 1

def _largest_prime_factor_sieve(n):
    """
    return largest prime divisor of n <= sqrt(n), 1 if n is prime or 1. for 0 <= k <= n
    details: http://oeis.org/A217581
    """

    fac = np.ones(n+1, dtype=np.int64)
    for p in range(2, n):
        if p * p > n:
            break
        if fac[p] == 1:
            fac[p*p::p] = p
    return fac

if largest_prime_factor_sieve is None:
    largest_prime_factor_sieve = _largest_prime_factor_sieve

def prime_counting(x):
    """
    a simple implementation of extended Meissel-Lehmer algorithm
    return number of prime numbers <= x, with both time and space complexity O(x^2/3)
    details see: http://oeis.org/A006880
    """

    y = int(x**(1./3))
    x_sqrt = int(isqrt(x))
    ub = x // y

    primes = primes_list(x_sqrt+1)
    mudelta = np.ones((2, y+1), dtype=np.int)  # [mu, smallest prime factor]
    pi_list = np.ones((3, ub+1), dtype=np.int) # [pi, fenwick tree of phi, sieved-or-not for phi]

    # tool functions for Fenwick Tree
    def add(i, val, fenwick):
        while i <= ub:
            fenwick[i] += val
            i += i & -i
        return i

    def sum_range(i, fenwick):
        s = 0
        while i:
            s += fenwick[i]
            i -= i & -i
        return int(s)

    # initialize pi_list
    pi_list[0, :2] = 0
    for n in range(2, ub+1):
        # linear sieve of pi(n) for n <= ub
        pi_list[0, n] += pi_list[0, n-1]
        for p in primes:
            if p > n or p > ub // n or p * p > ub:
                break
            else:
                pi_list[0, n * p] = 0
                if n % p == 0:
                    break

        # initialize Fenwick Tree for phi(m, b)
        if n & 1 == 0:
            m = 2
            while n % m == 0:
                m <<= 1
            pi_list[1, n] = m >> 1

    # initialize mu and delta
    for p in range(2, y+1):
        if mudelta[1, p] == 1:
            # mu
            mudelta[0, p::p] *= -1
            m = p * p
            if m <= y:
                mudelta[0, m::m] = 0

            # delta
            mudelta[1, p] = p
            m = p * p
            while m <= y:
                if mudelta[1, m] == 1:
                    mudelta[1, m] = p
                m += p + p if p > 2 else p

    res = 0

    # handle a - 1 - P2(x, a)
    a = int(pi_list[0, y])
    b = int(pi_list[0, x_sqrt])
    res += a - 1 - a*(a-1)//2 + b*(b-1)//2
    for p in primes:
        if p > x_sqrt:
            break
        elif p > y:
            res -= int(pi_list[0, x//p])

    # handle phi(x, a)
    # handle S0
    for m in range(1, y+1):
        if mudelta[0, m]:
            res += int(mudelta[0, m]) * (x // m)

    # handle S
    p = 0
    for q in primes:
        if p > y:
            break

        # sieve out p
        if p:
            pi_list[2, p] = 0
            add(p, -1, pi_list[1])

            # sieve out multiples of p
            m = p * p
            while m <= ub:
                if pi_list[2, m]:
                    add(m, -1, pi_list[1])
                    pi_list[2, m] = 0

                m += p + p if p > 2 else p  # acclererate a little
        
        for m in range(y//q+1, y+1):
            # m <= y < mq <= x
            if mudelta[0, m] and mudelta[1, m] > q:
                res += -int(mudelta[0, m]) * sum_range(x // m // q, pi_list[1])
        p = q

    return res

def atkin_prime_sieve(limit=1000000):
    """
    finding all prime numbers up to limit. O(n) time complexity, O(n) memory
    Sieve of Atkin, which is a much stronger version than sieve of Eratosthenes
    details see: https://en.wikipedia.org/wiki/Sieve_of_Atkin
    """

    plist = [0] * limit

    # n = 3x^2 + y^2 section
    x = 3
    for i in range(0, 12*int(sqrt((limit-1)/3)), 24):
        x += i
        y_limit = int(12*sqrt(limit-x)-36)
        n = x + 16
        for j in range(-12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

        n = x + 4
        for j in range( 12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

    # n = 4x^2 + y^2 section
    x = 0
    for i in range(4, 8*int(sqrt((limit-1)/4))+4, 8):
        x += i
        n = x + 1
        if x % 3:
            for j in range(0, 4*int(sqrt(limit-x))-3, 8):
                n += j
                plist[n] = not plist[n]
        else:
            y_limit = 12 * int(sqrt(limit-x)) - 36
            
            n = x + 25
            for j in range(-24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

            n = x + 1
            for j in range( 24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

    # n = 3x^2 - y^2 section
    x = 1
    for i in range(3, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x
        if n > limit:
            y = (int(sqrt(n-limit)) >> 2) << 2
            n -= y * y
            s = 4 * y + 4
        else:
            s = 4

        for j in range(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    x = 0
    for i in range(2, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x

        if n > limit:
            y = ((int(sqrt(n-limit)) >> 2) << 2) - 1
            n -= y * y
            s = 4 * y + 4
        else:
            n -= 1
            s = 0

        for j in range(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    # eliminate squares        
    for n in range(5, int(sqrt(limit))+1, 2):
        if plist[n]:
            for k in range(n*n, limit, n*n):
                plist[k] = False
    return [2,3] + list(filter(plist.__getitem__, range(5,limit,2)))