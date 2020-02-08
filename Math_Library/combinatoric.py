'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-06 22:08:13
'''
# coding: utf-8

from base import gcd, fac_mod, fac
from prime import primes_list, euler_phi

"""
combinatoric related functions:
    C(n, k): return C(n,k) = n!/k!(n-k)!, example: C(3, 1) = 3, C(4, 2) = 12
    C_mod(n, k, m): return C(n, k) % m.
    permutations_multiset(n): Calculate number of permutations of multiset, which defined by
                            n = [n1, n2, ...] or n = [n1, n1, n2...]
    list_multiset_permutations(multiset): List out all permutations of a multiset [a, a, ..., b, b, ..., c, c, ...]
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

def permutations_multiset(n): # n is a list
    """
    Calculate number of permutations of multiset, which defined by
    n = [n1, n2, ...] or n = [n1, n1, n2...]
    ps. 多重集的排列，n1 为某一类的重数, n = n1 + n2 + ... + nk
        P(n; n1*a1, n2*a2, ..., nk*ak) = n! / (n1! * n2! * ...* nk!), 其中，ai 是元素
        另外，选取 r 个元素的多重集的组合是：C(n; n1*a1, n2*a2, ..., nk*ak) = C(k+r-1, r)，k 是元素类型种数
    """
    s, p = 0, 1
    for v in n:
        s += v
        p *= fac(v)
    return fac(s) // p

def list_multiset_permutations(multiset):
    """
    List out all permutations of a multiset [a, a, ..., b, b, ..., c, c, ...]
    """

    class Node(object):
        def __init__(self, v, to=None):
            self.val = v
            self.to = to

    def visit(n):
        vlist = [n.val]
        while n.to:
            n = n.to
            vlist.append(n.val)
        return vlist

    if len(multiset) == 1:
        yield multiset
    elif len(multiset) == 2:
        if multiset[0] == multiset[1]:
            yield multiset
        else:
            yield multiset
            yield multiset[::-1]
    else:
        E = [Node(n) for n in sorted(multiset, reverse=True)]
        for i, n in enumerate(E[:-1]):
            n.to = E[i+1]

        head, i, afteri = E[0], E[-2], E[-1]
        yield visit(head)

        while afteri.to or afteri.val < head.val:
            if afteri.to and i.val >= afteri.to.val:
                beforek = afteri
            else:
                beforek = i
            k = beforek.to
            beforek.to = k.to
            k.to = head
            if k.val < head.val:
                i = k
            head, afteri = k, i.to
            yield visit(head)
