'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-02 22:16:04
'''
# coding: utf-8
"""
Tree Generation functions:
    pythagorean_triple_tree(triple=(3, 4, 5), forward=True, trust=False):
                    Given a PPT, it can generate three more PPTs. 
                    And if recursively applying this branching function from (3, 4, 5), 
                    all PPTs in form (odd, even, odd) can be generated in a trinary tree, 
                    which covers the entire set of PPTs completely and uniquely.
                    When forward is True, return all three children of current PPT.
                    When forward is False, return its parent in the PPT tree.
                    details see: https://en.wikipedia.org/wiki/Tree_of_primitive_Pythagorean_triples
    
"""
import numpy as np
from math import gcd, sqrt
from collections import deque


def pythagorean_triple_tree(triple=(3, 4, 5), forward=True, trust=False):
    """
    Primitive Pythagorean Triple (PPT) (a, b, c) is integer triple satisfying
        a**2 + b**2 = c**2
        gcd(a, b) = gcd(b, c) = gcd(a, c) = 1
    Given a PPT, it can generate three more PPTs. 
    And if recursively applying this branching function from (3, 4, 5), 
    all PPTs in form (odd, even, odd) can be generated in a trinary tree, 
    which covers the entire set of PPTs completely and uniquely.
    When forward is True, return all three children of current PPT.
    When forward is False, return its parent in the PPT tree.
    """

    a, b, c = triple
    if not trust:
        if a**2 + b**2 != c**2:
            raise ValueError("Invalid Primitive Pythagorean Triple")
        if gcd(a, b) * gcd(a, c) * gcd(b, c) != 1:
            raise ValueError("Invalid Primitive Pythagorean Triple")

    if forward:
        return ((a-2*b+2*c,   2*a-b+2*c,  2*a-2*b+3*c),
                (a+2*b+2*c,   2*a+b+2*c,  2*a+2*b+3*c),
                (-a+2*b+2*c, -2*a+b+2*c, -2*a+2*b+3*c))
    else:
        if triple == (3, 4, 5):
            return triple
        else:
            return (abs(-a-2*b+2*c), abs(-2*a-b+2*c), -2*a-2*b+3*c)

def co_prime_tree(pair=(0, 0), trust=False):
    """
    All co-prime pairs can be generated from (2, 1) (for (odd, even) and (even, odd) pairs) and (3, 1) (for (odd, odd) pairs).
    It follows a trinary tree, from co-prime pair (a, b), we get: (2*a - b, a), (2*a + b, a), (a + 2*b, b)
    It can be shown that the co-prime pairs in the tree are disjoint complete.
    """

    if pair == (0, 0):
        return ((2, 1), (3, 1))

    a, b = pair
    if not trust:
        if gcd(a, b) != 1:
            raise ValueError("Invalid co-prime pair!")
    return ((2*a - b, a), (2*a + b, a), (a + 2*b, b))

def stern_brocot_tree():
    """
    Stern-Brocot Tree, an infinite complete binary tree in which the vertices 
    correspond one-for-one to the positive rational numbers, 
    whose values are ordered from left to right as in a search tree. It related to Farey series closely.
    details see: https://en.wikipedia.org/wiki/Stern%E2%80%93Brocot_tree
    """

    sbt = deque([1, 1])
    while True:
        sbt += [sbt[0] + sbt[1], sbt[1]]
        sbt += [sbt[1] + sbt[2], sbt[2]]
        yield (sbt.popleft(), sbt.popleft())