'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-06 22:44:54
'''

"""
functions implementing maths with polynomials.
function list:
    poly_truncate(poly, direction="left"): truncate zero terms on left or right
                                           @poly (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
                                           @direction (str): left or right, left is for lower powers, right is for higher term
    poly_add(poly1, poly2): poly1 + poly2
                            @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_mul(poly1, poly2): poly1 * poly2
                            @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    

"""
import math
import numpy as np

def poly_truncate(poly, direction="left"):
    """
    truncate zero terms on left or right
    
    @poly (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    @direction (str): left or right, left is for lower powers, right is for higher terms
    """

    if direction == "left":
        i = 0
        while i + 1 < len(poly) and poly[i] == 0:
            i += 1
        return poly[i:]
    elif direction == "right":
        i = len(poly)
        while i > 1 and poly[i-1] == 0:
            i -= 1
        return poly[:i]

def poly_add(poly1, poly2):
    """
    poly1 + poly2
    
    @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    """

    l1 = len(poly1)
    l2 = len(poly2)
    if l1 > l2:
        res = poly1.copy()
        res[:l2] += poly2
    else:
        res = poly2.copy()
        res[:l1] += poly1
    return res

def poly_mul(poly1, poly2):
    """
    poly1 * poly2
    
    @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    """

    l1 = len(poly1)
    l2 = len(poly2)
    res = np.zeros(l1 + l2 - 1, dtype=poly1.dtype)
    if l1 > l2:
        for i, n in enumerate(poly2):
            res[i:l1+i] += poly1 * n
    else:
        for i, n in enumerate(poly1):
            res[i:l2+i] += poly2 * n
    return res
