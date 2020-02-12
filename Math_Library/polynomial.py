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
    poly_add(poly1, poly2): return poly1 + poly2
                            @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_mul(poly1, poly2): return poly1 * poly2
                            @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_divmod(poly1, poly2): return poly1 / poly2, poly1 % poly2
                               @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_mul_mod_p(poly1, poly2, p): return (poly1 * poly2) % p
                                     @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_divmod_mod_p(poly1, poly2, p): return poly1 / poly2, poly1 % poly2 in Z/Zp
                                        @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    poly_Zp_pow_mod(poly, n, polymod, p): return poly^n % polymod in Z/Zp
                                          @poly, polymod (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
                                          @n (int): power
                                          @p (int): prime number of field Z/Zp
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
    return poly1 + poly2
    
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
    return poly1 * poly2
    
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

def poly_divmod(poly1, poly2):
    """
    return poly1 / poly2, poly1 % poly2
    
    @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    """

    l1 = len(poly1)
    l2 = len(poly2)
    
    if l1 < l2:
        return np.zeros(1, dtype=poly1.dtype), poly1

    quotient = np.zeros(l1 - l2 + 1, dtype=poly1.dtype)
    remaining = poly1.copy()
    for i in range(l1-1, l2-2, -1):
        quotient[i-l2+1] = remaining[i] / poly2[-1]
        remaining[i-l2+1:i+1] -= poly2 * quotient[i-l2+1]
    return quotient, poly_truncate(remaining, direction="right")

def poly_mul_mod_p(poly1, poly2, p):
    """
    return (poly1 * poly2) % p
    
    @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    """

    l1 = len(poly1)
    l2 = len(poly2)
    res = np.zeros(l1 + l2 - 1, dtype=poly1.dtype)
    if l1 > l2:
        for i, n in enumerate(poly2):
            res[i:l1+i] = (res[i:l1+i] + (poly1 * n) % p) % p
    else:
        for i, n in enumerate(poly1):
            res[i:l2+i] = (res[i:l2+i] + (poly2 * n) % p) % p
    return res


def poly_divmod_mod_p(poly1, poly2, p):
    """
    return poly1 / poly2, poly1 % poly2 in Z/Zp
    
    @poly1, poly2 (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    """

    l1 = len(poly1)
    l2 = len(poly2)

    if l1 < l2:
        return np.zeros(1, dtype=poly1.dtype), poly1

    inv = pow(int(poly2[-1]), p-2, p)

    quotient = np.zeros(l1 - l2 + 1, dtype=poly1.dtype)
    remaining = poly1.copy()
    for i in range(l1-1, l2-2, -1):
        quotient[i-l2+1] = (remaining[i] * inv) % p
        remaining[i-l2+1:i+1] = (remaining[i-l2+1:i+1] - (poly2 * quotient[i-l2+1]) % p) % p
    return quotient, poly_truncate(remaining, direction="right")


def poly_Zp_pow_mod(poly, n, polymod, p):
    """
    return poly^n % polymod in Z/Zp
    
    @poly, polymod (np.array): polynomial coefficients, representing sum_i poly[i] * x^i
    @n (int): power
    @p (int): prime number of field Z/Zp
    """

    res = np.ones(1, dtype=poly.dtype)
    while n:
        if n & 1:
            res = poly_mul_mod_p(res, poly, p)
            _, res = poly_divmod_mod_p(res, polymod, p)

        n >>= 1
        if n:
            poly = poly_mul_mod_p(poly, poly, p)
            _, poly = poly_divmod_mod_p(poly, polymod, p)

    return poly_truncate(res, direction="right")