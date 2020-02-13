
"""
Functions implementing maths in linear algebra.

function list:
    dot_mod(A, B, mod=0): return matrix multiplication of A * B % mod
                         matrix multiplication, avoid overflow in numpy
                         @A, B(np.array)
    dot_mod_as_list(A, B, mod=0): return matrix multiplication of A * B % mod
                                  matrix multiplication defined as list , avoid overflow in numpy
                                  @A, B(list)
    mat_pow_mod(mat, n, mod=0): return (mat^n) % mod
                                @mat(np.array)
"""
#coding: utf-8

from copy import deepcopy

import numpy as np
from sympy import Symbol, Rational
from math import gcd
from base import  inv_mod


def dot_mod(A, B, mod=0):
    """
    return matrix multiplication of A * B % mod
    matrix multiplication, avoid overflow in numpy
    @A, B(np.array)
    """

    a = len(A)
    l = len(B)
    b = len(B[0])

    C = np.zeros((a, b), dtype=np.int64)
    for i in range(a):
        for j in range(b):
            cij = 0
            for k in range(l):
                if mod:
                    cij = (cij + A[i, k] * B[k, j]) % mod
                else:
                    cij += A[i, k] * B[k, j]
            C[i, j] = cij
    return C

def dot_mod_as_list(A, B, mod=0):
    """
    return matrix multiplication of A * B % mod
    matrix multiplication defined as list , avoid overflow in numpy
    @A, B(list)
    """
    a = len(A)
    l = len(B)
    b = len(B[0])

    C = [[0] * b for _ in range(a)]
    for i in range(a):
        for j in range(b):
            cij = 0
            for k in range(l):
                if mod:
                    cij = (cij + A[i][k] * B[k][j]) % mod
                else:
                    cij += A[i][k] * B[k][j]
            C[i][j] = cij
    return C

def mat_pow_mod(mat, n, mod=0):
    """
    return (mat^n) % mod
    """

    if n < 0:
        raise ValueError("power must be positive!")

    d = len(mat)
    res = np.eye(d, dtype=np.int64)
    while n:
        if n & 1:
            if mod:
                res = np.mod(np.dot(res, mat), mod)
            else:
                res = np.dot(res, mat)
        if mod:
            mat = np.mod(np.dot(mat, mat), mod)
        else:
            mat = np.dot(mat, mat)
        n >>= 1
    return res