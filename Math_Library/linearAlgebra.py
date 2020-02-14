'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-14 21:44:54
'''

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
    mat_pow_mod_as_list(mat, n, mod=0): return (mat^n) % mod
                           mat is defined as list, avoid overflow in numpy
                           @mat(list)
    mat_sum_pow_mod(A0, Q, n, mod=0): return (A0 + Q A0 + Q^2 A0 + ... + Q^n A0) % mod
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

def mat_pow_mod(mat, n, Mod=0):
    """
    return (mat^n) %  Mod
    """

    if n < 0:
        raise ValueError("power n must be positive!")

    d = len(mat)
    l = len(mat[0])
    res = np.eye(d, dtype=np.int64)
    while n:
        if n & 1:
            if Mod:
                res = np.mod(np.dot(res, mat), Mod)
            else:
                res = np.dot(res, mat)
        if Mod:
            mat = np.mod(np.dot(mat, mat), Mod)
        else:
            mat = np.dot(mat, mat)
        n >>= 1
    return res

def mat_pow_mod_as_list(mat, n, mod=0):
    """
    return (mat^n) % mod
    mat is defined as list, avoid overflow in numpy
    @mat(list)
    """

    if n < 0:
        raise ValueError("power n must be positive!")

    d = len(mat)
    res = [[0] * d for _ in range(d)]
    for i in range(d):
        res[i][i] = 1

    while n:
        if n & 1:
            res = dot_mod_as_list(res, mat, mod)
        mat = dot_mod_as_list(mat, mat, mod)
        n >>= 1
    return res

def mat_sum_pow_mod(A0, Q, n, mod=0):
    """
    return (A0 + Q A0 + Q^2 A0 + ... + Q^n A0) % mod
    """

    if n < 0:
        raise ValueError("power n must be positive!")
    if n == 0:
        return A0
    assert len(A0) == len(Q[0])

    if m:
        A0 = A0 % mod
        Q = Q % mod

    d = len(Q)
    O = np.zeros((d, d), dtype=np.int64)
    I = np.eye(d, dtype=np.int64)
    Q_ext = np.concatenate([np.concatenate([Q, I], axis=1), np.concatenate([O, I], axis=1)])
    Q_ext_pow = mat_pow_mod(Q_ext, n, mod)
    I2 = np.concatenate([I, I], axis=0)
    res = np.dot(Q_ext_pow, I2)
    if mod:
        res %= mod
    res = np.dot(res, A0)
    if mod:
        res %= mod
    return res[:d]