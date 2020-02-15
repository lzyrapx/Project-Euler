'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:29:08
'''
import time
import numpy as np
from base import _is_square, _factorial, fac
from base import _iroot, cprod, ggcd
from base import extended_gcd

print(_is_square(8))
print(_is_square(9))
print(_factorial(4))
print(fac(10))
print(_iroot(10,1)) # (10, True)
print(_iroot(10,2)) # (3, False)
print(cprod([2,3,5,6]))
print(cprod([2,3,5,7,11,12]))

print(ggcd([11,22,33]))

print(extended_gcd(3,6)) #  3*x + 6*y = gcd(3, 6) => (gcd(3,6), x, y) => (3,1,0)

from base import padic_base_p
print(padic_base_p(10,2)) # 1010
print(padic_base_p(10,3)) # 101

from base import iter_associate, sum_power_series_mod
def f(a, b):
    s = 1
    for i in range(1, a):
        s *= i
    print("a = {} b = {} s = {} ".format(a, b, s))
    return s
print(iter_associate(f, 4, 3)) # 120
print(sum_power_series_mod(1,100,19260817))

from base import rational_continous_frac
print(rational_continous_frac(415,93,10000)) # (4,2,6,7)

from base import continous_frac_convergent
print(continous_frac_convergent([4,2,6,7])) # [Fraction(4, 1), Fraction(9, 2), Fraction(58, 13), Fraction(415, 93)]
print(float(continous_frac_convergent([4,2,6,7])[2]))


from treeGeneration import pythagorean_triple_tree
print(pythagorean_triple_tree())
print(pythagorean_triple_tree((5,12,13),forward=False, trust=True))

from treeGeneration import stern_brocot_tree
res = stern_brocot_tree()
print("test...")
# while True:
#     g = next(res)
#     print(g)
#     if g[0] == 1 and g[1] == 4:
#         print(g)
#         break
#     time.sleep(1)
# print("stern_brocot_tree: ", next(res))

#########################################
from prime import _primes_list, _is_prime, _mobius_list, atkin_prime_sieve
p = _primes_list(100)
print(p)
p = atkin_prime_sieve(100)
print("atkin_prime_sieve = ", p)
Is_p = _is_prime(97)
print(Is_p)

mobius_list = _mobius_list(100)
print(len(mobius_list))
print(mobius_list[0])
print(mobius_list[100])

from prime import _pollard_rho,prime_divisor_decomposition, all_divisors
print(_pollard_rho(1000))
print(prime_divisor_decomposition(10**12+2000))
print(all_divisors(10 ** 9))

from prime import euler_phi, mobius, _largest_prime_factor_sieve, prime_counting
for i in range(1, 100):
    print(euler_phi(i), end=" ")
print()
for i in range(1, 100):
    print(mobius(i), end=" ")
print()
print(_largest_prime_factor_sieve(18))
print(_largest_prime_factor_sieve(25))
print(prime_counting(10 ** 5))

from prime import is_coprime
print("is_coprime: ", is_coprime(4, 5))
print("is_coprime: ", is_coprime(4, 2))

###############################################
from combinatoric import C, C_mod, permutations_multiset
from combinatoric import list_multiset_permutations
print(C(4,2))
print(C_mod(4,2,11))
print(C_mod(1000000, 123123, 10 ** 9 + 7))
print(permutations_multiset([1,2]))
print(permutations_multiset([1,1]))


print(permutations_multiset([1,1,1]))
for lists in list_multiset_permutations([1,2,3]):
    print(lists, end=" ")

from combinatoric import limited_combinations
from combinatoric import all_subsets, all_partitions, sequence_partitions

for lists in limited_combinations([[2, 2],[3, 4]]):
    print(lists, end=" ")

print("\nall_subsets testing...")
for lists in all_subsets([1, 2, 3, 4], 1, None):
    print(lists, end=" ")

print("\nall_partitions testing...")
for lists in all_partitions(3, 5):
    print(lists)

print(sequence_partitions([1, 2, 3], [1, 2]))
print(sequence_partitions([1, 2, 3, 4], [2, 2]))
print(sequence_partitions([2, 3, 4, 5, 6], [2, 3]))

##################################
from polynomial import poly_truncate, poly_add
print(poly_truncate(np.array([0,1,2,3,4,5,0]), "left")) # [1,2,3,4,5,0]
print(poly_truncate(np.array([0,1,2,3,4,5,0]), "right")) # [0,1,2,3,4,5]
print(poly_add(np.array([1,2,3,4,5]), np.array([5,4,3,2,1])))
print(poly_add(np.array([1,2,3,4,5]), np.array([1,2,3])))
from polynomial import poly_mul
print(poly_mul(np.array([1,2,3,4,5]),np.array([1,2,3])))


##################################
from linearAlgebra import dot_mod, dot_mod_as_list
from linearAlgebra import mat_pow_mod_as_list, mat_pow_mod
A = np.array([[1,2,3],[1,4,5]])
B = np.array(([2,3,4],[5,2,3],[0,0,0]))
print(dot_mod(A, B, 1997))
print(dot_mod(A, B, 10))
print("mat_pow_mod = ", mat_pow_mod(B, 10, 1997))
A = [[1,2,3],[1,4,5]]
B = [[2,3,4],[5,2,3]]
print(dot_mod_as_list(A, B, 1997))
print("mat_pow_mod = ", mat_pow_mod_as_list(B, 10, 1997))

from linearAlgebra import gauss_jordan_elimination, gauss_jordan_modular_elimination
A = [[3.0 ,2.0, 1.0, 6.0], [2.0, 2.0, 2.0, 4.0], [4.0, -2.0, -2.0, 2.0]] # ans = [1, 2, -1]
print(gauss_jordan_elimination(A))
A = [[3 ,2, 1, 6], [2, 2, 2, 4], [4, -2, -2, 2]] # ans = [1, 2, -1]
print(gauss_jordan_modular_elimination(A, 1997))
from linearAlgebra import get_integer_matrix_inverse_as_list
A = [[1, 2, 3], [2, 2, 1], [3, 4, 3]]
B = [[1, 2, 3], [3, 4, 6], [4, 6, 5]]
print(get_integer_matrix_inverse_as_list(A))
print(get_integer_matrix_inverse_as_list(B))




