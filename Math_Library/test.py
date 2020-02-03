'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:29:08
'''
import time
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
while True:
    g = next(res)
    print(g)
    if g[0] == 1 and g[1] == 4:
        print(g)
        break
    time.sleep(1)
print("stern_brocot_tree: ", next(res))
