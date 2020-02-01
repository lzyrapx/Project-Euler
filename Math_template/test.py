'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:29:08
'''
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
