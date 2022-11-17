import itertools
import sympy
import numpy

from sympy.abc import x, y, k
from scipy.special import zeta

def create_zetas_dict(start, end):
    zetas = {}
    for i in range(start, end):
        zetas['zeta({})'.format(i)] = zeta(i)
    return zetas

def create_zig_zag(n):
    retval = 0
    for i in range(n, 1, -1):
        retval += ((-1) ** (n - i)) * zeta(i)
    retval += (-1) ** (n-1)
    return retval
    return {'zig_zag({})'.format(n): retval}

def create_zig_zag_dict(n):
    return {'zig_zag({})'.format(n): create_zig_zag(n)}


Mx = sympy.Matrix([[(x ** k + (x + 1) ** k + y * ((x + 1) ** (k - 1))), ((-x ** (k - 1)) * (x + y))], [(x + 1) ** k, 0]])

My = sympy.Matrix([[(x ** k + (x + 1) ** k), -(x ** k + (x + 1) ** k)], [x ** (k + 1), -x ** (k + 1)]])

VARIABLES = [x, y]
MATRICES = [Mx]
START = numpy.array([1, 0])

def matrices(n):
    return [Mx.subs([[k, n]])]

def iterations(k):
    if k <= 3:
        return 10000
    if k <= 5:
        return 1000
    return 200
