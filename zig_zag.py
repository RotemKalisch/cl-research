import itertools
import sympy as sp
import numpy as np
import mpmath as mp
from sympy.abc import x, y, k

def zig_zag(k):
    mp.mp.dps = 50
    retval = 0
    for i in range(k, 1, -1):
        retval += ((-1) ** (k - i)) * mp.mpf(mp.zeta(i))
    retval += (-1) ** (k-1)
    return retval

def create_zig_zag_dict(k):
    return {'zig_zag({})'.format(k): zig_zag(k)}


Mx = sp.Matrix([[(x ** k + (x + 1) ** k + y * ((x + 1) ** (k - 1))), ((-x ** (k - 1)) * (x + y))], [(x + 1) ** k, 0]])

My = sp.Matrix([[(x ** k + (x + 1) ** k), -(x ** k + (x + 1) ** k)], [x ** (k + 1), -x ** (k + 1)]])

VARIABLES = [x, y]
MATRICES = [Mx]
START = np.array([1, 0])

def matrices(n):
    return [Mx.subs([[k, n]])]

def iterations(k):
    if k <= 3:
        return 10000
    if k <= 5:
        return 1000
    return 200
