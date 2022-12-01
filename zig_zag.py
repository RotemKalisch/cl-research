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

Mx = sp.Matrix([[(x ** k + (x + 1) ** k + y * ((x + 1) ** (k - 1))), ((-x ** (k - 1)) * (x + y))], [(x + 1) ** k, 0]]) # original matrix given by CL document
Mx_tag = sp.Matrix([[2 * (x * ((x-1)**k + x**k + x**(k-1))), 1], [-4*(x+1)**2 * x**(2*k), 0]]) # A nice permutation that gets to 1/zig_zag(k)
Mx_tag2 = sp.Matrix([[(x+1)*(x)**k + (x+1)**k + (x+1)**(k+1), -(x)**(2*k) * (x+1)**2], [1, 0]]) # Second attempt at ofir's research

"""
The following matrices only work for k=1
"""
def a(x, y, k):
    return 2*(2 + 4*x + 2*x**2 - y + y**2)
def b(x, k):
    return -4*(x**2)*(x + 1)**2
def f(x, y, k):
    return 2*x**2 + 2*x*y + y**2 + 2*x + y

Mx_1 = sp.Matrix([[a(x, y, k), b(x,k)], [1, 0]])
My_1 = sp.Matrix([[f(x, y, k), b(x,k)], [1, -f(x, -y, k)]])
print(Mx_1)
print(My_1)

VARIABLES = [x, y]
MATRICES = [Mx_1, My_1]
START = np.array([1, 1])

def matrices(n):
    return [matrix.subs([[k, n]]) for matrix in MATRICES]
