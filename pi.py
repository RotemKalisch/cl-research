import itertools
import sympy as sp
import numpy as np
import math

from sympy.abc import x, y

Mx = sp.Matrix([[2*x + y + 1, -x], [-y, x]])
My = sp.Matrix([[x + 2*y + 1, -x], [-1 - y, 1 + y]])

VARIABLES = [x, y]
MATRICES = [Mx, My]
START = np.array([1, 0])
CONSTANTS = {'pi': math.pi}
