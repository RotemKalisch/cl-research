import itertools
import sympy
import math
import numpy

from sympy.abc import x, y

Mx = sympy.Matrix([[2*x + y + 1, -x], [-y, x]])
My = sympy.Matrix([[x + 2*y + 1, -x], [-1 - y, 1 + y]])

VARIABLES = [x, y]
MATRICES = [Mx, My]
START = numpy.array([1, 0])
CONSTANTS = {'pi': math.pi}
