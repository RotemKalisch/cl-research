import itertools
import sympy
import numpy
import math

from sympy.abc import x, y

Mx = sympy.Matrix([[x + y + 2, -y - 1], [-1, 1]])
My = sympy.Matrix([[x + y + 1, -y - 1], [-1, 0]])

VARIABLES = [x, y]
MATRICES = [Mx, My]
START = numpy.zeros(len(VARIABLES))
CONSTANTS = {'e': math.e}
