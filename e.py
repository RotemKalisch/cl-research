import itertools
import sympy as sp
import numpy as np
import math

from sympy.abc import x, y

Mx = sp.Matrix([[x + y + 2, -y - 1], [-1, 1]])
My = sp.Matrix([[x + y + 1, -y - 1], [-1, 0]])

VARIABLES = [x, y]
MATRICES = [Mx, My]
START = np.zeros(len(VARIABLES))
CONSTANTS = {'e': math.e}
