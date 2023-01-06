from measure_time import measure_time

import logging
import mpmath as mp
import numpy as np
from sympy.abc import x, y, k

from lattice import *
from known_lattices import *

if __name__ == '__main__':
    print(eval(repr(degree_high_degenerated(7))))
