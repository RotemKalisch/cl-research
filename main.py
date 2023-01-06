from measure_time import measure_time
from lattice import *
import zig_zag 
import e
import pi

import logging
import mpmath as mp
import numpy as np
from sympy.abc import x, y, k

def search_zig_zag():
    for n in range(3, 15, 2):
        print(n)
        try:
            constants = zig_zag.create_zetas(n, n+1)
            print(identify_rams(constants, zig_zag.VARIABLES, zig_zag.matrices(n), zig_zag.START, zig_zag.iterations(n)))
            print(check_symmetric_limit(zig_zag.VARIABLES, zig_zag.matrices(n)[0], 0, zig_zag.START, zig_zag.iterations(n)))
        except Exception as ex:
            print(str(ex))

def check_symmetry():
    print(identify_rams(pi.CONSTANTS, pi.VARIABLES, pi.MATRICES, pi.START))
    print(check_symmetric_limits(pi.VARIABLES, pi.MATRICES, pi.START))

def check_reduce(reduce):
    constants = zig_zag.create_zetas(3, 4)
    powered = mat_pow(zig_zag.VARIABLES, zig_zag.matrices(3)[0], step(0, 2), zig_zag.START, zig_zag.iterations(3), reduce)
    print(ram(powered), co_ram(powered))

if __name__ == '__main__':
    k = 1
    y = 1
    results = calc_zig_zag(k, y, zig_zag.create_zig_zag_dict(k))
    for (value, result) in results:
        print(value, result)
