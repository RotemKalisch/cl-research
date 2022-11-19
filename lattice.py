import itertools
import math
import logging
import mpmath as mp
import numpy as np
import sympy as sp

EPSILON = 1e-5
DEFAULT_ITERATIONS=200

def dimension(variables):
    return len(variables)

def step(index, n):
    retval = np.zeros(n, dtype=int)
    retval[index] = 1
    return retval

def substitutions(variables, values):
    subs = []
    for i in range(len(values)):
        subs.append([variables[i], values[i]])
    return subs

def matrix_gcd(matrix):
    return math.gcd(*matrix)

def remove_gcd(matrix):
    gcd = matrix_gcd(matrix)
    for i in range(len(matrix)):
        matrix[i] //= gcd

# Make sure that len(start) = len(step) = amount of vars in matrix
def mat_pow(variables, matrix, step, start, iterations=DEFAULT_ITERATIONS, reduce=True, initial_matrix=sp.eye(2)):
    curr = np.copy(start)
    retval = initial_matrix
    while iterations > 0:
        retval = matrix.subs(substitutions(variables, curr)) * retval
        iterations -= 1
        curr += step
        if reduce == True and iterations % 200 == True:
            remove_gcd(retval)
    if reduce:
        remove_gcd(retval)
    return retval


def mat_pow_precise(variables, matrix, step, start, dps=mp.mp.dps, initial_matrix=sp.eye(2)):
    """
    Function multiples matrix by itself until the difference between consequent ram values is smaller than dps.
    """
    logger = logging.getLogger('lattice.mat_pow_precise')
    prev_dps = mp.mp.dps
    mp.mp.dps = 2 * dps
    curr = np.copy(start)
    retval = initial_matrix
    prev_ram = mp.mpf(0)
    next_ram = mp.mpf(math.inf)
    iterations = 0
    while abs(next_ram - prev_ram) >= 10 ** (-1 * dps):
        retval = matrix.subs(substitutions(variables, curr)) * retval
        curr += step
        prev_ram = next_ram
        next_ram = ram(retval)
        iterations += 1
        if iterations % 1000 == 0:
            logger.debug("mat_pow_precise: iteration {}, diff={}".format(iterations, abs(next_ram - prev_ram)))
    # restore dps
    mp.mp.dps = prev_dps
    logger.debug("mat_pow_precise: finished after {} iterations".format(iterations))
    return retval


def assert_ram(variables, matrix, step, start, limit, dps=mp.mp.dps, initial_matrix=sp.eye(2)):
    """
    Function multiples matrix by itself until |limit - ram| < dps
    """
    logger = logging.getLogger('lattice.assert_ram')
    prev_dps = mp.mp.dps
    mp.mp.dps = 2 * dps
    curr = np.copy(start)
    retval = initial_matrix
    ram_value = mp.mpf(math.inf)
    iterations = 0
    while abs(limit - ram_value) >= 10 ** (-1 * dps):
        retval = matrix.subs(substitutions(variables, curr)) * retval
        curr += step
        ram_value = ram(retval)
        iterations += 1
        if iterations % 1000 == 0:
            logger.debug("assert_ram: iteration {}, diff={}".format(iterations, abs(limit - ram_value)))
    # restore dps
    mp.mp.dps = prev_dps
    logger.debug("assert_ram: finished after {} iterations".format(iterations))
    return retval, ram_value


def ram(matrix):
    if mp.mpf(matrix[1]) == mp.mpf(0):
        return math.inf
    return mp.mpf(matrix[0]) / mp.mpf(matrix[1])

def co_ram(matrix):
    if mp.mpf(matrix[3]) == mp.mpf(0):
        return math.inf
    return mp.mpf(matrix[2]) / mp.mpf(matrix[3])

def check_symmetry(matrix):
    ram_val = ram(matrix)
    co_ram_val = co_ram(matrix)
    if abs(ram(matrix) - co_ram(matrix)) < EPSILON:
        return True
    print("Ram:{}, CoRam:{}".format(ram_val, co_ram_val))
    return False

def check_symmetric_limit(variables, matrix, index, start=None, iterations=DEFAULT_ITERATIONS):
    powered = mat_pow(variables, matrix, step(index, dimension(variables)), start, iterations)
    return check_symmetry(powered)


def check_symmetric_limits(variables, matrices, start=None, iterations=DEFAULT_ITERATIONS):
    return [check_symmetric_limit(variables, matrices[i], i, start, iterations) for i in range(dimension(variables))]

def evaluate_ram(variables, matrix, step, start=None, iterations=DEFAULT_ITERATIONS):
    powered = mat_pow(variables, matrix, step, start, iterations)
    return ram(powered)

def format_couple(a, b, constant_str):
    if b == 0:
        return "{}".format(a)
    if a == 0:
        if b == 1:
            return "{}".format(constant_str)
        return "({}*{})".format(b, constant_str)
    sign = '-' if b < 0 else '+'
    return "({} {} {}*{})".format(a, sign, abs(b), constant_str)

def format_mobius(eq, constant_str):
    numerator = format_couple(eq[0], eq[1], constant_str)
    denominator = format_couple(-1*eq[2], -1*eq[3], constant_str)
    return "{} / {}".format(numerator, denominator)

def identify_mobius(ram_value, constant, tol):
    """
    identifies a relation such that ram_value = mobius(constant)
    """
    assert len(constant) == 1
    constant_value = list(constant.values())[0]
    eq = mp.pslq([1, constant_value, ram_value, constant_value * ram_value], maxcoeff=10000, tol=tol)
    if eq is None:
        return None
    constant_str = list(constant.keys())[0]
    return format_mobius(eq, constant_str)
    
def identify_ram(constants, variables, matrix, step, start=None, iterations=DEFAULT_ITERATIONS):
    ram_value = evaluate_ram(variables, matrix, step, start, iterations)
    return mp.identify(ram_value, constants)

def identify_rams(constants, variables, matrices, start=None, iterations=DEFAULT_ITERATIONS):
    dim = min(dimension(variables), len(matrices))
    return [identify_ram(constants, variables, matrices[i], step(i, dimension(variables)), start, iterations) for i in range(dim)]
    
