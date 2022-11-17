import itertools
import math
import mpmath
import numpy
import sympy

EPSILON = 1e-5
DEFAULT_ITERATIONS=200

def dimension(variables):
    return len(variables)

def origin(n):
    return numpy.zeros(n, dtype=int)

def step(index, n):
    retval = origin(n)
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
def mat_pow(variables, matrix, step, start, iterations=DEFAULT_ITERATIONS, reduce=True, initial_matrix=sympy.eye(2)):
    if start is None:
        start = origin(dimension(variables))
    curr = numpy.copy(start)
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

def ram(matrix):
    retval = matrix[0] / matrix[1]
    #print(matrix)
    #print(retval)
    return float(retval)

def co_ram(matrix):
    retval = matrix[2] / matrix[3]
    #print(matrix)
    return float(retval)

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
    if a == 0:
        return "{}*{}".format(b, constant_str)
    sign = '-' if b < 0 else '+'
    return "{} {} {}*{}".format(a, sign, abs(b), constant_str)

def format_mobius(eq, constant_str):
    numerator = format_couple(eq[0], eq[1], constant_str)
    denominator = format_couple(eq[2], eq[3], constant_str)
    return "{} / ({})".format(numerator, denominator)

def identify_mobius(ram_value, constant):
    """
    identifies a relation such that ram_value = mobius(constant)
    """
    assert len(constant) == 1
    constant_value = list(constant.values())[0]
    eq = mpmath.pslq([1, constant_value, ram_value, constant_value * ram_value])
    if eq is None:
        return None
    constant_str = list(constant.keys())[0]
    return format_mobius(eq, constant_str)
    
def identify_ram(constants, variables, matrix, step, start=None, iterations=DEFAULT_ITERATIONS):
    ram_value = evaluate_ram(variables, matrix, step, start, iterations)
    return mpmath.identify(ram_value, constants)

def identify_rams(constants, variables, matrices, start=None, iterations=DEFAULT_ITERATIONS):
    dim = min(dimension(variables), len(matrices))
    return [identify_ram(constants, variables, matrices[i], step(i, dimension(variables)), start, iterations) for i in range(dim)]
    
