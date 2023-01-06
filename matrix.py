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


