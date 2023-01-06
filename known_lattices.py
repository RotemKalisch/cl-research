import itertools
import sympy as sp
import numpy as np
import math

from sympy.abc import x, y

c0 = sp.Symbol('c0')
c1 = sp.Symbol('c1')
c2 = sp.Symbol('c2')
c3 = sp.Symbol('c3')

def degree1():
    f = c0 + c1*(x+y)
    fbar = c2 + c3*(x-y)
    return (f, fbar)
    
def degree2():
    f = ((2*c1+c2)*(c1+c2)-c3*c0) + c3*((2*c1+c2)(x+y)+(c1+c2)*(2*x+y)) + c3**2 * (2*x**2+2*x*y+y**2)
    fbar = c3*(c0+c2*x+c1+y-c3*(2*x**2-2*x*y+y**2))
    return (f,fbar)

def degree_high_degenerated(d):
    f = -c0*(x**d) + c1*(y**d)
    fbar = c0*(x**d) + c1*(y**d)
    return (f,fbar)
