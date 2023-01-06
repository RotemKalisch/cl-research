import itertools
import sympy as sp
import numpy as np
import math

from sympy.abc import x, y

def a(f, fbar):
    return f - fbar.subs(x, x + 1)

def b(f, fbar):
    ff = f*fbar
    return ff.subs([[y, 0]]) - ff.subs([[x, 0], [y, 0]])

def Mx(f, fbar):
    return sp.Matrix([[0, b(f,fbar)], [1, a(f,fbar)]]);

def My(f, fbar):
    return sp.Matrix([[fbar, b(f,fbar)], [1, f]]);

class Lattice:
    def __init__(self, f, fbar):
        self.f = f
        self.fbar = fbar
        self.Mx = Mx(self.f, self.fbar)
        self.My = My(self.f, self.fbar)

    def __repr__(self):
        return "Lattice({},{})".format(self.f, self.fbar)

    def __str__(self):
        return str([self.Mx, self.My])

    def Mx(self):
        return self.Mx

    def My(self):
        return self.My

    def subs(substitutions):
        return Lattice(f.subs(substitutions), fbar.subs(substitutions))
