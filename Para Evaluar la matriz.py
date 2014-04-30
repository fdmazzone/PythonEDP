# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 10:07:48 2014

@author: fernando
"""
from sympy import *
p,u=symbols('p,u')
f=p**0.5
A = Matrix([[-u, f.diff(p)],[-p, -u]])
Ap=lambdify([u,p],A.diff(p),'numpy')
Au=lambdify([u,p],A.diff(u),'numpy')
