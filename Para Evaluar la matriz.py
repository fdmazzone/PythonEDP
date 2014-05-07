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


>>> f=lambda u, ro:array([[u*ro, u*ro], [u*ro, u*ro]])
>>> f(x,x)
array([[[ 0.  ,  0.01,  0.04,  0.09,  0.16,  0.25,  0.36,  0.49,  0.64,
          0.81,  1.  ],
        [ 0.  ,  0.01,  0.04,  0.09,  0.16,  0.25,  0.36,  0.49,  0.64,
          0.81,  1.  ]],

       [[ 0.  ,  0.01,  0.04,  0.09,  0.16,  0.25,  0.36,  0.49,  0.64,
          0.81,  1.  ],
        [ 0.  ,  0.01,  0.04,  0.09,  0.16,  0.25,  0.36,  0.49,  0.64,
          0.81,  1.  ]]])
>>> borrar=f(x,x)
>>> borrar[:,:,0]
array([[ 0.,  0.],
       [ 0.,  0.]])
