# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 10:12:38 2014

@author: Noelia
"""
from numpy import *
#from sympy import *

#p,u=symbols('p,u')
#f=p**0.5
#A = matrix([[-u, f.diff(p)],[-p, -u]])
#Ap=lambdify([u,p],A.diff(p),'numpy')
#Au=lambdify([u,p],A.diff(u),'numpy')
#A=array([[lambdify([u, ro], A[0,0], "numpy"), lambdify([u, ro], A[0,1], "numpy")], \
#              [lambdify([u, ro], A[1,0], "numpy"), lambdify([u, ro], A[1,1], "numpy")]])
A=array([[lambda u:-u, lambda ro: ro**0.5], [lambda ro: -ro, lambda u:-u]])
#Au=array([[-1, 0], ])
dx=.1
dt=.1
tmax=1
RoValorInic=lambda x:x
uValorInic=lambda x:x
A=lambda u, ro:array([[-u, ro**2], [-ro, -u]])
Au=lambda u, ro:array([[-1, 0], [0, -1]])
Aro=lambda u, ro:array([[0, 2*ro], [-1, 0]])


#def esquemalaxwendroff(dx,tmax,dt,RoValorInic,uValorInic,Ap,Au):
"""
RoValorInic,uValorInic son funciones de dos variables en los reales
A es una array con las funciones de variables u y ro escritas como cadena
"""    

x=linspace(0,1,int(1/dx)+1)
t=linspace(0,1,int(tmax/dt)+1)
   

V=np.zeros((size(t),2,size(x)))
V[0,0,:]=uValorInic(x)
V[0,1,:]=RoValorInic(x)
#vect=np.zeros([size(x), 2, 1])
V0= V[0,:,:]
AA=A(V[0,0,:],V[0,1,:])
AAu=Au(V[0,0,:],V[0,1,:])
AAro=Aro(V[0,0,:],V[0,1,:])

for t in range(0,int(tmax/dt)):
    Vx[1:int(1/dx),:]=(V0[2:int(1/dx)+1,:]-V0[0:int(1/dx)-1])/(2*dx)
    Vxx[1:int(1/dx),:]=(V0[2:int(1/dx)+1,:]-2*V0[1:int(1/dx),:]+\
        V0[0:int(1/dx)-1])/(dx**2)
        
    Vt[1:int(1/dx),:]=dot(A,Vx[1:int(1/dx),:])
    
    vect[:,0,0]=Vx[0,:,0]
    vect[:,1,0]=Vx[0,:,1]
    Mx=Mu*vect[:,0,0]+Mro*vect[:,1,0]
    vect[:,0,0]=Vt[0,:,0]
    vect[:,1,0]=Vt[0,:,1]
    Mt=Mu*vect[:,0,0]+Mro*vect[:,1,0]
    
    Vxt[t,1:int(1/dx),:]=Mx[1:int(1/dx),:,:]*Vx[t,1:int(1/dx),:]+\
        M[1:int(1/dx),:,:]*Vx[t,1:int(1/dx),:]
    Vtt[t,1:int(1/dx),:]=Mt[1:int(1/dx),:,:]*Vx[t,1:int(1/dx),:]+\
        M[1:int(1/dx),:,:]*Vxt[t,1:int(1/dx),:]
    V[t+1,1:int(1/dx),:]=V[t,1:int(1/dx),:]+Vt[t,1:int(1/dx),:]*dt+\
        Vtt[t,1:int(1/dx),:]*(dt**2)/2
        
    V0= V[t+1,:,:]
    AA=A(V[t+1,:,0],V[t+1,:,1])
    AAu=Au(V[t+1,:,0],V[t+1,:,1])
    AAro=Aro(V[t+1,:,0],V[t+1,:,1])