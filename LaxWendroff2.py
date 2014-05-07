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
#A=array([[lambda u:-u, lambda ro: ro**0.5], [lambda ro: -ro, lambda u:-u]])
#Au=array([[-1, 0], ])
dx=.1
dt=.1
tmax=1
RoValorInic=lambda x:x
uValorInic=lambda x:x
A=lambda u, ro:array([[-u, ro**2], [-ro, -u]])
Au=lambda u, ro:array([[-ones(shape(ro)), zeros(shape(ro))], [zeros(shape(ro))\
, -ones(shape(ro))]])
Aro=lambda u, ro:array([[zeros(shape(ro)), 2*ro], [-ones(shape(ro)), \
zeros(shape(ro))]])


#def esquemalaxwendroff(dx,tmax,dt,RoValorInic,uValorInic,Ap,Au):
"""
RoValorInic,uValorInic son funciones de dos variables en los reales
A es una array con las funciones de variables u y ro escritas como cadena
"""    


nFila=int(1/dx)
x=linspace(0,1,nFila+1)
tt=linspace(0,1,int(tmax/dt)+1)
V=np.zeros((size(tt),2,size(x)))
V[0,0,:]=uValorInic(x)
V[0,1,:]=RoValorInic(x)
#vect=np.zeros([size(x), 2, 1])
V0= V[0,:,:]
AA=A(V[0,0,:],V[0,1,:])
AAu=Au(V[0,0,:],V[0,1,:])
AAro=Aro(V[0,0,:],V[0,1,:])
Vx=zeros(shape(V0[:,1:-1]))
Vxx=Vx
for t in range(0,int(tmax/dt)):
    Vx=(V0[:,2:nFila+1]-V0[:,0:nFila-1])/(2*dx)
    Vxx=(V0[:,2:nFila+1]-2*V0[:,1:nFila]+V0[:,0:nFila-1])/(dx**2)
        
    Vt[1:nFila,:]=dot(AA,Vx)
    
    vect[:,0,0]=Vx[0,:,0]
    vect[:,1,0]=Vx[0,:,1]
    Mx=Mu*vect[:,0,0]+Mro*vect[:,1,0]
    vect[:,0,0]=Vt[0,:,0]
    vect[:,1,0]=Vt[0,:,1]
    Mt=Mu*vect[:,0,0]+Mro*vect[:,1,0]
    
    Vxt[t,1:nFila,:]=Mx[1:nFila,:,:]*Vx[t,1:nFila,:]+\
        M[1:nFila,:,:]*Vx[t,1:nFila,:]
    Vtt[t,1:nFila,:]=Mt[1:nFila,:,:]*Vx[t,1:nFila,:]+\
        M[1:nFila,:,:]*Vxt[t,1:nFila,:]
    V[t+1,1:nFila,:]=V[t,1:nFila,:]+Vt[t,1:nFila,:]*dt+\
        Vtt[t,1:nFila,:]*(dt**2)/2
        
    V0= V[t+1,:,:]
    AA=A(V[t+1,:,0],V[t+1,:,1])
    AAu=Au(V[t+1,:,0],V[t+1,:,1])
    AAro=Aro(V[t+1,:,0],V[t+1,:,1])