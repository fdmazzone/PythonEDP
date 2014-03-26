# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 19:20:06 2013

@author: Noelia
"""

from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

CantSubintEsp=10.0
dx=1/CantSubintEsp
CantIteracTiempo=10
#dt=dx**2/16
dt=dx**2/2**CantIteracTiempo
nu=1
ValorInic=lambda x: x**2
TermNoHom=lambda x,t:2*t-2
ContIzq=lambda t:t**2
ContDer=lambda t:1+t**2

tmax=2
u=zeros([CantSubintEsp+1,int(tmax/dt)+1])
u[:,0]=ValorInic(linspace(0,1,CantSubintEsp+1))
u[0,:]=ContIzq(linspace(0,tmax,tmax/dt+1))
u[CantSubintEsp,:]=ContDer(linspace(0,tmax,tmax/dt+1))
[Y,X]=meshgrid(linspace(0,tmax,tmax/dt+1),\
    linspace(0,1,CantSubintEsp+1))
NoHom=TermNoHom(X,Y)
factor=2*nu*dt/(dx**2)
u[1:CantSubintEsp,1]=u[1:CantSubintEsp,0]+factor/2*\
    (u[2:CantSubintEsp+1,0]-2*u[1:CantSubintEsp,0]+u[0:CantSubintEsp-1,0])\
    +dt*NoHom[1:CantSubintEsp,0]
for n in range(1,int(tmax/dt)):
    u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n-1]\
        +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
        +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
    linewidth=0, antialiased=False)
plt.show()
