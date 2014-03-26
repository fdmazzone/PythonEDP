# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 22:05:02 2013
Modificado 26-3-2014
Modificado 26-3-2014b

@author: usuario
"""
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from pylab import *

dx=0.07
nu=1
tmax=1
dt=0.5*dx**2
TermNoHom=lambda X,Y: X**2+Y**2
ValorIni=lambda x: zeros(x.shape)
ContIzq=lambda t: zeros(t.shape)
ContDer=lambda t: zeros(t.shape)
CantSubintEsp=int(1/dx)
CantSubintTiempo=int(tmax/dt)
u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
u[:,0]=ValorIni(linspace(0,1,CantSubintEsp+1))
u[0,:]=ContIzq(linspace(0,tmax,CantSubintTiempo+1))
u[CantSubintEsp,:]=ContDer(linspace(0,tmax,CantSubintTiempo+1))
[Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),\
linspace(0,1,CantSubintEsp+1))
NoHom=TermNoHom(X,Y);
factor=2*nu*dt/(dx**2);
u[1:CantSubintEsp-1,1]=u[1:CantSubintEsp-1,0]+factor/2*\
(u[2:CantSubintEsp,0]-2*u[1:CantSubintEsp-1,0]+u[0:CantSubintEsp-2,0])\
+dt*NoHom[1:CantSubintEsp-1,0]
for n in range(1,CantSubintTiempo):
    u[1:CantSubintEsp-1,n+1]=u[1:CantSubintEsp-1,n-1]\
    +factor*(u[2:CantSubintEsp,n]-2*u[1:CantSubintEsp-1,n]\
    +u[0:CantSubintEsp-2,n])+dt*NoHom[1:CantSubintEsp-1,n]
Z=u
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
plt.show()
savefig('test')
