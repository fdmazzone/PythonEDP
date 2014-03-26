# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:38:25 2013

@author: fernando
"""

from numpy import *
dx=0.1
nu=1
tmax=.1
dt=.5*dx**2
CondContDer=lambda t: zeros(t.shape)
CondContIzq=lambda t: zeros(t.shape)
CondValorInic=lambda x: zeros(x.shape)
TermNoHom=lambda X,Y: X**2+Y**2
CantSubintEsp=int(1/dx) 

 
CantSubintTiempo=int(tmax/dt) 
   
u=zeros([CantSubintEsp+1,CantSubintTiempo+1]) 
u[0:CantSubintEsp+1,0]=CondValorInic(linspace(0,1,CantSubintEsp+1))
u[0,:]=CondContIzq(linspace(0,tmax,CantSubintTiempo+1))
u[CantSubintEsp,:]=CondContDer(linspace(0,tmax,CantSubintTiempo+1)) 


[Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),linspace(0,1,\
CantSubintEsp+1)) 
NoHom=TermNoHom(X,Y) 
factor=2*nu*dt/(dx**2) 
u[1:CantSubintEsp-1,1]=u[1:CantSubintEsp-1,0]+factor/2*\
(u[2:CantSubintEsp,0]-2*u[1:CantSubintEsp-1,0]+u[0:CantSubintEsp-2,0])+\
dt*NoHom[1:CantSubintEsp-1,0] 
for n in range(2,CantSubintTiempo):
    u[1:CantSubintEsp-1,n]=u[1:CantSubintEsp-1,n-2]+factor*\
    (u[2:CantSubintEsp,n-1]-2*u[1:CantSubintEsp-1,n-1]+\
    u[0:CantSubintEsp-2,n-1])+dt*NoHom[1:CantSubintEsp-1,n-1] 
    
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from pylab import *

Z=u
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
plt.show()
savefig('test')    