# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/usuario/.spyder2/.temp.py
"""
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

"""  
Este programa resuelve la ecuacin del calor Vt=nu*Vxx+TermNoHom en el 
intervalo
[0,1] con condiciones de contorno de Neumann (mixto)
dx=espaciado en x 
nu=factor de la ecuacin
TermNoHom=termino no homogeneo de tipo funcn que acepte matrices
ValorIni=valor inicial (t=0) de tipo funcin que acepte vectores
CondContIzq, COndContDer=valores de contorno de tipo funcin que acepte
vectores
tmax=valor maximo de t
dt=espaciado en t (para obtener una mejor aproximacin de la solucin es 
conveniente tomar dt<0.5*dx^2)
"""    
dx=.1
nu=1 
TermNoHom=lambda X,Y: X**2+Y**2
CondValorInic=lambda X: zeros(X.shape)
ValorInicDeriv=lambda T: zeros(T.shape)
CondContDer=lambda T: zeros(T.shape)
tmax=1
dt=.5*dx**2


CantSubintEsp=int(1/dx)
CantSubintTiempo=int(tmax/dt)
u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
u[:,0]=CondValorInic(linspace(0,1,CantSubintEsp+1))
ValorDeriv=ValorInicDeriv(linspace(0,tmax,CantSubintTiempo+1))
u[CantSubintEsp,:]=CondContDer(linspace(0,tmax,CantSubintTiempo+1))
[Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),linspace(0,1,CantSubintEsp+1))
print 2
NoHom=TermNoHom(X,Y)
factor=nu*dt/(dx**2)
for n in range(CantSubintTiempo):
    u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n]\
        +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
        +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
    u[0,n+1]=u[1,n+1]-dx*ValorDeriv[n]
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
    linewidth=0, antialiased=False)
plt.show()
savefig('test')
