# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:56:57 2013

@author: Stefania Demaria y Noelia Cardarelli
"""
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def esquemadesaltos(dx,nu,TermNoHom,ValorInic,ContIzq,ContDer,tmax,dt):
    
    """  
    Este programa resuelve la ecuacin del calor Vt=nu*Vxx+TermNoHom en el 
    intervalo
    [0,1] con condiciones de contorno de Dirichlet
    dx=espaciado en x
    nu=factor de la ecuacin
    TermNoHom=termino no homogeneo de tipo funcn que acepte matrices
    ValorIni=valor inicial (t=0) de tipo funcin que acepte vectores
    CondContIzq, COndContDer=valores de contorno de tipo funcin que acepte
    vectores
    tmax=valor maximo de t
    dt=espaciado en t (para obtener una mejor aproximacin de la solucin es 
    conveniente tomar dt<0.5*dx^2)
    --- control de datos ---
    """
    CantSubintEsp=int(1/dx)
    CantSubintTiempo=int(tmax/dt)
    u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
    u[:,0]=ValorInic(linspace(0,1,CantSubintEsp+1))
    u[0,:]=ContIzq(linspace(0,tmax,CantSubintTiempo+1))
    u[CantSubintEsp,:]=ContDer(linspace(0,tmax,CantSubintTiempo+1))
    [Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),\
        linspace(0,1,CantSubintEsp+1))
    NoHom=TermNoHom(X,Y)
    factor=2*nu*dt/(dx**2)
    u[1:CantSubintEsp,1]=u[1:CantSubintEsp,0]+factor/2*\
        (u[2:CantSubintEsp+1,0]-2*u[1:CantSubintEsp,0]+u[0:CantSubintEsp-1,0])\
        +dt*NoHom[1:CantSubintEsp,0]
    for n in range(1,CantSubintTiempo):
        u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n-1]\
            +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
            +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
        linewidth=0, antialiased=False)
    plt.show()
    return X,Y,u
    
    
def esquemaexplicito(dx,nu,TermNoHom,ValorInic,ContIzq,ContDer,tmax,dt):
    
    """  
    Este programa resuelve la ecuacin del calor Vt=nu*Vxx+TermNoHom en el 
    intervalo
    [0,1] con condiciones de contorno de Dirichlet
    dx=espaciado en x
    nu=factor de la ecuacin
    TermNoHom=termino no homogeneo de tipo funcn que acepte matrices
    ValorIni=valor inicial (t=0) de tipo funcin que acepte vectores
    CondContIzq, COndContDer=valores de contorno de tipo funcin que acepte
    vectores
    tmax=valor maximo de t
    dt=espaciado en t (para obtener una mejor aproximacin de la solucin es 
    conveniente tomar dt<0.5*dx^2)
    --- control de datos ---
    """
    CantSubintEsp=int(1/dx)
    CantSubintTiempo=int(tmax/dt)
    u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
    u[:,0]=ValorInic(linspace(0,1,CantSubintEsp+1))
    u[0,:]=ContIzq(linspace(0,tmax,CantSubintTiempo+1))
    u[CantSubintEsp,:]=ContDer(linspace(0,tmax,CantSubintTiempo+1))
    [Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),\
        linspace(0,1,CantSubintEsp+1))
    NoHom=TermNoHom(X,Y)
    factor=nu*dt/(dx**2)
    for n in range(CantSubintTiempo):
        u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n]\
            +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
            +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
        linewidth=0, antialiased=False)
    plt.show()
    return X,Y,u
    
    
def esquemaexplicito2orden(dx,nu,TermNoHom,ValorInic,ContIzq,ContDer,tmax):
    
    """  
    Este programa resuelve la ecuacin del calor Vt=nu*Vxx+TermNoHom en el 
    intervalo
    [0,1] con condiciones de contorno de Dirichlet
    dx=espaciado en x
    nu=factor de la ecuacin
    TermNoHom=termino no homogeneo de tipo funcn que acepte matrices
    ValorIni=valor inicial (t=0) de tipo funcin que acepte vectores
    CondContIzq, COndContDer=valores de contorno de tipo funcin que acepte
    vectores
    tmax=valor maximo de t
    dt=espaciado en t (para obtener una mejor aproximacin de la solucin es 
    conveniente tomar dt<0.5*dx^2)
    --- control de datos ---
    """
    CantSubintEsp=int(1/dx)
    dt=dx**2/(6*nu)
    CantSubintTiempo=int(tmax/dt)
    u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
    u[:,0]=ValorInic(linspace(0,1,CantSubintEsp+1))
    u[0,:]=ContIzq(linspace(0,tmax,CantSubintTiempo+1))
    u[CantSubintEsp,:]=ContDer(linspace(0,tmax,CantSubintTiempo+1))
    [Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),\
        linspace(0,1,CantSubintEsp+1))
    NoHom=TermNoHom(X,Y)
    factor=nu*dt/(dx**2)
    for n in range(CantSubintTiempo):
        u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n]\
            +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
            +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
        linewidth=0, antialiased=False)
    plt.show()
    return X,Y,u      


def esquemaneumann(dx,nu,TermNoHom,CondValorInic,ValorInicDeriv,CondContDer,\
tmax,dt):
    CantSubintEsp=int(1/dx)
    CantSubintTiempo=int(tmax/dt)
    u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
    u[:,0]=CondValorInic(linspace(0,1,CantSubintEsp+1))
    ValorDeriv=ValorInicDeriv(linspace(0,tmax,CantSubintTiempo+1))
    u[CantSubintEsp,:]=CondContDer(linspace(0,tmax,CantSubintTiempo+1))
    [Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),linspace(0,1,CantSubintEsp+1))
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
    return X,Y,u

def esquemaneumann2orden(dx,nu,TermNoHom,CondValorInic,ValorInicDeriv,\
CondContDer,tmax,dt):
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
    CantSubintEsp=int(1/dx)
    CantSubintTiempo=int(tmax/dt)
    u=zeros([CantSubintEsp+1,CantSubintTiempo+1])
    u[:,0]=CondValorInic(linspace(0,1,CantSubintEsp+1))
    ValorDeriv=ValorInicDeriv(linspace(0,tmax,CantSubintTiempo+1))
    u[CantSubintEsp,:]=CondContDer(linspace(0,tmax,CantSubintTiempo+1))
    [Y,X]=meshgrid(linspace(0,tmax,CantSubintTiempo+1),\
        linspace(0,1,CantSubintEsp+1))
    NoHom=TermNoHom(X,Y)
    factor=nu*dt/(dx**2)
    factor2=2*factor
    factor3=2*nu*dt/(dx)
    for n in range(CantSubintTiempo):
        u[1:CantSubintEsp,n+1]=u[1:CantSubintEsp,n]\
            +factor*(u[2:CantSubintEsp+1,n]-2*u[1:CantSubintEsp,n]\
            +u[0:CantSubintEsp-1,n])+dt*NoHom[1:CantSubintEsp,n]
        u[0,n+1]=u[0,n]+factor2*(u[1,n]-u[0,n])-factor3*ValorDeriv[n];
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,\
        linewidth=0, antialiased=False)
    plt.show()
    return X,Y,u 