# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 19:49:22 2013

@author: Noelia
"""
from numpy import *
from pylab import *

CantSubintEsp=10.0
dx=1/CantSubintEsp
CantIteracTiempo=1000
#dt=(dx**2)/(2**CantIteracTiempo)
dt=dx**2/128

mu=4/dx**2*sin(pi*dx/2)**2
k=dt*mu
c1=(1-1/sqrt(k**2+1))/2
c2=(1+1/sqrt(k**2+1))/2
T=zeros(CantIteracTiempo)
T[0]=1
T[1]=1-k
for m in range(2,CantIteracTiempo):
    T[m]=c1*((-k+sqrt(k**2+1))**m)+c2*((-k-sqrt(k**2+1))**m)
x=linspace(0,CantIteracTiempo,CantIteracTiempo)
figure()
plot(x,T)
show()