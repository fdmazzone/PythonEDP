# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 18:38:53 2014

@author: Noelia
"""
"""GRAFICO 2D"""

from matplotlib import cm
import matplotlib.pyplot as plt

X,Y=meshgrid(x,tt)

fig, ax= subplots()
#p=ax.pcolor(X,Y,V[:,0,:], cmap=cm.RdBu, vmin=abs(V[:,0,:]).min(), vmax=abs(V[:,0,:]).max())
p=ax.pcolor(X,Y,V[:,0,:])
cb=fig.colorbar(p)

solucion=np.zeros(shape(V))
n=0
for s in x:
    solucion[:,:,n]=np.transpose(np.array([1-np.cos(4*s-tt), \
np.cos(4*s-tt)-1]),(1,0))
    n=n+1

fig, ax= subplots()
#p=ax.pcolor(X,Y,solucion[:,0,:], cmap=cm.RdBu, vmin=abs(solucion[:,0,:]).min(), vmax=abs(solucion[:,0,:]).max())
p=ax.pcolor(X,Y,solucion[:,0,:])
cb=fig.colorbar(p)



"""GRAFICO 3D"""

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig1 = plt.figure()
ax = Axes3D(fig1)
X,Y=meshgrid(x,tt)
ax.plot_surface(X, Y, V[:,0,:], rstride=1, cstride=1, cmap=cm.coolwarm,\
linewidth=0, antialiased=False)
plot(x,x)


solucion=np.zeros(shape(V))
n=0
for s in x:
    solucion[:,:,n]=np.transpose(np.array([1-np.cos(4*s-tt), \
np.cos(4*s-tt)-1]),(1,0))
    n=n+1

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, solucion[:,0,:], rstride=1, cstride=1, cmap=cm.coolwarm,\
linewidth=0, antialiased=False)
plot(x,x)

