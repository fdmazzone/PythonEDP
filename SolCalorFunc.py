from numpy import *
def SolCalor(m,v,a,b,f,tmax):
    dx=1.0/m
    dt=.5*dx**2
    factor=dt/(dx**2)
    M=int(tmax/dt)
    u=zeros([m+1,M+1])
    u[:,0]=f(linspace(0,1,m+1,endpoint=True))
    u[0,:]=a(linspace(0,M*dt,M+1,endpoint=True))
    u[m,:]=b(linspace(0,M*dt,M+1,endpoint=True))
    for n in range(M):
        u[2:m,n+1]=u[2:m,n]+v*(u[3:m+1,n]-2.0*u[2:m,n]+u[1:m-1,n])*factor
    X,Y = meshgrid(linspace(0,1,m+1,endpoint=True), linspace(0,M*dt,M+1,endpoint=True))
    Z=u.T
    return X,Y,Z
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#import matplotlib.pyplot as plt
#from pylab import *
#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)
#plt.show()
#savefig('test')