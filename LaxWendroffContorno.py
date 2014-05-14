# -*- coding: utf-8 -*-
"""
Created on Wed May 14 03:49:56 2014

@author: Noelia
"""
import numpy as np

dx=.1
dt=.1
tmax=1
RoValorInic=lambda x:x
uValorInic=lambda x:x
A=lambda u, ro:np.array([[-u, ro**2], [-ro, -u]])
Au=lambda u, ro:np.array([[-np.ones(np.shape(ro)), np.zeros(np.shape(ro))], \
[np.zeros(np.shape(ro)), -np.ones(np.shape(ro))]])
Aro=lambda u, ro:np.array([[np.zeros(np.shape(ro)), 2*ro], [-np.ones(np.shape(ro)), \
np.zeros(np.shape(ro))]])






""" CONDICIONES DE CONTORNO DE DIRICHLET """
#en el borde izquierdo es una funcuón conocinda ContIzq (dada en forma vectorial)
#en el borde derecho es una funcuón conocinda ContDer (dada en forma vectorial)

nPuntos=int(1/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,1,nPuntos)
tt=linspace(dt,1,nPasos)
V=np.zeros((nPasos,2,nPuntos)) #Inicializo V
V[0,0,:]=uValorInic(x) #Cálculo primer nivel para la variable u
V[0,1,:]=RoValorInic(x) #Cálculo primer nivel para la variable ro
V[:,:,0]=ContIzq(tt)
V[:,:,-1]=ContDer(tt)
V0=V[0,:,:] #Guardo el primer nivel para construir el siguiente
AA=A(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo A en u y ro en el tiempo t=0
AAu=Au(V[0,0,1:-1],V[0,1,1:-1]) ##Evalúo Au en u y ro en el tiempo t=0
AAro=Aro(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t=0
Vx=Vxx=Vt=Vxt=Vtt=np.zeros(np.shape(V0[:,1:-1])) #Inicializo Vx, Vxx, Vt, Vxt y Vtt
AAx=AAt=np.zeros(np.shape(AA)) #Inicializo AAx y AAt
cte1=2*dx
cte2=dx**2
cte3=cte2/2
for t in range(1,nPasos):
    Vx=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vxx=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
    
    Vt=sum(AA*Vx,1) #Cálculo de Vt por medio de la ecuación diferencial donde \
    # sum(AA*Vx,1) es un producto matricial
    
    AAx=AAu*Vx[0]+AAro*Vx[1]
    AAt=AAu*Vt[0]+AAro*Vt[1]
    
    Vxt=sum(AAx*Vx,1)+sum(AA*Vxx,1) #Cálculo de Vxt por medio de derivar la \
    # ecuación diferencial respecto de x
    Vtt=sum(AAt*Vx,1)+sum(AA*Vxt,1) #Cálculo de Vtt por medio de derivar la \
    # ecuación diferencial respecto de t
    
    V[t,:,1:-1]=V0[:,1:-1]+Vt*dt+Vtt*cte3 #Cálculo de la solución por Taylor \
    # salvo en los bordes
        
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t



""" CONDICIONES DE CONTORNO DE NEUMANN (MIXTO) """

#en el borde izquierdo la derivada respecto de x es cero
#en el borde derecho es cero

"Aproximacion de primer orden"

nPuntos=int(1/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,1,nPuntos)
tt=linspace(dt,1,nPasos)
V=np.zeros((nPasos,2,nPuntos)) #Inicializo V
V[0,0,:]=uValorInic(x) #Cálculo primer nivel para la variable u
V[0,1,:]=RoValorInic(x) #Cálculo primer nivel para la variable ro
V0=V[0,:,:] #Guardo el primer nivel para construir el siguiente
AA=A(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo A en u y ro en el tiempo t=0
AAu=Au(V[0,0,1:-1],V[0,1,1:-1]) ##Evalúo Au en u y ro en el tiempo t=0
AAro=Aro(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t=0
Vx=Vxx=Vt=Vxt=Vtt=np.zeros(np.shape(V0[:,1:-1])) #Inicializo Vx, Vxx, Vt, Vxt y Vtt
AAx=AAt=np.zeros(np.shape(AA)) #Inicializo AAx y AAt
cte1=2*dx
cte2=dx**2
cte3=cte2/2
for t in range(1,nPasos):
    Vx=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vxx=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
    
    Vt=sum(AA*Vx,1) #Cálculo de Vt por medio de la ecuación diferencial donde \
    # sum(AA*Vx,1) es un producto matricial
    
    AAx=AAu*Vx[0]+AAro*Vx[1]
    AAt=AAu*Vt[0]+AAro*Vt[1]
    
    Vxt=sum(AAx*Vx,1)+sum(AA*Vxx,1) #Cálculo de Vxt por medio de derivar la \
    # ecuación diferencial respecto de x
    Vtt=sum(AAt*Vx,1)+sum(AA*Vxt,1) #Cálculo de Vtt por medio de derivar la \
    # ecuación diferencial respecto de t
    
    V[t,:,1:-1]=V0[:,1:-1]+Vt*dt+Vtt*cte3 #Cálculo de la solución por Taylor \
    # salvo en los bordes
    V[t,:,0]=V[t,:,1] #Completo la solución en el borde izquierdo
        
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t

"Aproximacion de segundo orden (pto fantasma)"

nPuntos=int(1/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,1,nPuntos)
tt=linspace(dt,1,nPasos)
V=np.zeros((nPasos,2,nPuntos)) #Inicializo V
V[0,0,:]=uValorInic(x) #Cálculo primer nivel para la variable u
V[0,1,:]=RoValorInic(x) #Cálculo primer nivel para la variable ro
V0=V[0,:,:] #Guardo el primer nivel para construir el siguiente
AA=A(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo A en u y ro en el tiempo t=0
AAu=Au(V[0,0,1:-1],V[0,1,1:-1]) ##Evalúo Au en u y ro en el tiempo t=0
AAro=Aro(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t=0
Vx=Vxx=Vt=Vxt=Vtt=np.zeros(np.shape(V0[:,1:-1])) #Inicializo Vx, Vxx, Vt, Vxt y Vtt
AAx=AAt=np.zeros(np.shape(AA)) #Inicializo AAx y AAt
cte1=2*dx
cte2=dx**2
cte3=cte2/2
cte4=2/cte2
for t in range(1,nPasos):
    Vx=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vxx=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
    
    Vt=sum(AA*Vx,1) #Cálculo de Vt por medio de la ecuación diferencial donde \
    # sum(AA*Vx,1) es un producto matricial
    
    AAx=AAu*Vx[0]+AAro*Vx[1]
    AAt=AAu*Vt[0]+AAro*Vt[1]
    
    Vxt=sum(AAx*Vx,1)+sum(AA*Vxx,1) #Cálculo de Vxt por medio de derivar la \
    # ecuación diferencial respecto de x
    Vtt=sum(AAt*Vx,1)+sum(AA*Vxt,1) #Cálculo de Vtt por medio de derivar la \
    # ecuación diferencial respecto de t
    
    V[t,:,1:-1]=V0[:,1:-1]+Vt*dt+Vtt*cte3 #Cálculo de la solución por Taylor \
    # salvo en los bordes
    V[t,:,0]=V0[:,0]+sum(AA(V[t,0,0],V[t,1,0]),sum(AA(V[t,0,0],V[t,1,0]),\
    cte4*(Vx[:,1]-Vx[:,0]),1),1)#Completo la solución en el borde izquierdo
        
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t



#en el borde izquierdo la derivada respecto de x es una funcuón conocinda \ContIzq (dada en forma vectorial)
#en el borde derecho es una funcuón conocinda ContDer (dada en forma vectorial)

"Aproximacion de primer orden"

nPuntos=int(1/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,1,nPuntos)
tt=linspace(dt,1,nPasos)
V=np.zeros((nPasos,2,nPuntos)) #Inicializo V
V[0,0,:]=uValorInic(x) #Cálculo primer nivel para la variable u
V[0,1,:]=RoValorInic(x) #Cálculo primer nivel para la variable ro
V[:,:,-1]=ContDer(tt)
V0=V[0,:,:] #Guardo el primer nivel para construir el siguiente
AA=A(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo A en u y ro en el tiempo t=0
AAu=Au(V[0,0,1:-1],V[0,1,1:-1]) ##Evalúo Au en u y ro en el tiempo t=0
AAro=Aro(V[0,0,1:-1],V[0,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t=0
Vx=Vxx=Vt=Vxt=Vtt=np.zeros(np.shape(V0[:,1:-1])) #Inicializo Vx, Vxx, Vt, Vxt y Vtt
AAx=AAt=np.zeros(np.shape(AA)) #Inicializo AAx y AAt
cte1=2*dx
cte2=dx**2
cte3=cte2/2
for t in range(1,nPasos):
    Vx=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vxx=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
    
    Vt=sum(AA*Vx,1) #Cálculo de Vt por medio de la ecuación diferencial donde \
    # sum(AA*Vx,1) es un producto matricial
    
    AAx=AAu*Vx[0]+AAro*Vx[1]
    AAt=AAu*Vt[0]+AAro*Vt[1]
    
    Vxt=sum(AAx*Vx,1)+sum(AA*Vxx,1) #Cálculo de Vxt por medio de derivar la \
    # ecuación diferencial respecto de x
    Vtt=sum(AAt*Vx,1)+sum(AA*Vxt,1) #Cálculo de Vtt por medio de derivar la \
    # ecuación diferencial respecto de t
    
    V[t,:,1:-1]=V0[:,1:-1]+Vt*dt+Vtt*cte3 #Cálculo de la solución por Taylor \
    # salvo en los bordes
    V[t,:,0]=V[t,:,1]-ValorDeriv(tt[t])*dx #Completo la solución en el borde izquierdo
        
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t


"Aproximación de segundo orden"

nPuntos=int(1/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,1,nPuntos)
tt=linspace(dt,1,nPasos)
V=np.zeros((nPasos,2,nPuntos)) #Inicializo V
V[0,0,:]=uValorInic(x) #Cálculo primer nivel para la variable u
V[0,1,:]=RoValorInic(x) #Cálculo primer nivel para la variable ro
V[:,:,-1]=ContDer(tt)
V0=V[0,:,:] #Guardo el primer nivel para construir el siguiente
AA=A(V[0,0,:-1],V[0,1,:-1]) #Evalúo A en u y ro en el tiempo t=0
AAu=Au(V[0,0,:-1],V[0,1,:-1]) ##Evalúo Au en u y ro en el tiempo t=0
AAro=Aro(V[0,0,:-1],V[0,1,:-1]) #Evalúo Aro en u y ro en el tiempo t=0
Vx=Vxx=Vt=Vxt=Vtt=np.zeros(np.shape(V0[:,:-1])) #Inicializo Vx, Vxx, Vt, Vxt y Vtt
AAx=AAt=np.zeros(np.shape(AA)) #Inicializo AAx y AAt
cte1=2*dx
cte2=dx**2
cte3=cte2/2
cte4=2/cte2
for t in range(1,nPasos):
    Vx[:,-1:]=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vx[:,0]=ValorDeriv(tt[t])
    Vxx[:,-1:]=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
    Vxx[:,0]=cte4*(V0[:,1]-V0[:,0])-Vx[:,0]
    
    Vt=sum(AA*Vx,1) #Cálculo de Vt por medio de la ecuación diferencial donde \
    # sum(AA*Vx,1) es un producto matricial
    
    AAx=AAu*Vx[0]+AAro*Vx[1]
    AAt=AAu*Vt[0]+AAro*Vt[1]
    
    Vxt=sum(AAx*Vx,1)+sum(AA*Vxx,1) #Cálculo de Vxt por medio de derivar la \
    # ecuación diferencial respecto de x
    Vtt=sum(AAt*Vx,1)+sum(AA*Vxt,1) #Cálculo de Vtt por medio de derivar la \
    # ecuación diferencial respecto de t
    
    V[t,:,:-1]=V0[:,:-1]+Vt*dt+Vtt*cte3 #Cálculo de la solución por Taylor \
    
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,:-1],V[t,1,:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,:-1],V[t,1,:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,:-1],V[t,1,:-1]) #Evalúo Aro en u y ro en el tiempo t


