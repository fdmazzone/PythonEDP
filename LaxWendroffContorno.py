# -*- coding: utf-8 -*-
"""
Created on Wed May 14 03:49:56 2014

@author: Noelia
"""
import numpy as np

dx=.01
dt=.01
tmax=2 #tamaño del intervalo para t
intervalo=np.pi/2 #tamaño del intervalo para x
RoValorInic=lambda x:np.cos(4*x)-1 #condición inicial para ro en t=0
uValorInic=lambda x:1-np.cos(4*x) #condición inicial para u en t=0
#matriz del sistema Vt=A(V)Vx
A=lambda u, ro:np.array([[np.zeros(np.shape(u)), np.ones(np.shape(ro))], \
[np.ones(np.shape(u)), np.zeros(np.shape(ro))]])
#derivada de la matriz del sistema respecto de u
Au=lambda u, ro:np.array([[np.zeros(np.shape(ro)), np.zeros(np.shape(ro))], \
[np.zeros(np.shape(ro)), np.zeros(np.shape(ro))]])
#derivada de la matriz del sistema respecto de ro
Aro=lambda u, ro:np.array([[np.zeros(np.shape(ro)), np.zeros(np.shape(ro))], \
[np.zeros(np.shape(ro)), np.zeros(np.shape(ro))]])



""" CONDICIONES DE CONTORNO DE DIRICHLET """

#en el borde izquierdo es una funcuón conocinda ContIzq (dada en forma vectorial)
#en el borde derecho es una funcuón conocinda ContDer (dada en forma vectorial)

ContIzq=lambda tt: np.transpose(np.array([np.zeros(np.shape(tt)), \
np.zeros(np.shape(tt))]),(1,0))

ContDer=lambda tt: np.transpose(np.array([1-np.cos(2*np.pi+tt), \
np.cos(2*np.pi-tt)-1]),(1,0))

nPuntos=int(intervalo/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,intervalo,nPuntos)
tt=np.linspace(0,tmax,nPasos)
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

"Primer caso"
#en el borde izquierdo la derivada respecto de x es cero
#en el borde derecho es cero

"Aproximacion de primer orden"

nPuntos=int(intervalo/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,intervalo,nPuntos)
tt=np.linspace(0,tmax,nPasos)
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

nPuntos=int(intervalo/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,intervalo,nPuntos)
tt=np.linspace(0,tmax,nPasos)
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
    V[t,:,0]=V0[:,0]+sum(A(V[t,0,0],V[t,1,0])*sum(A(V[t,0,0],V[t,1,0])*(cte4*(Vx[:,1]-Vx[:,0])),1),1)#Completo la solución en el borde izquierdo
        
    V0= V[t,:,:]  #Guardo el nivel contruido para construir el siguiente
    AA=A(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo A en u y ro en el tiempo t
    AAu=Au(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Au en u y ro en el tiempo t
    AAro=Aro(V[t,0,1:-1],V[t,1,1:-1]) #Evalúo Aro en u y ro en el tiempo t



"Segundo Caso"

#en el borde izquierdo la derivada respecto de x es una funcuón conocinda \
#ContIzq (dada en forma vectorial)
#en el borde derecho es una funcuón conocinda ContDer (dada en forma vectorial)

ValorDeriv=lambda tt: np.array([4*np.sin(-tt), -4*np.sin(-tt)])

ContDer=lambda tt: np.transpose(np.array([1-np.cos(2*np.pi+tt), \
np.cos(2*np.pi-tt)-1]),(1,0))

"Aproximacion de primer orden"

nPuntos=int(intervalo/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,intervalo,nPuntos)
tt=np.linspace(0,tmax,nPasos)
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

"Aproximación de segundo orden (pto fantasma)"

nPuntos=int(intervalo/dx)+1
nPasos=int(tmax/dt)+1
x=np.linspace(0,intervalo,nPuntos)
tt=np.linspace(0,tmax,nPasos)
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
    Vx[:,1:]=(V0[:,2:]-V0[:,0:-2])/cte1 #Cálculo de Vx por diferencia central
    Vx[:,0]=ValorDeriv(tt[t])
    Vxx[:,1:]=(V0[:,2:]-2*V0[:,1:-1]+V0[:,0:-2])/cte2 #Cálculo de Vxx por diferencia central
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


