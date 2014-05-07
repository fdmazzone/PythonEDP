# -*- coding: utf-8 -*-
"""
Created on Thu Apr 03 08:31:07 2014

@author: Noelia
"""

from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def esquemalaxwendroff(dx,dy,tmax,dt,RoValorInic,uValorInic,vValorInic,p):
    """
    RoValorInic,uValorInic,vValorInic,p son funciones de dos variables en los reales
    """
    x=linspace(0,1,int(1/dx)+1)
    y=linspace(0,1,int(1/dy)+1)
    
    ro[:,:,0]=RoValorInic(x,y)
    u[:,:,0]=uValorInic(x,y)
    v[:,:,0]=vValorInic(x,y)
    
    for t in range(0,int(tmax/dt)):
       ki:,t]=-(u[:,:,t]*d2ud2x[:,:,t]+dudx[:,:,t]**2+v[:,:,t]*d2udxdy[:,:,t]+dudy[:,:,t]*dvdx[:,:,t]+d2pd2x[:,:,t]/ro[:,:,t]-drodx[:,:,t]*dpdx[:,:,t]/(ro[:,:,t]**2))
        
        d2rodxdt[1:int(1/dx),:,t]=(drodt[2:int(1/dx)+1,:,t]-drodt[0:int(1/dx)-1,:,t])/(2*dx)
        d2rodydt[:,1:int(1/dy),t]=(drodt[:,2:int(1/dy)+1,t]-drodt[:,0:int(1/dy)-1,t])/(2*dy)
        
        d2vdxdy[:,1:int(1/dy),t]=(dvdx[:,2:int(1/dy)+1,t]-dvdx[:,0:int(1/dy)-1,t])/(2*dy)
        d2vd2y[:,1:int(1/dy),t]=(dvdy[:,2:int(1/dy)+1,t]-dvdy[:,0:int(1/dy)-1,t])/(2*dy)
        dpdy[:,1:int(1/dy),t]=(p[:,2:int(1/dy)+1,t]-p[:,0:int(1/dy)-1,t])/(2*dy)        
        d2pd2y[:,1:int(1/dy),t]=(dpdy[:,2:int(1/dy)+1,t]-dpdy[:,0:int(1/dy)-1,t])/(2*dy)        
        d2vdydt[:,:,t]=-(v[:,:,t]*d2vd2y[:,:,t]+dvdy[:,:,t]**2+u[:,:,t]*d2vdxdy[:,:,t]+dudy[:,:,t]*dvdx[:,:,t]+d2pd2y[:,:,t]/ro[:,:,t]-drody[:,:,t]*dpdy[:,:,t]/(ro[:,:,t]**2)
        
        dvdt[:,:,t]=-(u[:,:,t]*dvdx[:,:,t]+v[:,:,t]*dvdy[:,:,t]+dpdy[:,:,t]/ro[:,:,t])
        
        d2rod2t[:,:,t]=-(ro[:,:,t]*d2udxdt[:,:,t]+dudx[:,:,t]*drodt[:,:,t]+u[:,:,t]*d2rodxdt[:,:,t]+drodx[:,:,t]*dudt[:,:,t]+ro[:,:,t]*d2vdydt[:,:,t]+dvdy[:,:,t]*drodt[:,:,t]+v[:,:,t]*d2rodydt[:,:,t]+drody[:,:,t]*dvdt[:,:,t])
        
        ro[:,:,t+1]=ro[:,:,t]+drodt[:,:,t]*dt+d2rod2t[:,:,t]*(dt**2)/2
        
        dudt[:,:,t]=-(u[:,:,t]*dudx[:,:,t]+v[:,:,t]*dudy[:,:,t]+dpdx[:,:,t]/ro[:,:,t])
        
        d2ud2y[:,1:int(1/dy),t]=(dudy[:,2:int(1/dy)+1,t]-dudy[:,0:int(1/dy)-1,t])/(2*dy)
        d2pdxdy[:,1:int(1/dy),t]=(dpdx[:,2:int(1/dy)+1,t]-dpdx[:,0:int(1/dy)-1,t])/(2*dy)
        d2pdxdt[:,:,t]=(dpdx[:,:,t+1]-dpdx[:,:,t])/dt        
        d2udydt[:,:,t]=-(v[:,:,t]*d2ud2y[:,:,t]+dvdy[:,:,t]dudy[:,:,t]+u[:,:,t]*d2udxdy[:,:,t]+dudy[:,:,t]*dudx[:,:,t]+d2pdxdy[:,:,t]/ro[:,:,t]-drody[:,:,t]*dpdx[:,:,t]/(ro[:,:,t]**2)
        
        d2ud2t[:,:,t]=-(dudt[:,:,t]*dudx[:,:,t]+u[:,:,t]*d2udxdt[:,:,t]+dvdt[:,:,t]*dudy[:,:,t]+v[:,:,t]*d2udydt[:,:,t]-drodt[:,:,t]*dpdx[:,:,t]/(ro[:,:,t]**2)+d2pdxdt[:,:,t]/ro[:,:,t])
        
        u[:,:,t+1]=u[:,:,t]+dudt[:,:,t]*dt+d2ud2t[:,:,t]*(dt**2)/2
        
        d2vdxdt[:,:,t]=-(dudx[:,:,t]*dvdx[:,:,t]+u[:,:,t]*d2vd2x[:,:,t]+dvdx[:,:,t]*dvdy[:,:,t]+v[:,:,t]*d2vdxdy[:,:,t]-dpdy[:,:,t]*drodx[:,:,t]/(ro[:,:,t]**2)+d2pdxdy[:,:,t]/ro[:,:,t])
        d2pdydt[:,:,t]=(dpdy[:,:,t+1]-dpdy[:,:,t])/dt
        
        d2vd2t[:,:,t]=-(dudt[:,:,t]*dvdx[:,:,t]+u[:,:,t]*d2vdxdt[:,:,t]+dvdt[:,:,t]*dvdy[:,:,t]+v[:,:,t]*d2vdydt[:,:,t]-drodt[:,:,t]*dpdy[:,:,t]/(ro[:,:,t]**2)+d2pdydt[:,:,t]/ro[:,:,t])
        
        v[:,:,t+1]=v[:,:,t]+dvdt[:,:,t]*dt+d2vd2t[:,:,t]*(dt**2)/2
        
        return ro, u, v