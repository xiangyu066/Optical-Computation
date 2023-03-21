# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 21:17:56 2022

@author: XYZ
"""

#%%
import math
import numpy as np
import matplotlib.pyplot as plt
print('Running...')

#%% Define units
mm          =1

#%% Configuration
n_air       =1.0                # the refractive index of the air
n_lens      =2.5                # the refractive index of a lens
R1          =1000*(mm)          # the first interface
R2          =-3.6*(mm)          # the second interface
thickness   =5*(mm)             # the thickness of a lens
beamsize    =6*(mm)             # the diameter of bean
nLights     =6                  # the number of simulated lights

#%% Define ray matrix
def Refraction(n_t,n_i,R):
    R_matrix=np.array([[1,-(n_t-n_i)/R],[0,1]])
    return R_matrix
    
def Transfer(d,n_t):
    T_matrix=np.array([[1,0],[d/n_t,1]])
    return T_matrix

#%% the light propagate through a plano-convex lens
x1=np.zeros((nLights,6))
y1=np.zeros((nLights,6))

for nLight in range(nLights):
    # initial position of a source
    r=np.array([[n_air*0],[beamsize/2-nLight*beamsize/(nLights-1)]])
    l=np.array([-2*mm])
    x1[nLight,0]=l
    y1[nLight,0]=r[1]
    
    r=np.matmul(Transfer(2*mm, n_air),r)
    x1[nLight,1]=x1[nLight,0]+2*mm
    y1[nLight,1]=r[1]
    
    r=np.matmul(Refraction(n_lens, n_air,R1),r)
    x1[nLight,2]=x1[nLight,1]
    y1[nLight,2]=r[1]
    
    th_=np.arccos(r[1]/abs(R2))
    x_=R2+thickness+abs(R2)*np.sin(th_)
    r=np.matmul(Transfer(x_, n_lens),r)
    x1[nLight,3]=x1[nLight,2]+x_
    y1[nLight,3]=r[1]
    
    r=np.matmul(Refraction(n_air, n_lens,R2),r)
    x1[nLight,4]=x1[nLight,3]
    y1[nLight,4]=r[1]
    
    r=np.matmul(Transfer(13*mm, n_lens),r)
    x1[nLight,5]=x1[nLight,4]+13*mm
    y1[nLight,5]=r[1]

plt.figure()
th=np.linspace(0.5*math.pi,1.5*math.pi,1000)
plt.plot(R1+R1*np.cos(th),0+R1*np.sin(th),'k',lw=3)
th=np.linspace(-0.5*math.pi,0.5*math.pi,1000)
plt.plot(R2+thickness+abs(R2)*np.cos(th),0+abs(R2)*np.sin(th),'k',lw=3)
plt.plot(x1.T,y1.T,lw=1)
plt.grid(lw=0.3)
plt.xlabel('X [mm]',weight='bold')
plt.ylabel('Y [mm]', weight='bold')
plt.xlim((-2,12))
plt.ylim((-3.5,3.5))
plt.show()

plt.figure()
plt.plot(x1.T,y1.T,lw=1)
plt.grid(lw=0.3)
plt.xlabel('X [mm]',weight='bold')
plt.ylabel('Y [mm]', weight='bold')
plt.xlim((8,12))
plt.ylim((-1,1))
plt.show()

#%%
print('Done.')


