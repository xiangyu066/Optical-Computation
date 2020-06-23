import random
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
r_s=0.0                   
r_e=1.75                     
r_step=1000                
iters=100
nSeeds=100

###############################################################################
r=np.linspace(r_s,r_e,num=r_step)
data_r=np.zeros((len(r),nSeeds))
data_x=np.zeros((len(r),nSeeds))

###############################################################################
for rCount in range(0,len(r)-1):
    for nSeed in range(0,nSeeds-1):
        x_n=random.random()
        for iter in range(0,iters-1):
            x_n=r[rCount]*x_n+(1-x_n)**3
        data_r[rCount,nSeed]=r[rCount]
        data_x[rCount,nSeed]=x_n
plt.plot(data_r,data_x,'.',color='black')
plt.show()
