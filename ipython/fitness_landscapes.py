# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 10:22:39 2016

@author: colinsh
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(20,20))
ax = fig.gca(projection='3d')
r=6
X = np.arange(-r*np.pi, r*np.pi, .5*np.pi)
Y = np.arange(-r*np.pi, r*np.pi, .5*np.pi)
X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
R = np.sqrt(X**2 + Y**2)
coff=0
soff=.5

#Peaked
#Z = -(.05*R)**2+.5
#for i in xrange(len(Z)):
#    for j in xrange(len(Z[i])):
#        if Z[i][j]<0:
#            Z[i][j]*=0


#print Z
#Time 1
#Z = (np.sin(.25*X+soff*np.pi)*(1.3*np.cos(.25*Y+coff*np.pi)))

#Time 2
Z = abs((np.sin(.25*X+soff*np.pi))*(1.3*np.cos(.25*Y+coff*np.pi)))#/(.2*(R+2))#np.exp(-.1*R))#/(R+1))

#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

for i in xrange(len(Z)):
    for j in xrange(len(Z[i])):
        if Z[i][j]>3:
            Z[i][j]=1
        if Z[i][j]<0:
            Z[i][j]*=-.5*np.random.rand(1,1)
        if Z[i][j]<0.1:
            Z[i][j]+=.7*np.random.rand(1,1)
        Z[i][j]*=.5*np.random.rand(1,1)+.1*np.random.rand(1,1)   
#        
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='lightgray',
                       linewidth=.5, antialiased=False)
ax.set_zlim(0, 3)
plt.axis('off')
ax.set_axis_bgcolor('w')

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()