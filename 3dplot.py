from mpl_toolkits.mplot3d.axes3d import Axes3D
from pylab import *
import matplotlib.pyplot as plt
import numpy
	
#param = alpha, beta, gamma, delta1, delta2
param=[1.0,1.0,1.0,1.0,1.0]
psi = linspace(-1,1,100)
phi = linspace(-1,1,100)
v1=(psi**2-param[3])*(psi**2-1)**2
v2=(param[0]
v=v1+v2
#unos= numpy.ones(100)
X,Y = meshgrid(psi, phi)
un=X/X
Z=V1+V2

fig = plt.figure(figsize=(14,6))

# `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
ax = fig.add_subplot(1, 2, 1, projection='3d')

p = ax.plot_surface(X, Y, Z, rstride=4, cstride=4, linewidth=0)

# surface_plot with color grading and color bar
ax = fig.add_subplot(1, 2, 2, projection='3d')
p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)
#plt.show()
