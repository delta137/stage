import numpy as np
import math
import pylab
import scipy
from scipy.integrate import quad

graph=True

# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

xl=0.001
xr=8
n=100
x = np.linspace(xl, xr, n)

eps=0.1
lam=5
mu=2
def rhs(y,x):
    """ RHS function.  ici y[0,:] = psi, y[1,:] = psi' 
        y[2,:]= phi, y[3,:] = phi' """
    d1=y[1]
    d2=(1/(4*x))*(math.sqrt(5)*eps*x-2*mu**2*y[0]*x+2*lam*x*y[0]**3-12*y[1])
    dydx = np.asfarray([d1,d2], dtype='float') 
    return dydx
yi=-0.910
yip=0
y1=scipy.integrate.odeint(rhs,[yi,yip],x)

pylab.scatter(x, y1[:,0])
pylab.show()
