import numpy as np
import math
import pylab
import scipy
from scipy.integrate import ode
from scipy.optimize import fsolve,root

def rhs(y):
    dydx = np.array([y[1],
                    -y[0]]) 
    return dydx

xl, xr, n = 0,2*math.pi,100
q=2

init=np.array(0.0) #valeur connue, y(0)
cible=np.array(0) #y(2pi) est aussi connue, c'est aussi 0
guess=np.array(0.95) #on prend un guess pour y'(0)
V=0

x = np.linspace(xl, xr, n)
h = x[1] - x[0]  # stepsize

#def pente(x): return (x-xl)/(xr-xl)
#y1 = map(pente, range(xl, xr))
y=np.zeros((q,n))
    # left boundary initialization
y[0,0] = init
y[1,0] = V
k = 0
while (k < n-1):
    
    dydx1= rhs(y[:,k])
    dydx2= rhs(y[:,k] + 0.5*h*dydx1)
    dydx3= rhs(y[:,k] + 0.5*h*dydx2)
    dydx4= rhs(y[:,k] + h*dydx3)

    y[:,k+1] = y[:,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)                        
    k += 1
    
x = np.linspace(xl, xr, n)

pylab.scatter(x, y[0,:])
pylab.show()

