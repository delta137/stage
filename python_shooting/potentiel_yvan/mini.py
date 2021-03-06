#librairies
import numpy as np
import math
import pylab
import scipy
from scipy.integrate import odeint
from scipy.optimize import broyden1,fsolve as broy,fsolve

#range d'integration, bornes et stepsize
xl, xr, n = 0,2,100
init=0
guess=0
y0 = [init,guess]
x = np.linspace(xl,xr,n)
#systeme dequations differentielles
#y est un vecteur, avec y= (y1,y2)
#y1 est la position, y2=y1'
def rhs(y,x):
    dy1dx=y[1]
    dy2dx=y[0]
    return np.array([[dy1dx],[dy2dx]])
#def jac(y,x):
#    return np.array([[0,1],[1,0]])
    
y=odeint(rhs,y0,x)
pylab.scatter(x, y[:,0])
pylab.show()


    
