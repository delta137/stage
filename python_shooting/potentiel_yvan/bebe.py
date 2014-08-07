import scipy
from scipy.integrate import ode
from scipy.optimize import fsolve,root
import numpy as np
def F(x):
    return lambda x: x**2+2*x+1
sol=root(F,np.array([0.5,0.6]))
