#!/home/delta137/anaconda/bin/python
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
xl=0
xr=0.1
Pdrop = -100
mu = 1

def odefun(U, y):
    u1, u2 = U
    du1dy = u2
    du2dy = 1.0 / mu * Pdrop
    return [du1dy, du2dy]

u1_0 = 0 # known
xspan = np.linspace(xl, xr)

def objective(u2_0):
    xspan = np.linspace(xl, xr)
    U = odeint(odefun, [u1_0, u2_0], xspan)
    u1 = U[:,0]
    return u1[-1]

u2_0, = fsolve(objective, 1.0)

# now solve with optimal u2_0
U = odeint(odefun, [u1_0, u2_0], xspan)

plt.plot(xspan, U[:,0], label='Numerical solution')
plt.plot([xr],[0], 'ro')

# plot an analytical solution
u = -(Pdrop) * xr**2 / 2 / mu * (xspan / (xr-xl) - (xspan / (xr-xl))**2)
plt.plot(xspan, u, 'r--', label='Analytical solution')
plt.xlabel('x')
plt.ylabel('$u_1$')
plt.legend(loc='best')
plt.savefig('shoot3.png')
plt.show()
