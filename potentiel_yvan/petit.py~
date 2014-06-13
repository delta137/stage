import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def fun(y, t, a):
    """Define the right-hand side of equation dy/dt = a*y""" 
    f = a * y
    return f


# Initial condition
y0 = 100.0

# Times at which the solution is to be computed.
t = np.linspace(0, 1, 51)

# Parameter value to use in `fun`.
a = -2.5

# Solve the equation.
y = odeint(fun, y0, t, args=(a,))

# Plot the solution.  `odeint` is generally used to solve a system
# of equations, so it returns an array with shape (len(t), len(y0)).
# In this case, len(y0) is 1, so y[:,0] gives us the solution.
plt.plot(t, y[:,0])
plt.xlabel('t')
plt.ylabel('y')
plt.show()
