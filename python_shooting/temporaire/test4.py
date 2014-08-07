#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)
#version 3: condition frontiere ajustable
#Yvan: j'ai implemente des arrays, plus faciles avec plus d'equations

#librairies
import numpy
import math
import pylab
import scipy
from scipy.integrate import ode

# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

#************Configurations********************#
graph=True #montrer le graph a la fin ou non

# param = alpha, beta, gamma, delta1, delta2 
param=[1,1]

q=2 #nombre de degres de liberts

#histogramme
xl, xr, n = 0,4,100
npts=n

#condition initiales
y1_0 = 0.0        #connu
y2_0 = [0.3, 0.9] #plage de valeurs qu'on va etudier

#cible du shooting (y1(xr))
cible = param[1]/param[0]**(0.5) *(1-math.exp(-2**(0.5)*xr)) 

#tolerance
err_min = 1e-5
#********************************#

#************Fonctions********************#
#integration par methode de Runge-Kutta 4
"""
def rk4(init, param, rhs, xl, xr, n,q):
           #R-K 4 integration:
           #y1_0 and y2_0 are y1(0) and y2(0)
           #rhs is the righthand side function
           #xl and xr are the domain limits
           #n is the number of integration points 

    x = numpy.linspace(xl, xr, n)
    h = x[1] - x[0]  # stepsize
    
    #def pente(x): return (x-xl)/(xr-xl)
    #y1 = map(pente, range(xl, xr))
    y=numpy.zeros((q,n))
        # left boundary initialization
    y[0,0] = init[0]
    y[1,0] = init[1]
    k = 0
    while (k < n-1):
        
        dydx1 = rhs(y[:,k],param)
        dydx2= rhs(y[:,k] + 0.5*h*dydx1,param)
        dydx3= rhs(y[:,k] + 0.5*h*dydx2,param)
        dydx4= rhs(y[:,k] + h*dydx3,param)

        y[:,k+1] = y[:,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)                        
        k += 1

    return y
    """

#definition du probleme differentiel
def f(t,y,arg1):
    """ RHS function.  ici y1 = psi, y2 = psi' 
        z1 = phi, z2 = phi' 
        Alors: y2' = phi'' = lambda*phi**3 - m**2*phi """

    dydx = numpy.array([y[1],
    param[0]*(y[0]**3)-(param[1]**2)*y[0]]) 

    return dydx

#soln analytique (kink)
def analytic(x,param):
    """ analytic solution """
    return (param[1]/param[0]**(0.5))*numpy.tanh((param[1]/2**(0.5)) * x )

#jacobien    
def jac(t, y,arg1):
    """ J_{i,j} = df_i/dy_j """

    y1 = y[0]
    y2 = y[1]
    
    
    df1dy1 = 0
    df1dy2 = 1 

    df2dy1 = 3*param[0]*y1**2-param[1]**2
    df2dy2 = 0
    

    return numpy.array([ [ df1dy1, df1dy2],
                         [ df2dy1, df2dy2]])
                         
x = numpy.linspace(xl,xr,n)
iter =1
itermax=10

    #1ere integration avec la borne inf de y2_0
init=[y1_0,y2_0[0]]

  
t0=0
r = ode(f, jac).set_integrator('vode', method='bdf', with_jacobian=True)
r.set_initial_value(init, t0).set_f_params(2.0).set_jac_params(2.0)
t1 = xr
dt = 1
xout = [r.t]
y1out = [init[0]]
y2out = [init[1]]
while r.successful() and r.t < t1:
    r.integrate(r.t+dt)
    xout.append(r.t)
    y1out.append(r.y[0])
    y2out.append(r.y[1])
    

"""
r = ode(f, jac).set_integrator("DVODE", method="bdf", 
                                     with_jacobian=True,
                                     atol=1.e-10, rtol=1.e-10,
                                     nsteps = 15000, order=5) #, min_step=dx)

x = 0.0
r.set_initial_value(init, x)

dt=xr/n
tout = [x]
y1out = [init[0]]
y2out = [init[1]]
#y3out = [init[2]]
while r.successful() and r.t < xr:
    r.integrate(r.t+dt)
    dt=10*dt

        xout.append(r.t)
        y1out.append(r.y[0])
        y2out.append(r.y[1])

        dt = dt

    return numpy.array(xout), \
        numpy.array([[y1out],[y2out]])     
"""
