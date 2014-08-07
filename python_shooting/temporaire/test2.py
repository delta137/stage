#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)
#version 3: condition frontiere ajustable

#librairies
import numpy
import math
import pylab

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

#definition du probleme differentiel
def rhs(y,param):
    """ RHS function.  ici y1 = psi, y2 = psi' 
        z1 = phi, z2 = phi' 
        Alors: y2' = phi'' = lambda*phi**3 - m**2*phi """

    dydx = [y[1],
    param[0]*(y[0]**3)-(param[1]**2)*y[0]] 

    return dydx
init=[y1_0,y2_0[0]]
 
x = numpy.linspace(xl, xr, n)
h = x[1] - x[0]  # stepsize
    
#def pente(x): return (x-xl)/(xr-xl)
#y1 = map(pente, range(xl, xr))
y=numpy.zeros((q,n))
# left boundary initialization
y[0,0] = init[0]
y[1,0] = init[1]
k = 0
