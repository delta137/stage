#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)
#version 3: condition frontiere ajustable
#Yvan: j'ai implemente des arrays, plus faciles avec plus d'equations
#Yvan3: apres multiples essais avec Jacobien et fonction 
#je reviens a mon code qui marchait bien lol

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
xl, xr, n = 0,2,100
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
    
#definition du probleme differentiel
def rhs(y,param):
    """ RHS function.  ici y[0,:] = psi, y[1,:] = psi' 
        y[2,:]= phi, y[3,:] = phi'"""

    dydx = numpy.array([y[1],
    param[0]*(y[0]**3)-(param[1]**2)*y[0]]) 

    return dydx

#soln analytique (kink)
def analytic(x,param):
    """ analytic solution """
    return (param[1]/param[0]**(0.5))*numpy.tanh((param[1]/2**(0.5)) * x )


                
#******************************************#
    
    
#//////////////////////PROBLEME/////////////////////////



#on va shooter dans un range de conditions intiales,
#la conditon pour la vitesse initiale est inconnue, mais 
#une reponse optimale est trouvee en comparant avec la 
#deuxieme condition frontiere (pas initiale)

#////#boucle d'integration, pour la victoire!/////
x = numpy.linspace(xl,xr,n)
iter =1
itermax=15
while (iter < itermax):
    #1ere integration avec la borne inf de y2_0
    init=[y1_0,y2_0[0]]
    x,y =rk4(init, param, rhs, xl, xr, n,q)
    err1= y[0,-1]-cible

