#librairies
import numpy as np
import math
import pylab
import scipy
from scipy.integrate import ode
from scipy.optimize import broyden1 as broy

colors = ["k", "r", "g", "b", "c"]
#************Configurations********************#
graph=True          #montrer le graph a la fin ou non
p=[1,1,1,1,1]       # param = alpha, beta, gamma, delta1, delta2 
q=4                 #nombre de degres de libertes
xl, xr, n = 0,10,100#histogramme
x = np.linspace(xl,xr,n)
err_min = 1e-5      #tolerance

#soit y le vecteur y=[psi,psi',phi,phi']
#condition frontieres
#initiales
y1_0,y4_0 = 0.0,0.0     #connues
V=np.array([1,1])   #devinettes

#a l'infini, eq diff qui decrit deux conditions
def cible(y):
    sig1=(8*(1-p[3]))**0.5
    sig2=((p[0]*(16+3*p[4]))/(2*p[1]*(1+p[2])))**0.5 
    cond1=y[1]+sig1*(y[0]-1)
    cond2=y[3]+sig2*(y[2]+1)
    return np.array([[cond1],[cond2]])

#********************************************#

#************Fonctions********************#
#integration par methode de Runge-Kutta 4

#definition du probleme differentiel
#j'ai confirme ce quil y avait sur la feuille de Richard
def rhs(y,x):
    a1=(y[0]**2-1)
    a2=(y[0]**2-p[3])
    a3=(y[2]**2-1)
    a4=(p[2]+y[0]**2)
    a5=y[2]+1
    a6=y[2]-2
    a7=(4*y[2]-3*p[4]/4)
    
    dy1dx=y[1]
    dy2dx=2*y[0]*(2*a1*a2+a1**2-p[0]*(a3**2-(p[4]/4)*a5*a6)/a4)
    dy3dx=y[3]
    dy4dx=p[1]*a3*a7/a4
    
    return np.array([[dy1dx],[dy2dx],[dy3dx],[dy4dx]])
init=[0,0]         
h = x[1] - x[0]            # stepsize
y=np.zeros((q,n))          #y=[psi,psi',phi,phi']        
y[0,0], y[3,0] = y1_0,y4_0 #conditions intiales
y[1,0], y[2,0] = init[0],init[1] #guess   
k = 0

while (k < n-1):
    dydx1 = rhs(y[:,k],k)
    dydx2= rhs(y[:,k] + 0.5*h*dydx1,k)
    dydx3= rhs(y[:,k] + 0.5*h*dydx2,k)
    dydx4= rhs(y[:,k] + h*dydx3,k)
    #y[:,k+1] = y[:,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)                        
    k += 1

