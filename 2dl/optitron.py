#librairies
import numpy as np
import math
import pylab
import scipy
from scipy.integrate import ode
from scipy.optimize import root
from scipy.optimize import broyden1 as broy
#from scipy.optimize import fsolve 

colors = ["k", "r", "g", "b", "c"]
#************Configurations********************#
graph=True          #montrer le graph a la fin ou non
p=[1,1,1,1,1]       # param = alpha, beta, gamma, delta1, delta2 
q=4+5            #nombre de degres de libertes
xl, xr, n = 0,2,100#histogramme

#soit y le vecteur y=[psi,psi',phi,phi']
#condition frontieres
#initiales
y1_0,y4_0 = 0.0,0.0     #connues
y2_0,y3_0=0.25,0.5
V=np.array([y2_0,y3_0,p[0],p[1],p[2],p[3],p[4]])   #devinettes
origine=V

#a l'infini, eq diff qui decrit deux conditions
def cible(y):
    p=y[4:y.size] # les parametres a la fin du vecteur
    sig1=(8*(1-p[3]))**0.5
    sig2=((p[0]*(16+3*p[4]))/(2*p[1]*(1+p[2])))**0.5 
    cond1=y[1]+sig1*(y[0]-1)
    cond2=y[3]+sig2*(y[2]+1)
    return np.array([cond1,cond2])

#********************************************#

#************Fonctions********************#
#integration par methode de Runge-Kutta 4

#definition du probleme differentiel
#j'ai confirme ce quil y avait sur la feuille de Richard
def rhs(y,p):
    #p=y[4:y.size] 
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
    
    return np.array([dy1dx,dy2dx,dy3dx,dy4dx])

#integration Runge-Kutta 4
def rk4(init): 
               
    h = x[1] - x[0]            # stepsize
    y=np.zeros((q,n))          #y=[psi,psi',phi,phi']        
    y[0,0], y[3,0] = y1_0,y4_0 #conditions intiales
    y[1,0], y[2,0] = init[0],init[1] #guess   
    k = 0
    p=init[2:init.size]
    while (k < n-1):
        dydx1 = rhs(y[0:4,k],p)
        dydx2= rhs(y[0:4,k] + 0.5*h*dydx1,p)
        dydx3= rhs(y[0:4,k] + 0.5*h*dydx2,p)
        dydx4= rhs(y[0:4,k] + h*dydx3,p)
        y[0:4,k+1] = y[0:4,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)
        y[4:y.shape[0],k]=init[2:init.size]                        
        k += 1
    if(opti):return y[:,-1]
    else: return y
def F(V):
    """V sont les conditions initiales devinees, pour les N-n2 composantes
    init sont les conditions initiales connues"""
    if(opti):return cible(rk4(V))
    else: print ("active opti mon champion")
                   
#******************************************#     
#//////////////////////PROBLEME/////////////////////////
#potentiel d'Yvan
opti = True
sol=V
for i in range(0,1):
    x = np.linspace(xl,xr,n)
    sol= root(F,sol)
    xr += 1
    n += 50
print(sol)
"""
err=cible(rk4(sol))


origine=sol

if(sol.all()!=origine.all()):
#if(1==1):
    f = open('/home/delta137/stage_2014/2dl/prospect.txt', 'a')
    f.write("\npar: alpha beta gamma delta1 delta2\n")
    f.write("       %d    %d     %d     %d      %d  \n" % (p[0],p[1],p[2],p[3],p[4]))
    f.write("ch:   psi'(0)           phi(0)          err(cond1)          err(cond2)\n")
    f.write("    %e      %e      %e         %e \n" % (sol[0],sol[1],err[0],err[1]))      
    f.close()
#x = np.linspace(xl,xr,n)
#pylab.scatter(x, y[0,:])
#pylab.show()

"""


