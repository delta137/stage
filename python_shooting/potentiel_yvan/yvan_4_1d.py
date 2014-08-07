#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)
#version 3: condition frontiere ajustable
#Yvan: j'ai implemente des arrays, plus faciles avec plus d'equations
#Yvan3: apres multiples essais avec Jacobien et fonction 
#je reviens a mon code qui marchait bien lol

#librairies
import numpy as np #tu peux nommer tes fonctions comme tu veux
import math     
import pylab    #pour les graphiques?
import scipy
from scipy.integrate import odeint
from scipy.optimize import broyden1,fsolve as broy,fsolve

# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

#-------Configurations-------
graph=True #montrer le graph a la fin ou non

# -------param = lambda, m  (on les retrouve dans le kink) -------
param=[1,1]

q=2 #nombre de degres de libertes, utile quand tu veux deux champs, q->4

#histogramme, range et discretisation d'integration
xl, xr, n = 0,5,100 #5 est notre infini
npts=n

opti = False 

#condition initiales, dimension N-n1
init=np.array(0.0)

#de dimension N-n2
cible = np.array(param[1]/param[0]**(0.5) *(1-math.exp(-2**(0.5)*xr))) 
#********************************#

#************Fonctions********************#

#integration par methode de Runge-Kutta 4
def rk4(V):
    """V sont les conditions initiales devinees, pour les N-n2 composantes
    init sont les conditions initiales connues"""
           #R-K 4 integration:
           #y1_0 and y2_0 are y1(0) and y2(0)
           #rhs is the righthand side function
           #xl and xr are the domain limits
           #n is the number of integration points 

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
        
        dydx1 = rhs(y[:,k])
        dydx2= rhs(y[:,k] + 0.5*h*dydx1)
        dydx3= rhs(y[:,k] + 0.5*h*dydx2)
        dydx4= rhs(y[:,k] + h*dydx3)

        y[:,k+1] = y[:,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)                        
        k += 1
   if (opti):
      return y[0,k]
   else:
      return y



#difference entre le resultat de l'integration et la cible        
#-------fonction a minimiser-------
def F(V):
    """V sont les conditions initiales devinees, pour les N-n2 composantes
    init sont les conditions initiales connues"""
    
    return cible- rk4(V)      
    
       

    
#------definition du probleme differentiel-------
def rhs(y,x):
    """ RHS function.  ici y[0,:] = psi, y[1,:] = psi' 
        y[2,:]= phi, y[3,:] = phi' """
    d1=y[1]
    d2=param[0]*(y[0]**3)-(param[1]**2)*y[0]
    dydx = np.asfarray([d1,d2], dtype='float') 

    return dydx
    
    
    
#-------soln analytique (kink)-------
def analytic(x,param):
    """ analytic solution """
    return (param[1]/param[0]**(0.5))*np.tanh((param[1]/2**(0.5)) * x )
               
#******************************************#
#//////////////////////PROBLEME/////////////////////////
#on va shooter dans un range de conditions intiales,
#la conditon pour la vitesse initiale est inconnue, mais 
#une reponse optimale est trouvee en comparant avec la 
#deuxieme condition frontiere (pas initiale)

#////#boucle d'integration, pour la victoire!/////
x = np.linspace(xl,xr,n)
y1=scipy.integrate.odeint(rhs,[0,0.707],x)
pylab.scatter(x, y1[:,0])
#  pylab.scatter(x, y2[0,:])
pylab.show()
    
"""
opti=True
sol= broy(F,V)
opti = False
y=rk4(sol)
pylab.scatter(x, y[:,0])
#pylab.show()
"""
"""
itermax=20
while (iter < itermax):
    #1ere integration avec la borne inf de y2_0
    y =rk4(V)
    err1= y[0,-1]-cible
    
    #si l'erreur est assez petite, on peut rentrer chez nous
    if (abs(err1) < err_min):
        err=err1
        y2_vrai=y2_0[0]
        print ("break")
        break

    #2eme integration avec la borne moy de y2_0
    init=[y1_0, (y2_0[0]+y2_0[1])/2]
    z =rk4(init, param, rhs, xl, xr, n,q)
    err2= z[0,-1]-cible

    #si l'erreur est assez petite, on peut rentrer chez nous
    if (abs(err2) < err_min):
        err=err2
        y2_vrai=(y2_0[0]+y2_0[1])/2
        print ("break")
        break
    
    #UNDERSHOOT:
    if (err1*err2) > 0: #decrit le cas de undershoot
        y2_0[0] = (y2_0[0]+y2_0[1])/2 #on eleve la borne inf

    #OVERSHOOT:
    else:
        y2_0[1] = (y2_0[0]+y2_0[1])/2 # on abaisse la borne sup

    iter += 1
    
    
    if(iter >= (itermax-5)): #trace les 5 dernieres iterations

    #plot moi ca Rejean, le meilleur des deux essais
        if (abs(err1) < abs(err2)):
            pylab.scatter(x, y[0,:], label="iteration %d" % (iter), marker="x", 
                    c=colors[iter%len(colors)])
        else:
            pylab.scatter(x, z[0,:], label="iteration %d" % (iter), marker="x",c=colors[iter%len(colors)])
	


pylab.plot(x, analytic(x,param), color="0.5", label="analytique($\\phi(x) = m/\\sqrt{\\lambda} tanh(m/\\sqrt{2} (x-x_0)$)")
#pylab.figtext(0.15,0.7,"$\\lambda = m =1$")

leg = pylab.legend(loc=2)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)


print (err1)
print (iter)
pylab.xlim(xl, xr)
pylab.savefig("shoot.png")
if (graph):pylab.show("shoot.png")
"""


