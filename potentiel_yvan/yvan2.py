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
param=[1,1,1,1,1]

q=4 #nombre de degres de liberts

#histogramme
xl, xr, n = 0,4,100
npts=n

#condition initiales psi: y, phi: z
y1_0 = 0.0        #connu
y2_0 = [0.3, 0.9] #plage de valeurs qu'on va etudier
z1_0=[phi_plus-eps, phi_plus]
z2_0=0            #on ne bouge pas initialement en phi
#cible du shooting (y1(xr))
cible = param[1]/param[0]**(0.5) *(1-math.exp(-2**(0.5)*xr)) 

#tolerance
err_min = 1e-5
#********************************#

#************Fonctions********************#
#integration par methode de Runge-Kutta 4
def rk4(init, param, rhs, xl, xr, n,q):
    """ R-K 4 integration:
           y1_0 and y2_0 are y1(0) and y2(0)
           rhs is the righthand side function
           xl and xr are the domain limits
           n is the number of integration points """

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
    """ RHS function.  ici y1 = psi, y2 = psi' 
        z1 = phi, z2 = phi' 
        Alors: y2' = phi'' = lambda*phi**3 - m**2*phi """
    #plein de facteurs un peu longs
    a1=(y[0]**2-1)
    a2=(y[0]**2-param[3])
    a3=(y[2]**2-1)
    a4=(param[2]+y[0]**2)
    a5=y[2]-1
    a6=y[2]+2
    a7=(4*y[2]+3*param[4]/4)
    
    #le sys deq diff
    dydx = numpy.array([y[1],
    [2*y[0]*(2*a1*a2+a1**2-param[0]*(a3+(param[4]/4)*a5*a6)/a4)],
    [y[3]],
    [param[0]*a3*a7/a4]]) 

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
itermax=10
while (iter < itermax):
    #--- premiere fonction ---#
    #1ere integration avec la borne inf de y2_0
    init=[y1_0,y2_0[0]]
    y = rk4(init, param, rhs, xl, xr, n,q)
    err1= y[0,-1]-cible
    
    #si l'erreur est assez petite, on peut rentrer chez nous
    if (abs(err1) < err_min):
        err=err1
        y2_vrai=y2_0[0]
        print ("break")
        break

    #2eme integration avec la borne moy de y2_0
    init=[y1_0, (y2_0[0]+y2_0[1])/2]
    z= rk4(init, param, rhs, xl, xr, n,q)
    err2= z[0,-1]-cible

    #si l'erreur est assez petite, on peut rentrer chez nous
    if (abs(err2) < err_min):
        err=err2
        y2_vrai=y2_0[1]
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
            pylab.scatter(x, z[0,:], label="iteration %d" % (iter), marker="x", 
                    c=colors[iter%len(colors)])
	


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


