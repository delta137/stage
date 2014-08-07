#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)

#librairies
import numpy
import math
import pylab

# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

#conditions du code
#output
graph=True
#1= methode eta, 2= methode marie-lou
mode=1


# param[0] =lam, param[1]=m
param=[1,1]

#************Fonctions********************#
#integration par methode de Runge-Kutta 4
def rk4(y1_0, y2_0, param, rhs, xl, xr, n):
    """ R-K 4 integration:
           y1_0 and y2_0 are y1(0) and y2(0)
           rhs is the righthand side function
           xl and xr are the domain limits
           n is the number of integration points """

    x = numpy.linspace(xl, xr, n)
    h = x[1] - x[0]  # stepsize
    
    #def pente(x): return (x-xl)/(xr-xl)
    #y1 = map(pente, range(xl, xr))
    y1 = numpy.zeros((n), dtype = float)
    y2 = numpy.zeros((n), dtype = float)

    # left boundary initialization
    y1[0] = y1_0
    y2[0] = y2_0

    k = 0
    while (k < n-1):
        
        dy1dx_1, dy2dx_1 = rhs(y1[k], y2[k],param)
        dy1dx_2, dy2dx_2 = rhs(y1[k] + 0.5*h*dy1dx_1, y2[k] + 0.5*h*dy2dx_1,param)
        dy1dx_3, dy2dx_3 = rhs(y1[k] + 0.5*h*dy1dx_2, y2[k] + 0.5*h*dy2dx_2,param)
        dy1dx_4, dy2dx_4 = rhs(y1[k] + h*dy1dx_3, y2[k] + h*dy2dx_3,param)

        y1[k+1] = y1[k] + (h/6.0)*(dy1dx_1 + 2.0*dy1dx_2 + 2.0*dy1dx_3 + dy1dx_4)
        y2[k+1] = y2[k] + (h/6.0)*(dy2dx_1 + 2.0*dy2dx_2 + 2.0*dy2dx_3 + dy2dx_4)
    
        k += 1

    return y1, y2

#definition du probleme differentiel
def rhs(y1, y2,param):
    """ RHS function.  ici y1 = phi, y2 = phi' 
        Alors: y2' = phi'' = lambda*phi**3 - m**2*phi """

    dy1dx = y2
    dy2dx = param[0]*(y1**3)-(param[1]**2)*y1
  

    return dy1dx, dy2dx

#soln analytique (kink)
def analytic(x,param):
    """ analytic solution """
    return (param[1]/param[0]**(0.5))*numpy.tanh((param[1]/2**(0.5)) * x )
#******************************************#
    
# shoot de x = 0 a x -> infini. 
# value for y2 and use a secant method to adjust it until we reach the
# desired boundary condition at y1(infini)

#///////////////////////////////////////////////
#parametres du probleme
#integration
# bornes et discretisation
xl, xr, n = 0,10,100
npts=n


y1_0 = 0.0        #connu
y2_0 = [0.3, 0.9] #plage de valeurs qu'on va etudier
beta=1
cible = beta*param[1]/param[0]**(0.5) #cible du shooting (y1(xr))

#tolerance
err_min = 1e-5

#////#boucle d'integration, pour la victoire!/////
x = numpy.linspace(xl,xr,n)
iter =1
while (iter < 30):

    #1ere integration avec la borne inf de y2_0
    y1, y2 = rk4(y1_0, y2_0[0], param, rhs, xl, xr, n)
    err1= y1[-1]-cible

    if (abs(err1) < err_min):
        err=err1
        y2_vrai=y2_0[0]
        print ("break")
        break

    #2eme integration avec la borne moy de y2_0
    y3, y4 = rk4(y1_0, (y2_0[0]+y2_0[1])/2, param, rhs, xl, xr, n)
    err2= y3[-1]-cible

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
    if(iter >= 25):
    #plot moi ca Rejean
        if (abs(err1) < abs(err2)):
            pylab.scatter(x, y1, label="iteration %d" % (iter), marker="x", 
                    c=colors[iter%len(colors)])
        else:
            pylab.scatter(x, y3, label="iteration %d" % (iter), marker="x", 
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

