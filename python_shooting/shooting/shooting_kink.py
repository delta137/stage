#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)

import numpy
import math
import pylab
graph=True

# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

#integration
#methode de Runge-Kutta 4
def rk4(y1_0, y2_0, lam, m, rhs, xl, xr, n):
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
        
        dy1dx_1, dy2dx_1 = rhs(y1[k], y2[k],lam,m)
        dy1dx_2, dy2dx_2 = rhs(y1[k] + 0.5*h*dy1dx_1, y2[k] + 0.5*h*dy2dx_1,lam,m)
        dy1dx_3, dy2dx_3 = rhs(y1[k] + 0.5*h*dy1dx_2, y2[k] + 0.5*h*dy2dx_2,lam,m)
        dy1dx_4, dy2dx_4 = rhs(y1[k] + h*dy1dx_3, y2[k] + h*dy2dx_3,lam,m)

        y1[k+1] = y1[k] + (h/6.0)*(dy1dx_1 + 2.0*dy1dx_2 + 2.0*dy1dx_3 + dy1dx_4)
        y2[k+1] = y2[k] + (h/6.0)*(dy2dx_1 + 2.0*dy2dx_2 + 2.0*dy2dx_3 + dy2dx_4)
    
        k += 1

    return y1, y2

def rk4nous(y1_0, y2_0, lam, m, rhs, xl, xr, n):
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
        
        dy1dx_1, dy2dx_1 = rhs(y1[k], y2[k],lam,m)
        dy1dx_2, dy2dx_2 = rhs(y1[k] + 0.5*h*dy1dx_1, y2[k] + 0.5*h*dy2dx_1,lam,m)
        dy1dx_3, dy2dx_3 = rhs(y1[k] + 0.5*h*dy1dx_2, y2[k] + 0.5*h*dy2dx_2,lam,m)
        dy1dx_4, dy2dx_4 = rhs(y1[k] + h*dy1dx_3, y2[k] + h*dy2dx_3,lam,m)

        y1[k+1] = y1[k] + (h/6.0)*(dy1dx_1 + 2.0*dy1dx_2 + 2.0*dy1dx_3 + dy1dx_4)
        y2[k+1] = y2[k] + (h/6.0)*(dy2dx_1 + 2.0*dy2dx_2 + 2.0*dy2dx_3 + dy2dx_4)
    
        k += 1

    return y1, y2

def rhs(y1, y2,lam,m):
    """ RHS function.  ici y1 = phi, y2 = phi' 
        Alors: y2' = phi'' = lambda*phi**3 - m**2*phi """

    dy1dx = y2
    dy2dx = lam*(y1**3)-(m**2)*y1
  

    return dy1dx, dy2dx



def analytic(x,lam,m):
    """ analytic solution """
    return (m/lam**(0.5))*numpy.tanh((m/2**(0.5)) * x )

    
# shoot from x = 0 to x -> infini. Pour ca, on prendra x =50. 
#We will do this by selecting a boundary
# value for y2 and use a secant method to adjust it until we reach the
# desired boundary condition at y1(infini)

#///////////////////////////////////////////////
#parametres du probleme
m=1
lam =1
#integration
# bornes et discretisation
xl, xr, n = 0,10,100 
npts=n

# desired tolerance (on va commencer moins rough)
eps = 1e-5

# initial guess
y1_0 = 0.0   # this is the correct boundary condition a x = 0
y2_0 = 0.6  # on veut trouver l'asymptote, eventuellement

# desired right BC, y1(infini)
y1_1_true = m/lam**(0.5)
#//////////////////////////////////////////////

# integrate
y1_old, y2_old = rk4(y1_0, y2_0, lam, m, rhs, xl, xr, n)

x = numpy.linspace(xl,xr,n)
pylab.scatter(x, y1_old, label="premier essai", marker="x",c=colors[0])

# new guess -- we don't have any info on how to compute this yet, so
# just choose something
y2_m1 = y2_0   # store the old guess
y2_0 = 1/(2**(0.5))+1e-5




# Secant loop
dy = 1000.0   # fail first time through

# keep track of iteration for plotting
iter = 1

while (dy > eps) and (iter < 4):
    
    # integrate
    y1, y2 = rk4nous(y1_0, y2_0, lam, m, rhs, xl, xr, n)

    pylab.scatter(x, y1, label="iteration %d" % (iter), marker="x", 
                  c=colors[iter%len(colors)])
    
    # do a Secant method to get 
    # let eta = our current y2(0) -- this is what we control
    # we want to zero f(eta) = y1_1_true(1) - y1_1(eta)
    
    # derivative (for Secant)
    dfdeta = ( (y1_1_true - y1_old[npts-1]) - (y1_1_true - y1[npts-1]) ) / (y2_m1 - y2_0)

    # correction by f(eta) = 0 = f(eta_0) + dfdeta deta 
    deta = -(y1_1_true - y1[npts-1])/dfdeta
    

    y2_m1 = y2_0
    y2_0 += deta

    dy = abs(deta)
    
    y1_old = y1
    y2_old = y2

    iter += 1
    if (iter > 8):
        break

if (dy > eps): print ("youppi, optimisation reussie")

pylab.plot(x, analytic(x,lam,m), color="0.5", label="analytique($\\phi(x) = m/\\sqrt{\\lambda} tanh(m/\\sqrt{2} (x-x_0)$)")
#pylab.figtext(0.15,0.7,"$\\lambda = m =1$")

leg = pylab.legend(loc=2)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)


print (eps)
print (dy)
print (iter)
pylab.xlim(xl, xr)
pylab.savefig("shoot.png")
if (graph):pylab.show("shoot.png")

