#librairies
import numpy as np
import math
import pylab
import scipy
from scipy.integrate import ode
from scipy.optimize import broyden1 as broy
from scipy.optimize import fsolve 

colors = ["k", "r", "g", "b", "c"]
#************Configurations********************#
graph=True          #montrer le graph a la fin ou non
p=[1,1,1,1,1]       # param = alpha, beta, gamma, delta1, delta2 
q=4                 #nombre de degres de libertes
xl, xr, n = 0,6,1000#histogramme
x = np.linspace(xl,xr,n)
err_min = 1e-5      #tolerance

#soit y le vecteur y=[psi,psi',phi,phi']
#condition frontieres
#initiales
y1_0,y4_0 = 0.0,0.0     #connues
V=np.array([  1.03251504e-83  , 9.99551219e-01])   #devinettes
origine=V

#a l'infini, eq diff qui decrit deux conditions
def cible(y):
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
    
    return np.array([dy1dx,dy2dx,dy3dx,dy4dx])

#integration Runge-Kutta 4
def rk4(init): 
               
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
        y[:,k+1] = y[:,k] + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4)                        
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
sol= fsolve(F,V)
#opti = False
#y=rk4(sol)
err=cible(rk4(sol))
origine=sol
print(sol)
if(sol.all()!=origine.all()):
#if(1==1):
    f = open('/home/delta137/stage_2014/2dl/prospect.txt', 'a')
    f.write("\npar: alpha beta gamma delta1 delta2\n")
    f.write("       %d    %d     %d     %d      %d  \n" % (p[0],p[1],p[2],p[3],p[4]))
    f.write("ch:   psi'(0)           phi(0)          err(cond1)          err(cond2)\n")
    f.write("    %e      %e      %e         %e \n" % (sol[0],sol[1],err[0],err[1]))    
  
    f.close()
#pylab.scatter(x, y[0,:])
#pylab.show()
"""
iter =1
itermax=20
while (iter < itermax):
    mem=V
    ""shooting selon psi l=(0) et phi l=(1)""
    for i in range(0,2):        
        #1ere integration avec la borne inf de y2_0
        y=rk4(V[:,0]) 
        if(i==0):err1,test=cible(y[:,-1]),cible(y[:,-1])
        else:err1_1,test=cible(y[:,-1]),cible(y[:,-1]
        #si l'erreur est assez petite, on peut rentrer chez nous
        if (abs(test) < err_min):
            err=test
            vraies_ini=V[:,0]
            print ("break")
            break

        #2eme integration avec la borne moy de y2_0
        V[i,0]=(V[i,0]+V[i,1])/2
        z =rk4(V[:,0])
        if(i==0):err1=cible(y[:,-1])
        else:err1_1=cible(y[:,-1])

        #si l'erreur est assez petite, on peut rentrer chez nous
        if (abs(err2) < err_min):
            err=err2
            vraies_ini=V[:,0]
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


