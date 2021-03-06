#librairies
import numpy as np
import math
import pylab
import scipy
from scipy.integrate import ode
from scipy.optimize import broyden1 as broy
from scipy.optimize import root

colors = ["k", "r", "g", "b", "c"]
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
    
#fonction a minimiser    
def F(V):
    """V sont les conditions initiales devinees, pour les N-n2 composantes
    init sont les conditions initiales connues"""
    if(opti):return cible(rk4(V))
    else: print ("active opti mon champion")
    
#a l'infini, eq diff qui decrit deux conditions
def cible(y):
    sig1=(8*(1-p[3]))**0.5
    sig2=((p[0]*(16+3*p[4]))/(2*p[1]*(1+p[2])))**0.5 
    cond1=y[1]+sig1*(y[0]-1)
    cond2=y[3]+sig2*(y[2]+1)
    return np.array([cond1,cond2])
                   
#******************************************#     
#************Configurations********************#
graph=True          #montrer le graph a la fin ou non
p=np.zeros(5)
# param = alpha, beta, gamma, delta1, delta2 
p[0],p[1],p[3]=1,1,0.3
p[2],p[4]=0.3,0.3

      
q=4                 #nombre de degres de libertes

#soit y le vecteur y=[psi,psi',phi,phi']
#condition frontieres
#initiales
y1_0,y4_0 = 0.0,0.0     #connues
br=0
youppi=np.zeros([25,1])
youp=0
l=[[0]]
psi=[]
phi=[]
xpp=[]
f_psi_=[0, 0.03]
f_phi =[0.9, 0.995] 
r_g=50
r_w=50

f = open('/home/eric/prospect1.txt', 'a')
f.write("\npar: alpha beta gamma delta1 delta2\n")
f.write("       %d    %d     %d     %d      %d  \n" % (p[0],p[1],p[2],p[3],p[4]))
for s in range(0,r_g*r_w):
    l.append([s+1])
for g in range(0,r_g):
    for w in range(0,r_w):
        prout=False     #est-ce que la recherche de racine a fail?
        xl, xr, n = 0,10,100#histogramme
        V=np.array([f_psi_[0]+w*(f_psi_[1]-f_psi_[0])/(r_w), \
                    f_phi[0]+g*(f_phi[1]-f_phi[0])/(r_g)])   #devinettes
        origine=V


        #//////////////////////PROBLEME/////////////////////////
        opti = True             #rk4 va renvoyer y final
        #iteration: une fois qu'une solution interessante est trouvee, je verifie si 
        #elle va converger quand j'augmente le range de x. Je plugge comme guess initial
        #la solution trouvee avec la premiere iteration (si elle existe) 
       
        x = np.linspace(xl,xr,n)
        sol = root(F,V)
        print(sol.x)                #work in progress :v
        if(sol.x==V).all():
            br += 1
            prout=True      #fail de la recherche de racines
            print("breaky %d" %(br))
        else:   
            err=cible(rk4(sol.x))        
            print("er %e, %e" %(err[0],err[1]))
            print("-----")
            if((err<=1).all()) : 
                l[youp].append(V)
                youp+=1          
                xpp.append(err[0])
                psi.append(sol.x[0])
                phi.append(sol.x[1])
                f.write("====================================\n" )  
                f.write("ch:   psi'(0)           phi(0)          err(cond1)          err(cond2)\n")
                f.write("    %e      %e      %e         %e \n" % (sol.x[0],sol.x[1],err[0],err[1]))
  
