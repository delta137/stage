
%initialisations
xl=0;
xr=6;
n =100 ;
init=0.0  ;         
V=1;              
cible =1-exp(-2^(0.5)*xr);


x0=1/sqrt(2);
sol=fzero(F,x0)

