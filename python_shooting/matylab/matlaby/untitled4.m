parametres=0;
conv=1e-2;
nombrediterations=100;
Range_x=[0 10];
Bound_Cond=[0 0 1 -1]; %psi(0), phi'(0), psi(inf), phi(inf)
Range_psip_initial = [-0.5 0.5];
Range_phi_initial = [-0.5 0.5];


[x y]=shooting_method(parametres,conv,nombrediterations,...
                                             Range_x,Bound_Cond,Range_y1_initial )                                         
plot(x,y)