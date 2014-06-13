parametres=0;
conv=1e-2;
nombrediterations=100;
Range_x=[0 4.2];
Bound_Cond=[0 1];
Range_y2_initial = [0.76 15];

[x y]=shooting_method3(parametres,conv,nombrediterations,...
                                             Range_x,Bound_Cond,Range_y2_initial )                                         
plot(x,y)