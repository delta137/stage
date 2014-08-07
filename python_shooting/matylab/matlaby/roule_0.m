parametres=0;
conv=1e-2;
nombrediterations=500;
Range_x=[0 5];
Bound_Cond=[0 1];
Range_y2_initial = [0.70 0.73];

[x y]=shooting_method3(parametres,conv,nombrediterations,...
                                             Range_x,Bound_Cond,Range_y2_initial )                                         
plot(x,y)