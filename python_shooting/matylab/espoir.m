conv=1e-2;
parametres=0;
options = odeset('RelTol',conv);

[x,y] = ode45(@(x,y)bvp_0(x,y,parametres),...
              [0 4.2],...
              [0   0.76],...
              options ) ;
