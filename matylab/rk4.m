function [rep] = rk4(V)
xl=0;
xr=6;
n =100; 
init=0.0 ;
x=xl:xr/n:xr;
h=xr/n;
y=zeros(2,n);
y(1,1)=init;
y(2,1)=V;
for k=1:n-1
        dydx1= rhs(y(:,k));
        dydx2= rhs(y(:,k) + 0.5*h*dydx1);
        dydx3= rhs(y(:,k) + 0.5*h*dydx2);
        dydx4= rhs(y(:,k) + h*dydx3);

        y(:,k+1) = y(:,k) + (h/6.0)*(dydx1 + 2.0*dydx2 + 2.0*dydx3 + dydx4) ;                       
        k += 1;

end
   rep=y(1,k);
end