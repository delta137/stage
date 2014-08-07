function [rep] = F(x0)
xr=6;
cible =1-exp(-2^(0.5)*xr);
rep= cible - rk4(x0);
end