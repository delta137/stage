function [res] = pw(x,pref,haut)
m=2.5;
L=m*pref;
if(x<-L)
    res=0
end
if(abs(x)<L)
    res=haut*(-2*(x/L))
end

if(x>L)
    res=0
end

