function [res] = pw(x,pref,haut)
m=2.5;
L=m*pref;
if(x<-L)
    res=-1
end
if(abs(x)<L)
    res=haut*(-(x/L).^2+1)-1
end

if(x>L)
    res=-1
end

