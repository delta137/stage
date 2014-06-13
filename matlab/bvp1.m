%système d'équations différentielles soit y''=-y, équation harmonique
%S1: y2=y1'
%S2:-y1=y2'
function [ dydx ] = bvp(x,y,p)
    %legende p
    %1:alpha
    %2: beta
    %
    %
    %
    %plein de facteurs un peu longs
    a1=(y(1).^2-1);
    a2=(y(1).^2-p(4));
    a3=(y(3).^2-1);
    a4=(p(3)+y(1).^2);
    a5=y(3)+1;
    a6=y(3)-2;
    a7=(4.*y(3)-3.*p(5)/4);
    
    dy1dx=y(2);
    dy2dx=2.*y(1).*(2.*a1.*a2+a1.^2-p(0).*(a3.^2-(p(5)/4).*a5.*a6)/a4);
    dy3dx=y(4);
    dy4dx=p(1).*a3.*a7/a4;
    %le sys deq diff
    dydx = [dy1dx; dy2dx ;dy3dx ;dy4dx ];
end