%système d'équations différentielles soit y''=-y, équation harmonique
%S1: y2=y1'
%S2:-y1=y2'
function [ dydx ] = bvp_0(x,y,p)
    %legende p
    %1:alpha
    %2: beta
    %
    %
    %
    %plein de facteurs un peu longs
    dy1dx=y(2);
    dy2dx=y(1).^3-y(1);
    %le sys deq diff
    dydx = [dy1dx ;dy2dx];
end