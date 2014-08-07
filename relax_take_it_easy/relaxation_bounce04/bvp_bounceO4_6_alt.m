%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour r�soudre le probl�me aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                fronti�res bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
function [g,y,bc] = bvp_bounceO4_6(parametres,Y,X,x)
%--------------------------------------------------------------------------
%Construction du probl�me aux conditions fronti�res bvp_vortices_6.m
%--------------------------------------------------------------------------
%%%Variables � l'entr�e:

  %-> Variables symboliques:
  %   [X,Y]
  %-> Vecteur x, variables ind�pendantes: 
  %   [x]
%--------------------------------------------------------------------------
%%%Variables � la sortie:
  
  %-> Les conditions aux fronti�res (bvp.m): 
  %   [bc] 
  %-> Le syst�me d'�quations diff�rentielles � r�soudre (bvp.m): 
  %   [g]  
  %-> Le guess initial (bvp.m): 
  %   [y] 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
 
%Construction de la matrice NX1 g : d(Y)/dX = g(X,Y)
%Y(0+1)=phi // %Y(1+1)=phi' //%Y(2+1)=chi  //%Y(3+1)=chi' 
%Le systeme ODE est:

g = [ Y(1+1);
      (-1/X)*Y(1+1)+2*Y(0+1)*(Y(0+1)^2-1)*(3*Y(0+1)^2-2*parametres(2)-1)-2*Y(0+1)*((Y(2+1)^2-1)^2-(parametres(3)/4)*(Y(2+1)-2)*(Y(2+1)+1)^2)/(Y(0+1)^2+parametres(1))^2;
      Y(3+1);
      (-1/X)*Y(3+1)+(Y(2+1)^2-1)*(4*Y(2+1)-(3*parametres(3)/4))/(Y(0+1)^2+parametres(1))];


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Guess initial y_{j,k}, for j=0..4 and k=0..M-1 (matrice 4M X 1)
g1=ones(size(x));
g2=zeros(size(x));
del=0.8;
pref=0.1/parametres(2);
h=0.5/(0.2+parametres(2));

y = [ (del/2)*tanh(argu(pref,h,x))-(1-del/2);
    0.5*pref*del*sech(argu(pref,h,x)).*tanh(argu(pref,h,x));
    -g1;
    g2];

% y = [ tanh( x*2*sqrt(1-parametres(2)) ).^(parametres(3)) ; 
%       parametres(3)*(tanh( x*2*sqrt(1-parametres(2)) ).^(parametres(3)-1))*2*sqrt(1-parametres(2)).*sech( x*2*sqrt(1-parametres(2)) ).^2 ;
%       tanh( x*sqrt(2)*parametres(1) ).^(2*parametres(3));
%       2*parametres(3)*( tanh( x*sqrt(2)*parametres(1) ).^(2*parametres(3)-1) ).*( sqrt(2)*parametres(1) ).*sech( x*sqrt(2)*parametres(1) ).^2 ];
% 
%    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Conditions fronti�res du probl�me ( matrice 2NX1 bc):
%Il y a n1 conditons � l'origine x(k=0),         �� n1 premi�res lignes
%et n2=N-n1 conditions � la fin de l'intervale   �� n2 derni�res lignes
%le reste est des z�ros...

bc=[ 0;0;
     Y(1+1);
     Y(2+1)+1.10282;
     Y(1+1)+((1+Y(0+1)).*(1+4*sqrt(2-2*parametres(2))*X))./(2*X);
     Y(2+1)+1;
     0;0];
 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
