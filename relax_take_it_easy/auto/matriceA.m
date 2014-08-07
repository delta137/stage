%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour résoudre le problème aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                frontières bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
function [ A ] = matriceA(x,y,h,N,M,n1,n2,parametres,s,S)

%--------------------------------------------------------------------------
%Construction de la matrice A de dimensions (4Mx4M)
%--------------------------------------------------------------------------
%%%Variables à l'entrée:

  %-> Variables x et y:
  %   [x,y]
  %-> Espacement entre chaque x: 
  %   [h]
  %-> Nombre d'équations différentielles:
  %   [N]
  %-> Nombre de points:
  %   [M]
  %-> Nombre de conditions intitiales:
  %   [n1]
  %-> Nombre de conditions finales:
  %   [n2]
  %-> Parametres du problèmes bvp.m : 
  %   [parametres]
  %->  Matrice S : 
  %   [S] et [s] 
%--------------------------------------------------------------------------
%%%Variables à la sortie:
  
  %-> Matrice A
  %   [A]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
A(N*M,N*M)=0;              %preallocating
%newvectory(N,N*M)=0;       %preallocating
SS(N,2*N)=0;               %preallocating


%%% The first n1XN block of nonzero entries in the top left-hand corner of 
%%% A comes from the boundary condition S_j,n at point k = 0:

J=1:N;

SS=S(x(1),y(J,1)'); %remplacer les x %remplacer les y

    j=n2+1:1:N;
for n=1:1:N
    A(j-n2,n)=SS(j,n+N); %construire A
end


%%% we will have a last n2XN block corresponding to the second 
%%% boundary condition :
clear  SS 
%newvectory(N,N*M)=0;
SS(N,2*N)=0;


SS=S(x(M),y(J,M)'); %remplacer les x %remplacer les y

    j=1:1:n2;
for n=1:1:N
    A(N*M-n2+j,N*(M-1)+n)=SS(j,n); %construire A
end

 
%%% The centers blocks will be the 4X8 S_j^n , there will be M-1 blocks 
%%% like that, just like there is M-1 interior points.
%%% We need to construct this S_j,n block first: fonction matriceS.m

clear S newvectory
newvectory(1,N)=0;
S(N,N*M)=0;

n=1:1:2*N;
for i=0:1:M-2;     
    newvectory(J)=(1/2)*(y(J,i+2)+y(J,i+2-1)); %remplacer les y
    S=s((1/2)*(x(i+2)+x(i+2-1)),newvectory(J)); %remplacer les x 
    A(n1+J+i*N,n+i*N)= S(J,n); 
end
 
      


