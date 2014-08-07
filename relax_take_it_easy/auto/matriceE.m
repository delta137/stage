%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour résoudre le problème aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                frontières bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
function [ E ] = matriceE(x,y,k,h,N,parametres,g)
%--------------------------------------------------------------------------
%Construction de la matrice E_{n,k}
%--------------------------------------------------------------------------
%%%Variables à l'entrée:

  %-> Variables x et y:
  %   [x,y]
  %-> indice k: 
  %   [k]
  %-> Espacement entre chaque x: 
  %   [h]
  %-> Nombre d'équations différentielles:
  %   [N]
  %-> Parametres du problèmes bvp.m : 
  %   [parametres]
  %-> Le système d'équations différentielles à résoudre : 
  %   [g]  
%--------------------------------------------------------------------------
%%%Variables à la sortie:
  
  %-> Matrice E
  %   [E]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
newvectory(1,N)=0;

n=1:N;

newvectory(n)=(1/2)*(y(n,k+1)+y(n,k+1-1));  %remplacer les y

g=g((1/2)*(x(k+1)+x(k+1-1)),newvectory);    %remplacer les x

E(n,k)=y(n,k+1)- y(n,k+1-1)- h*g(n);        %construire la matrice E

