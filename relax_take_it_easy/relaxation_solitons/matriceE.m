%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour r�soudre le probl�me aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                fronti�res bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
function [ E ] = matriceE(x,y,k,h,N,parametres,g)
%--------------------------------------------------------------------------
%Construction de la matrice E_{n,k}
%--------------------------------------------------------------------------
%%%Variables � l'entr�e:

  %-> Variables x et y:
  %   [x,y]
  %-> indice k: 
  %   [k]
  %-> Espacement entre chaque x: 
  %   [h]
  %-> Nombre d'�quations diff�rentielles:
  %   [N]
  %-> Parametres du probl�mes bvp.m : 
  %   [parametres]
  %-> Le syst�me d'�quations diff�rentielles � r�soudre : 
  %   [g]  
%--------------------------------------------------------------------------
%%%Variables � la sortie:
  
  %-> Matrice E
  %   [E]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
newvectory(1,N)=0;

n=1:N;

newvectory(n)=(1/2)*(y(n,k+1)+y(n,k+1-1));  %remplacer les y

g=g((1/2)*(x(k+1)+x(k+1-1)),newvectory);    %remplacer les x

E(n,k)=y(n,k+1)- y(n,k+1-1)- h*g(n);        %construire la matrice E

