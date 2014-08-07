%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour résoudre le problème aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                frontières bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
function [ B ] = matriceB(x,y,h,N,M,n1,n2,parametres,bc,g)
%--------------------------------------------------------------------------
%Construction de la matrice B
% ( matrice -E(j,k) où j donne le # de la variable et k le x(k) )
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
  %-> Les conditions aux frontières de bvp.m : 
  %   [bc] 
  %-> Le système d'équations différentielles à résoudre : 
  %   [g]  
%--------------------------------------------------------------------------
%%%Variables à la sortie:
  
  %-> Matrice B
  %   [B]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
B(N*M,1)=0; %preallocating

%the right 4MX1 matrix can be calculated by giving an initial guess
%for each y_{j,k}.
%--------------------------------------------------------------------------
%%%Pour les conditions frontières:

newvectory(1,N)=0; %preallocating
BC(2*N,1)=0;       %preallocating

j=1:N;

  %-> n1 conditions frontières:
newvectory(j)=y(j,1); %remplacer les y
BC=bc(x(1),newvectory);%remplacer les x

n=N-n1+1:1:N;
B(n-N+n1)=-BC(n); %construire B


  %-> n2 conditions frontières:
newvectory(j)=y(j,M); %remplacer les y
BC=bc(x(M),newvectory);%remplacer les x

n=N+1:1:N+n2;
B(N*(M-1) +n-N +n1)=-BC(n); %construire B

%--------------------------------------------------------------------------
%%%Pour un point intérieur:

E(N,M-1)=0; %preallocating

for k=0:1:M-2
    
E=matriceE(x,y,k+1,h,N,parametres,g) ;

B(n1+j+k*N)=-E(j,k+1) ; %construire B
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%