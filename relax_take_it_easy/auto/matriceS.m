function [ S ] = matriceS(X,Y,h,N,parametres,g)
%--------------------------------------------------------------------------
%Construction de la matrice NX2N S_{j,n}
%--------------------------------------------------------------------------
%%%Variables à l'entrée:

  %-> Variables symboliques:
  %   [X,Y]
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
  
  %-> Matrice S
  %   [S]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for j=1:1:N;
for n=1:1:N;
 S(j,n)   =-KronD(j,n)-(h/2)*diff(g(j),Y(n));
 S(j,n+N) = KronD(j,n)-(h/2)*diff(g(j),Y(n));
end
end


    

