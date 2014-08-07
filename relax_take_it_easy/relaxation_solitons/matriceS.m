function [ S ] = matriceS(X,Y,h,N,parametres,g)
%--------------------------------------------------------------------------
%Construction de la matrice NX2N S_{j,n}
%--------------------------------------------------------------------------
%%%Variables � l'entr�e:

  %-> Variables symboliques:
  %   [X,Y]
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


    

