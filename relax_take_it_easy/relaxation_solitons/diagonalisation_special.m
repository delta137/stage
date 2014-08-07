function [ diagA,omega ] = diagonalisation_special(A,B,omega,n1)

% This function receives an matrix, a vector B and a vector omega (
%it will keep track of the col switches).
% It returns the factored matrix A calculated by Gaussian elimination with
% complete pivoting.

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%SPECIAL:
%the pivots elements can only be found amoung the first n1 cols of the matrix
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

A=[A,B];
[nbrrow,nbrcol]=size(A);

for k = 1:nbrrow            % zero out column k
% select the pivot element
amax =  A(k,k);    %le pivot est d'abord choisit comme l'élément sur la diagonale  
	imax = k;
	jmax = k;
    
	for i = k : nbrrow     %Recherche du pivot dans la «sous-matrice» (n-k+1)*(n-k+1).
	for j = k : n1
			temp = A(i,j) ;%possible pivot...
            
			if  abs(temp) > abs(amax)   %s'il est > que les autres, on le choisit,
				amax = temp;   %on l'appel amax,
				imax = i;     %imax et jmax sont ses indices.
				jmax = j;
            end
    end
    end
    %Pour mettre le pivot sur la diagonale:
    
	% switch rows imax and k  »» n'affecte pas le vecteur solution
    %                         »» mais B doit aussi changer (d'où le n+1)
	    j = 1 : nbrcol;
		temp = A(k,j);
		A(k,j) = A(imax,j);
		A(imax,j) = temp;
    
	% switch column jmax and k  »» affecte le vecteur solution dont il faut
	%                           »» changer la row.
    %                           »» mais n'affecte pas B
	    i = 1 : nbrrow;
		temp = A(i,k);
		A(i,k) = A(i,jmax);
		A(i,jmax) = temp;
   
	% change omega to remember the pivots
    %!!s'il y a eu changement de pivot
    %donc seulement pour les n1 premières col
    if k <= n1  %(less than or equal to)
	temp = omega(k);
	omega(k) = omega(jmax);
	omega(jmax) = temp;
    end
    
	% now we start the elimination:
     if  abs(A(k,k))~=0 %is not equal  to 0)
        % when a(k,k) is zero, column k already 0 below the diagonal
 
        % calculate the multipliers and store in a(i,k)
        for j = k : nbrcol
		   A(k,j)=A(k,j) / amax;
        end
      
        temp=A;
		for i = k+1 : nbrrow % pour toutes les col, sauf la k ième
        for j = k : nbrcol % pour toutes les row
         A(i,j) = temp(i,j) - temp(i,k) * temp(k,j);    
        end
        end
			   
     end
end

%Now the matrix is upper-diaponal.
%Pour la rendre totalement diagonale:
 for k = nbrrow:-1: 2
     temp=A;
 for i = 1 : k-1      % pour toutes les row
 for j = 1 : nbrcol   % pour toutes les col   
     A(i,j) = temp(i,j) - temp(i,k) * temp(k,j)  ;
 end
 end
 end
diagA=A;