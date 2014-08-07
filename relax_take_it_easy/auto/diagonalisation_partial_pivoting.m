function [ diagA] = diagonalisation_partial_pivoting(A)

% where A=[A,B]
%
% This function receives an matrix A and a vector B. 
% It returns the factored matrix A calculated by Gaussian elimination with
% PARTIAL pivoting. (only rows interchanging)

[nbrrow,nbrcol]=size(A);

for k = 1:nbrrow            % zero out column k
% select the pivot element
amax =  A(k,k);    %le pivot est d'abord choisit comme l'élément sur la diagonale  
	imax = k;
	%jmax = k;
       for i = k : nbrrow     %Recherche du pivot dans la rangée.
			temp = A(i,k) ;%possible pivot...
          if  abs(temp) > abs(amax)   %s'il est > que les autres, on le choisit,
				amax = temp;   %on l'appel amax,
				imax = i;     %imax et jmax sont ses indices.
				%jmax = k;
          end
    
       end
    %Pour mettre le pivot sur la diagonale:
    
	% switch rows imax and k  »» n'affecte pas le vecteur solution
    %                         »» mais B doit aussi changer (d'où le nbrcol)
    
    if imax ~= k %(s'il y a bien eu un switch)
        temp = A(k,:);
		A(k,:) = A(imax,:);
		A(imax,:) = temp;
    end

    
	% now we start the elimination:
     if   abs(A(k,k))~=0 %is not equal  to 0
        % when a(k,k) is zero, column k already 0 below the diagonal
        % calculate the multipliers and store in a(i,k)
           A(k,:)=A(k,:) / amax; 
           temp=A;
      for i = k+1 : nbrrow; % pour toutes les col, sauf la k ième
                                % pour toutes les row (:)
            A(i,:) = temp(i,:) - temp(i,k) * temp(k,:);    
      end
     end
end


%Now the matrix is upper-diaponal.
%Pour la rendre totalement diagonale:
 for k = nbrrow:-1: 2
     temp=A;
 for i = 1 : k-1        % pour toutes les col 
                        % pour toutes les row (:)
     A(i,:) = temp(i,:) - temp(i,k) * temp(k,:)  ;
 end
 end
diagA=A;
