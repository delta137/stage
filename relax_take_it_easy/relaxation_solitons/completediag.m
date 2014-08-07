function [ DELTAY ] = completediag(A,B,N,n1,n2,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN ELIMINATION with
%COMPLETE PIVOTING

% This function receives a matrix A and a vector B. 
% It returns the solution x=deltaY calculated by a particular Gaussian 
% elimination with complete pivoting. If the matrix is singular x will 
% have components equal to infinity (inf) or not a number (NaN).

%The entire procedure, except the backsubstitution step, operates
%only on one block of the matrix at a time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize vectors to keep track of the pivot elements
for i = 1:N*M
	omega(i) = i  ;     % omega will keep track of the col switches
end

[ diagA,omega ] = diagonalisation(A,B,omega);
deltaY=diagA(:,N*M+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute the deltaY vector to get the solution DELTAY 
%(this step is required if there were any column switches)
for i = 1 : M*N
	DELTAY( omega(i) ) = deltaY(i);
end
