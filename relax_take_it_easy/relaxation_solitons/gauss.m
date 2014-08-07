function [ DELTAY ] = gauss(A,B,N,n1,n2,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN ELIMINATION with
%COMPLETE or PARTIAL PIVOTING  

% This function receives a matrix A and a vector B. 
% It returns the solution x=deltaY calculated by a particular Gaussian 
% elimination with complete pivoting. If the matrix is singular x will 
% have components equal to infinity (inf) or not a number (NaN).

%The entire procedure, except the backsubstitution step, operates
%only on one block of the matrix at a time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%réunir en une même matrice A et B
temp(N*M,N*M+1)=0;
%temp=[A,B];
i=1:N*M;
temp(:,N*M+1)=B;
temp(:,i)=A;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize vectors to keep track of the pivot elements
omega_bef = 1:N*M  ;     % omega will keep track of the col switches

%INITIAL CONDITION BLOCK: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ diagA,omega_aft ] = diagonalisation(temp(1:n1,1:N*M+1),omega_bef(1:N));

    %if there were any column switches during diagonalisation:  
        for i = 1:N
        if omega_bef(i)~=omega_aft(i) 
        [column] = find(omega_aft==omega_bef(i));
        %permute the col of the last temp
        temp( :,column) = temp(:,i);
        %fprintf('initial block column switche at %d\n',i)
        end    
        end
        
    %Reconstruire omega:
        i= 1:N;
        omega(N*M,1)=0;
        omega=omega_bef;
        omega(i)=omega_aft;
    %reconstruire A|B
        i=1:n1;
        temp(i,:)=diagA;


%The CENTERS BLOCKS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the 4X8 S_j^n , there is M-1 blocks like that, 
%just like there is M-1 interior points.

for m=0:M-2 %For each block: (il y en a M-1)
        
    %Reduce to zero the N*n1 left part of the S block using the
    %last n1*n1 identity block.
    
        i =  1 + n1+m*N : N  + n1+m*N ;   % pour toutes les rangées de S
        k =  1 + m*N    : n1  +m*N ;
        j =   N*M+1   :-1 : 1+m*N ;       % pour toutes les col de AB à droite de la partie N*n1.         
        temp(i, j ) = temp(i, j ) - temp(i, k)*temp(k,j) ;  
        
        
        
    %Diagonaliser le N*N block du milieu (à la droite des N*n1 zéros):
        jj=m*N+n1+1:M*N+1;
        ii=m*N+n1+1:(m+1)*N;
        ABtemp=temp(i,jj);
        %omega_bef=omega(ii);
      
        [ diagblock] = diagonalisation_partial_pivoting( ABtemp);

    %if there were any column switches during diagonalisation:  
    
      %NO COL SWITCHES->> using diagonalisation_partial_pivoting
        %for i = n1+m*N+1:(m+1)*N
        %if omega_bef(i-n1-m*N)~=omega_aft(i-n1-m*N) 
        %[column] = find(omega_aft==omega_bef(i-n1-m*N));
        %permute the col of the last temp
        %temp( :,i+ column-1) = temp(:,i);
        %fprintf('center block column switche at %d\n',i)
        %end    
        %end
        
    %Reconstruire AB:
        temp(i,jj)=diagblock;     
    %Reconstruire omega:
        %omega(ii)=omega_aft;

  
end

%FINAL CONDITION BLOCK: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reduce to zero the n2*n1 left part of the last n2*N block using the
%last n1*n1 identity block.
         i =  1 +n1+(M-1)*N : M*N ; % pour toutes les rangées
         k =  1 +(M-1)*N    : n1 +(M-1)*N ;
         j =  N*M +1: -1 :(M-1)*N+1;           % pour toutes les col 
         temp(i, j ) = temp(i, j ) - temp(i,k)*temp(k,j) ; 
       

    %Diagonaliser le n2*n2 block de droite:
        ABtemp    =temp ( (M-1)*N+n1+1:M*N ,(M-1)*N+n1+1:length(temp) );
        %omega_bef =omega( M*N-n2+1:M*N );
        [ diagblock] = diagonalisation_partial_pivoting( ABtemp);

    %NO COL SWITCHES->> using diagonalisation_partial_pivoting
    %if there were any column switches during diagonalisation: 
        %for i = M*N-n2+1:M*N
        %if omega_bef(i-M*N+n2)~=omega_aft(i-M*N+n2) 
        %[column] = find(omega_aft==omega_bef(i-M*N+n2));
        %permute the col of the last temp
        %temp( :,i+column-1) = temp(:,i); 
        %fprintf('final block column switche\n')
        %end    
        %end   
         
        j=(M-1)*N+n1+1:M*N+1;
    %Reconstruire omega:  
        %omega(i)=omega_aft;
    %Reconstruire AB:
    %Target structure of the Gaussian elimination
        temp(i,j)=diagblock;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The matrix is now factored, now solve the system.

%BACKSUBSTITUTION PROCEDURE:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaY=zeros(N*M,1);
%on commence par la fin:
%FINAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = N*M:-1:N*M-n2+1;
deltaY(i)=temp(i,N*M+1);

%CENTERS BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CENTERS BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1:n1;
for m=M-2:-1:0 %For each block: (il y en a M-1)
    i = N+m*N+n1:-1:1+m*N+n1; %for each rows
    somme     = temp(i,(m+1)*N+n1+j)*deltaY((m+1)*N+n1+j);
    deltaY(i) = temp(i,M*N+1) - sum(somme,2,'double');
           
end
%INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = n1:-1:1;
    somme     = temp(i,n1+j)*deltaY(n1+j);
    deltaY(i) = temp(i,M*N+1) - sum(somme,2,'double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute the deltaY vector to get the solution DELTAY 
%(this step is required if there were any column switches)
 i = 1 : M*N;
 DELTAY( omega(i) ) = deltaY(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ABtarget(:,1:N*M)*deltaY-ABtarget(:,N*M+1)%...backsubstitution step verification 

