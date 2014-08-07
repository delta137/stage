start_time=tic;
%profile on; 
%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour r�soudre le probl�me aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                fronti�res bvp.m                 %%%%%%%%%%%%%%%%%

%Cette fonction solutionne un syst�me de N �quations diff�rentielles � 
%conditions fronti�re bvp.m avec les param�tres pour la r�solution de
%parametres_resolution.m

% La m�thode:
%We have a linear system of the form A*X = B witch can be solve for
%X=Delta y_{j,k} with the Matlab program X = linsolve(A,B) . 
%The coefficients of the left 4MX4M matrix and the right 
%4X1 matrix can be calculated by giving an initial guess 
%for each y_{j,k}. The solution will be y_{j,k}+Delta y_{j,k}.

%--------------------------------------------------------------------------
%%%Variables � l'entr�e:

  %-> fonctions utilis�es:
  %   bvp.m
  %   matriceS.m
  %   matriceA.m
  %   matriceB.m
  %   matriceE.m
  %   gauss.m ou linsolve
  
  fprintf('RElAXATION METHOD\n')
  fprintf('Nbr d''�quations diff:        %d \n',N)
  fprintf('Nbr de cond intitiales:      %d \n',n1)
  fprintf('Nbr de cond finales:         %d \n',n2)
  fprintf('pr�cision du r�sultat:       %d \n',conv)
  fprintf('Nombre de points:            %d \n',M)
  fprintf('Longueur de l''intervalle:    %d \n',L)
  fprintf('Espacement entre les points: %d \n',h)
  fprintf('Nbr d''it�rations maximal:    %d \n',nombrediterations) 
  parametres
  %fprintf('Param�tres: %d \n',parametres)
  fprintf('........................................\n') 
%--------------------------------------------------------------------------
%%%Variables � la sortie:
  
  %-> La solution trouv�e:
  %   [solutionnum]
  %-> Le nombre de fois o� �slowc� est utilis� (avec les �!�), si la 
  %   m�thode a converg�e ou pas, le # d'int�ration et avec quel 
  %   programme (gauss ou linsolve):
  %   [en fprintf]
  %-> Graphiques:
  %   [Solution y vs x  et  Log10(err) vs le # d'it�rations]
  %-> Temps total de la m�thode:
  %   [total_time]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Le �Boundary Value Problem� est:
%[g,y,bc] = bvp(parametres,Y,X,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construction de la matrice S pour ce probl�me:
s=matriceS(X,Y,h,N,parametres,g);
s=matlabFunction(s,'vars',{X,Y});       %transforme l'expression symbolique 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prealocating:
solutionnum(N*(nombrediterations+1),M)=0;
Delta(N*(nombrediterations+1),M)=0;
%err(1,nombrediterations+1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quelques calculs symboliques de d�part (pour construire A et B):

%pour les conditions fronti�res: (bc)
clear S
%Construction de S aux fronti�res:
for j=n2:1:N-1 
    for n=1:1:N
    S(j+1,n+N)=diff(bc(j+1),Y(n));
    end
end

for j=N+1:1:N+n2
    for n=1:1:N
    S(j-N,n)=diff(bc(j),Y(n));
    end
end

S=matlabFunction(S,'vars',{X,Y});       %transforme l'expression symbolique

%Contruction de la matrice contenant les conditions fronti�res pour B:
bc=matlabFunction(bc,'vars',{X,Y});     %transforme l'expression symbolique 

%Contruction de la matrice contenant les �q diff pour B:
g=matlabFunction(g,'vars',{X,Y});       %transforme l'expression symbolique


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D�but des it�rations:
for j=1:nombrediterations+1 
    
    if j<1.5 % premi�re relaxation � partir du guess initial
     if option_graphique==1   
     %Graphique des guess initiaux:
      figure(1)
      hold on
      for indice =1:length(var_graphique)
      plot( (x), (y(var_graphique(indice),:)) ,'m--')
      end
      title('Solution num�rique de Y(x) avec Relaxation Method');
      xlabel('x');
      axis([ axis_graphique_x(1) axis_graphique_x(2) axis_graphique_y(1) axis_graphique_y(2) ])
     end
      
    else   % it�rations � partir de la solution num�rique pr�c�dente
      y(N,M)=0;
      i=1:N;
      y(i,:)=solutionnum((j-2)*N+i,:);
      
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Construction de la matrice A
%start_time_A=tic;
A(N*M,N*M)=0;
A=matriceA(x,y,h,N,M,n1,n2,parametres,s,S);
%Cond=cond(A);
%if Cond > 10^16
%    b-parametres(2)/5;reak
%end    
%maximum=max(max(A))
%[I,J] = find(A == max(max(A)))
%minimum=min(min(A))
%[I,J] = find(A == min(min(A)))
%time_A=toc(start_time_A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construction de la matrice B
%start_time_B=tic;
B(N*M,1)=0;
B=matriceB(x,y,h,N,M,n1,n2,parametres,bc,g);
%time_B=toc(start_time_B)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R�solution de A*delta(y)=B:

if M<=350
%Solve the linear system with the Matlab program X=linsolve(A,B)=delta(y)
%start_time_linsolve=tic;
deltay = linsolve(A,B)'; 
%time_linsolve=toc(start_time_linsolve)
method=1;

else
%Solve the linear system with X=gauss(A,B)=delta(y)
%start_time_gauss=tic;
deltay=gauss(A,B,N,n1,n2,M);
%time_gauss=toc(start_time_gauss)
method=2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcul de err:
%We compute a value for the average correction err by summing the absolute
%value of all corrections, weighted by a scale factor appropriate to each
%type of variable:
for n=1:N
k=0:1:M-1;
Delta(n+(j-1)*N,1+k)=deltay(n+N*k);
end
for n=1:N
eRR(n)=sum(abs(Delta(n+(j-1)*N,:)),2,'double')/scalv(n);
end
err(1,j)=sum(eRR,'double')/(M*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Utilisation de �slowc�:
%The number slowc modulates the application of corrections. After each
%iteration we apply only a fraction of the corrections found by matrix
%inversion: thus,when err>slowc, only a fraction of the corrections are
%used, but when err<slowc the entire correction gets applied.
slowc=1;
if err(1,j)<slowc %the entire correction gets applied
  for n=1:N
  k=0:1:M-1;
  solutionnum(n+(j-1)*N,1+k)=y(n,1+k)+deltay(n+N*k);
  end
else %err>slowc, only a fraction of the corrections are used
  for n=1:N 
  k=0:1:M-1;
  solutionnum(n+(j-1)*N,1+k)=y(n,1+k)+deltay(n+N*k)/(err(1,j));
  end
  fprintf('!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphiques:

if err(1,j)>conv 
    if option_graphique==1   
   %Graphique des it�rations, sauf la derni�re:
    if mod(j,nbr_iteration_graphique)==0 %(si j est un multiple entier de nbr_iteration_graphique )
    figure(1)
    hold on
      for indice =1:length(var_graphique)
      plot(x,solutionnum((j-1)*N+var_graphique(indice),:),'c:');
      end
    end
    end
   %Graphique de la derni�re it�ration s'il n'y a pas eu convergence de la m�thode:     
    if j> nombrediterations+0.5
    fprintf('\n........................................\nIt�ration #%d : la m�thode n''a pas converg�e... ',j-1)
    
    if option_graphique==1   
    figure(1)
    hold on
      for indice =1:length(var_graphique)
      plot(x,solutionnum((j-1)*N+var_graphique(indice),:),'r:');
      end
   %Graphique err vs it�ration
    J=1:j;
    figure(2)
    plot(J,log10(err(J)),'-*','MarkerSize',3)
    title(' Log10(err) vs le # d''it�rations ');
    xlabel('# it�ration');
    ylabel('log10(err)');
    end
    
    
    end
     
    
else
    if option_graphique==1   
   %Graphique de la derni�re it�ration s'il y a eu convergence de la m�thode:
    figure(1)
    hold on
     for indice =1:length(var_graphique)
     plot(x, (solutionnum((j-1)*N+var_graphique(indice),:)) ,'b-'  ); %h
     end
   %Graphique err vs it�ration
    J=1:j;
    figure(2)
    plot(J,log10(err(J)),'-*','MarkerSize',3)
    title(' Log10(err) vs le # d''it�rations ');
    xlabel('# it�ration');
    ylabel('log10(err)');
    end
    fprintf('\n........................................\nIt�ration #%d : la m�thode a converg�e! ',j-1)
    
 break %La boucle est bris�e
 
end 
%pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Et la m�thode utilis�e fut:   
if method==1
fprintf('(avec l''utilisation de �linsolve�)\n')
else
fprintf('(avec l''utilisation de �gauss.m�)\n')
end 

total_time=toc(start_time);
fprintf('En un temps total de: %d \n',total_time)

if option_graphique==1
figure(4)
plot(x,y(1,:))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile viewer; 