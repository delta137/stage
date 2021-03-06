%%%%%%%%%                RELAXATION METHOD                %%%%%%%%%%%%%%%%%
%%%%%%%%%      pour r�soudre le probl�me aux valeurs      %%%%%%%%%%%%%%%%%
%%%%%%%%%                fronti�res bvp.m                 %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
%%%%%%%%%   (Les param�tres pour la r�solution du bvp)    %%%%%%%%%%%%%%%%%

%Cette fonction donne les param�tres pour la r�solution d'un syst�me 
%de N �quations diff�rentielles � conditions fronti�re bvp.m avec 
%relaxation_method.m .
%--------------------------------------------------------------------------

  %-> Nombre maximal d'it�ration:
      nombrediterations=2000;
  %-> Pr�cision du r�sultat: 
      conv=5*10^-8;  
      
  %-> Nombre d'�quations diff�rentielles:
      N=4 ; 
  %-> Nombre de conditions intitiales:
      n1=2;
  %-> Nombre de conditions finales:
      n2=N-n1;
                
%On a donc ces N+1 variables symboliques:
syms y0 y1 y2 y3 X
Y=[y0, y1, y2, y3]; 
   
  %-> Nombre de points:
      M=100;
  %-> Longueur de l'intervalle en x:
      L=80;             
  %-> Alors l'espacement entre chaque x_k:
      h=L/(M-1);       
  %-> On a donc le vecteur x_k (M X 1) (variable ind�pendante):
      x=linspace(0,L,M);       

%Scale factor appropriate to each type of variables to mesures the typical 
%size of each variable:
      scalv=[  1,...
               2*sqrt(1-parametres(2)),...
               1,...
               sqrt(2)*parametres(1)   ];

%Param�tres graphiques:
  
%                *Option graphique: oui(1), non(0) 
                  option_graphique=1   ;                               

%Figure 1: 
  %-> subplot(2,1,1): Solution num�rique de Y(x) avec Relaxation Method

%                *on veut un gaphique des variables d�pendantes Y... 
                  var_graphique = [1,3];

%                *on veut un gaphique avec un range en y de... 
                  axis_graphique_y =[0,max(scalv(var_graphique))];
%                *on veut un gaphique avec un range en x de... 
                  axis_graphique_x = [x(1),x(M)];
                  
%                *on veut un gaphique tout les: ... it�rations 
                  nbr_iteration_graphique = 1;

  %-> subplot(2,1,2): Log10(err) vs le # d'it�rations


