%%%%%%%%%              Choix de guess initial             %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
%Cette fonction est utilis�e lorsque l'on cherche une solution � un bvp
%dont les param�tres changes. On se servira de la premi�re solution comme 
%guess initial pour trouver la prochaine; si cela est avantageux compar� au
%guess initial d�crit dans bvp.m

%--------------------------------------------------------------------------
%%%Variables � l'entr�e:
  %-> Le bvp:
  %   [bvp.m]
  %-> Les parametres pour la r�solution:
  %   [parametres_resolution.m]

%--------------------------------------------------------------------------
%%%Variables � la sortie:
  %-> si la m�thode est efficace ou non:
  %   [en fprintf]

%--------------------------------------------------------------------------
 fprintf('\n****Utilisation de choix_guess_initial: v�rification-->\n')
 
%--------------------------------------------------------------------------
%Renommer le guess initial normal pr�vu dans bvp.m:
guess_initial=y;
%Renomme la solution pr�c�dente obtenue:
K=(j-1)*N+(1:4);
solution_precedente=solutionnum(K,:);

%--------------------------------------------------------------------------
%Pour garder en m�moire le nombre d'it�ration de parametres_resolution.m:
nombrediterations_initial=nombrediterations;
%on veut tester qu'une seule it�ration:
nombrediterations=0;

%--------------------------------------------------------------------------
%Pour garder en m�moire les variables d�pendentes dont on veut un
%graphique par la suite (choix_guess_initial.m ne produit pas de graphique
%des deux it�rationx qu'ils calculent)
var_graphique_initial=var_graphique;
var_graphique = [];
%--------------------------------------------------------------------------
%Faire une premi�re it�ration avec la solution pr�c�dente:
clear err
clear solutionnum

G=g;
BC=bc;
y=solution_precedente;
relaxation_method;
err_solution_precedente=err;
K=(j-1)*N+(1:4);
solution_precedente_1=solutionnum(K,:);

%--------------------------------------------------------------------------
%Faire une premi�re it�ration avec le guess initial pr�vu:
clear err
clear solutionnum

g=G;
bc=BC;
y=guess_initial;
relaxation_method;
err_guess_initial=err;
K=(j-1)*N+(1:4);
guess_initial_1=solutionnum(K,:);
g=G;
bc=BC;

%--------------------------------------------------------------------------
%�valuer quel est le meilleur guess en comparant les erreurs apr�s le
%calcul d'une it�ration:

if err_guess_initial > err_solution_precedente
    
    y=solution_precedente_1;
    fprintf('****choix_guess_initial: Utile!\n')
    
else

    y=guess_initial_1;
    fprintf('****choix_guess_initial: Pas utile...\n')

end

%--------------------------------------------------------------------------
%Redonner le bon nombre d'it�ration pr�vu dans parametres_resolution.m:
nombrediterations=nombrediterations_initial;
%Redonner les variables d�pendantes dont on veut le graphique:
var_graphique=var_graphique_initial;

