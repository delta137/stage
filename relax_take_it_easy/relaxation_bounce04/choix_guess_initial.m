%%%%%%%%%              Choix de guess initial             %%%%%%%%%%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%
%Cette fonction est utilisée lorsque l'on cherche une solution à un bvp
%dont les paramètres changes. On se servira de la première solution comme 
%guess initial pour trouver la prochaine; si cela est avantageux comparé au
%guess initial décrit dans bvp.m

%--------------------------------------------------------------------------
%%%Variables à l'entrée:
  %-> Le bvp:
  %   [bvp.m]
  %-> Les parametres pour la résolution:
  %   [parametres_resolution.m]

%--------------------------------------------------------------------------
%%%Variables à la sortie:
  %-> si la méthode est efficace ou non:
  %   [en fprintf]

%--------------------------------------------------------------------------
 fprintf('\n****Utilisation de choix_guess_initial: vérification-->\n')
 
%--------------------------------------------------------------------------
%Renommer le guess initial normal prévu dans bvp.m:
guess_initial=y;
%Renomme la solution précédente obtenue:
K=(j-1)*N+(1:4);
solution_precedente=solutionnum(K,:);

%--------------------------------------------------------------------------
%Pour garder en mémoire le nombre d'itération de parametres_resolution.m:
nombrediterations_initial=nombrediterations;
%on veut tester qu'une seule itération:
nombrediterations=0;

%--------------------------------------------------------------------------
%Pour garder en mémoire les variables dépendentes dont on veut un
%graphique par la suite (choix_guess_initial.m ne produit pas de graphique
%des deux itérationx qu'ils calculent)
var_graphique_initial=var_graphique;
var_graphique = [];
%--------------------------------------------------------------------------
%Faire une première itération avec la solution précédente:
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
%Faire une première itération avec le guess initial prévu:
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
%Évaluer quel est le meilleur guess en comparant les erreurs après le
%calcul d'une itération:

if err_guess_initial > err_solution_precedente
    
    y=solution_precedente_1;
    fprintf('****choix_guess_initial: Utile!\n')
    
else

    y=guess_initial_1;
    fprintf('****choix_guess_initial: Pas utile...\n')

end

%--------------------------------------------------------------------------
%Redonner le bon nombre d'itération prévu dans parametres_resolution.m:
nombrediterations=nombrediterations_initial;
%Redonner les variables dépendantes dont on veut le graphique:
var_graphique=var_graphique_initial;

