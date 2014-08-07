

ou bien:




%%%%%%%%%                 SHOOTING METHOD                 %%%%%%%%%%%%%%%%%
%%%%%%%%% Méthode Canon + ode45 pour résoudre le problème %%%%%%%%%%%%%%%%%
%%%%%%%%%        aux valeurs frontières bvp.m             %%%%%%%%%%%%%%%%%
%
%Cette fonction solutionne un système de 2 équations différentielles à
%conditions frontière (bvp.m), y(1)(x) et y(2)(x) avec la méthode «canon».
%
%La solution numérique devra respecter les conditions frontières suivantes:
%  y(1)(x_min) = y1_initial    y(1)(x_max) = y1_final
%  y(2)(x_min) = y2_initial    y(2)(x_max) = y2_final
%
%Le problème à 2 équations diff donne 2 conditions aux frontières:
%  y(2)(x_min) = y2_initial    y(1)(x_max) = y1_final
%
%La méthode canon consiste à transformer le bvp en un problème aux valeurs
%initiales qui va servir de base pour le shooting method:
%
%On transforme les condtions frontières en condition initiale:
%  y(1)(x_min) = y1_initial
%  y(2)(x_min) = y2_initial
%
%La question est: quel est y1_initial?
%
%Il faudra choisir y1_initial, puis intégrer l'équation différentielle (ici
%nous prendrons ode45). Si la condition frontière y(1)(x_max) = y1_final
%n'est pas respectée, on modifie y1_initial jusqu'à ce qu'elle soit satisfaite.

%--------------------------------------------------------------------------
%%%Variables à l'entrée:

  %-> Paramètres spécifique aux bvp du problème:
  %   [parametres]
  %-> L'erreur (y1_final_num-y1_final) qui est demandée:
  %   [conv]
  %-> nombre d'itérations maximum:
  %   [nombrediterations]
  %-> Le range en x:
  %   [x_min   x_max] = Range_x
  %-> Les conditions frontières qui doit être respectée:
  %   [y2_initial   y1_final] = Bound_Cond
  %-> Le range de y1_initial qui sera testé:
  %   [y1_initial_min   y1_initial_max] = Range_y1_initial
 
%--------------------------------------------------------------------------
%%%Variables à la sortie:
 
  %-> Est-ce que la méthode a fonctionnée? (nombre d'itérations)
  %   (en fprintf)
  %-> La condition initiale manquante:
  %   [y1_initial] (en fprintf)
  %-> L'erreur de y(1) à la frontière:
  %   [y1_final_num-y1_final] = [err] (en fprintf)
  %-> La solution numérique:
  %   [x,y]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x y ] = shooting_method( parametres,conv,nombrediterations,...
                                             Range_x,Bound_Cond,Range_y1_initial )
%où:

x_min = Range_x(1);
x_max = Range_x(2);

y2_initial = Bound_Cond(1);
y1_final = Bound_Cond(2);

y1_initial_min = Range_y1_initial(1);
y1_initial_max = Range_y1_initial(2);    

for i=1:1:nombrediterations
   
%Première intégration avec ode45:------------------------------------------

options = odeset('RelTol',conv);
[x,y] = ode45(@(x,y)bvp(x,y,parametres),...
              [x_min x_max],...
              [y1_initial_min   y2_initial],...
              options ) ;
%où:

  %-> bvp est le système d'équations différentielles à résoudre avec des
  %   paramètres prédéfinis.
  %-> La variable indépendante x passe de x_min à x_max
  %-> y1_initial_min est la valeur initiale de y(1) à x=0
  %-> y2_initial     est la valeur initiale de y(2) à x=0
  %-> L'option permet de choisir la tolerance relative voulue, c'est-à-dire
  %   «conv», prédéfinie.
 
%Valeur de y(1) à la frontière:
y1_final_num=y(length(x),1);
%L'erreur est donc de:
err1=y1_final_num-y1_final;

if abs(err1)<conv
   err=err1;
   y1_initial=y1_initial_min;
   break
end

%Deuxière intégration avec ode45:------------------------------------------

[x,y] = ode45(@(x,y)bvp(x,y,parametres),...
              [x_min x_max],...
              [y1_initial_min+(y1_initial_max-y1_initial_min)/2   y2_initial],...
              options ) ;
         
%Valeur de y(1) à la frontière:
y1_final_num=y(length(x),1);
%L'erreur est donc de:
err2=y1_final_num-y1_final;

if abs(err2)<conv
   err=err2;
   y1_initial=y1_initial_min+(y1_initial_max-y1_initial_min)/2;
   break
end

%--------------------------------------------------------------------------

%UNDERSHOOT:
if err1 * err2 > 0
%Si err1 et err2 SONT du même signe (alors leur produit est POSITIF)...
%Redéfinir y1_initial_min et y1_initial_max:
y1_initial_min = y1_initial_min +(y1_initial_max - y1_initial_min)/2;
y1_initial_max = y1_initial_max;

%OVERSHOOT:
else
%Si err1 et err2 NE sont PAS du même signe (alors leur produit est NÉGATIF)...
%Redéfinir y1_initial_min et y1_initial_max:
y1_initial_min = y1_initial_min;
y1_initial_max = y1_initial_min +(y1_initial_max - y1_initial_min)/2;

end


end %Fin des itérations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%La méthode canon a-t-elle fonctionnée?
if i<nombrediterations
   fprintf('La méthode canon a fonctionnée à l''itération # %d !!! \n',i)
   fprintf('La condition initiale manquante trouvée numériquement: %d \n',y1_initial)
   fprintf('L''erreur de y(1) à la frontière: %d \n',err)
else
   fprintf('La méthode canon n''a pas fonctionnée pour ce nombre maximum d''itération (%d).\n',i)
   fprintf('La condition initiale manquante trouvée numériquement: %d \n',y1_initial_max)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%La solution trouvée est: [x,y];

