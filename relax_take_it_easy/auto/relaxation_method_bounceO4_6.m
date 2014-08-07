clear;
clc;
close all %(close all figure)
format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auteur: Marie-Lou Gendron Marsolais
% Date: 25 juin 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programme pour solutionner f(r) et a(r) par la �relaxation method�

%Valeurs des param�tres du bvp.m:
  %-> parametres=[gamma,delta1,delta2]
iter=0;
[gamma,d1,d2]=choisis()
for ii=1:1
  %parametres=[0.4,0.7-ii/10,0.6];
  parametres=[gamma(ii),d1(ii),d2(ii)]
%Param�tres pour la r�solution du probl�me:
parametres_resolution_bounceO4_6

%Le bvp:
[g,y,bc] = bvp_bounceO4_6(parametres,Y,X,x);

%Relaxation method:
if(iter>0)
    choix_guess_initial
end
relaxation_method

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcul des densit�s d'�nergie:
%  rhograd = solutionnum((j-1)*N+2,:).^2 + (parametres(3)./x).^2 .*( 1-solutionnum((j-1)*N+3,:) ).^2 .*solutionnum((j-1)*N+1,:).^2;
%  rhomag = (1/2)*(parametres(3)*solutionnum((j-1)*N+4,:)./( parametres(1)*x ) ).^2;
%  rhopot  = ( solutionnum((j-1)*N+1,:).^2 -parametres(2) ).*( solutionnum((j-1)*N+1,:).^2 -1 ).^2;
%  rhotot  = rhograd + rhomag + rhopot;
% 
% %Calcul de l'�nergie totale:   
%  Etot =2*pi* h*trapz(x(2:M).*rhotot(2:M));
%  fprintf('\n�nergie totale: %d \n',Etot)
%    
% %Graphique de la densit� d'�nergie en fonction de r
%  figure(3)
%  plot(x,rhotot,'k-', x,rhograd,'r:', x,rhomag,'b:', x,rhopot,'g:')
%     title('Densit� d''�nergie en fonction de r');
%     xlabel('r');
%     ylabel('Densit� d''�nergie');
%     legend('\rho_{tot}','\rho_{grad}','\rho_{mag}','\rho_{pot}');
%     
Sact=2*action(y(1,:),y(3,:),y(2,:),y(4,:),x,parametres)
iter=iter+1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

