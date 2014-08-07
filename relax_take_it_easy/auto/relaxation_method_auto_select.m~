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
   iter1=0;
  iter2=0;
 %[gamma,delta1,delta2]=choisis;
%ii=1;
gamma=[0.001 0.01 0.1 0.2 0.4 0.6]
for ii=1:6  
for gamma=1:1
    for d1=1:1
        for d2=1:1
if(mod(iter1,2)==0)
    d1eff=d1;
else
    d1eff=8-d1;
end
if(mod(iter2,2)==0)
    d2eff=d2;
else
    d2eff=8-d2;
end

  parametres=[0.1*gamma,0.1*d1eff,0.1*d2eff];    
%Param�tres pour la r�solution du probl�me:
parametres_resolution_bounceO4_6;

%Le bvp:
[g,y,bc] = bvp_bounceO4_6(parametres,Y,X,x);

%Relaxation method:
if(iter2>0)
    choix_guess_initial;
end
relaxation_method;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcul des densit�s d'�nergie:
%  rhograd = solutionnum((j-1)*N+2,:).^2 + (parametres(3)./x).^2 .*( 1-solutionnum((j-1)*N+3,:) ).^2 .*solutionnum((j-1)*N+1,:).^2;
%  rhomag = (1/2)*(parametres(3)*solutionnum((j-1)*N+4,:)./( parametres(1)*x ) ).^2;
%  rhopot  = ( solutionnum((j-1)*N+1,:).^2 -parametres(2) ).*( solutionnum((j-1)*N+1,:).^2 -1 ).^2;
%  rhotot  = rhograd + rhomag + rhopot;

%Calcul de l'�nergie totale:   
%  Etot =2*pi* h*trapz(x(2:M).*rhotot(2:M));
%  fprintf('\n�nergie totale: %d \n',Etot)
   
%Graphique de la densit� d'�nergie en fonction de r
%  figure(3)
%  plot(x,rhotot,'k-', x,rhograd,'r:', x,rhomag,'b:', x,rhopot,'g:')
%     title('Densit� d''�nergie en fonction de r');
%     xlabel('r');
%     ylabel('Densit� d''�nergie');
%     legend('\rho_{tot}','\rho_{grad}','\rho_{mag}','\rho_{pot}');
%     
%Sact(ii)=2*action(y(1,:),y(3,:),y(2,:),y(4,:),x,parametres);
if err(1,j)>conv 
    liste(gamma,d1eff,d2eff)=1;
%iter=iter+1;    
else
    liste(gamma,d1eff,d2eff)=2;
end
phi(gamma,d1eff,d2eff,:)=y(1,:);
chi(gamma,d1eff,d2eff,:)=y(3,:);
        end 
  iter2=iter2+1;  end
iter1=iter1+1; end
save('res.mat','phi','chi','x','liste');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

