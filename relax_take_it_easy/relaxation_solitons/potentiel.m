function [ res ] = potentiel( phi,chi,t,p)
u=ones(size(phi));
v1=(phi.^2-p(2)*u).*(phi.^2-u).^2;
v2=((chi.^2-u).^2-(p(3)/4).*(chi-2*u).*(chi+u).^2)./(phi.^2+p(1));
res=v1+v2;
end

