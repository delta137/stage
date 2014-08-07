clear;
C=load('res_meilleur2.txt')
iter=1;
vec=abs(C(:,4,1));
a=size(vec)
for ii= 1:a(1)
    if vec(ii,1)<0.001
        choix(iter)=ii
        iter=iter+1
    end
end
gamma=C(choix,1,1)
delta1=C(choix,2,1)
delta2=C(choix,3,1)