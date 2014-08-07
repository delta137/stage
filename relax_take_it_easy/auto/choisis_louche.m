voir=open('res_eric.mat')
iter=0;
for gamma=1:7
    for d1=1:7
        for d2=1:7
            if voir.liste(gamma,d1,d2) == 1
               iter=iter+1;
               %lista(iter)=[gamma d1 d2];
            end
        end
    end
end
