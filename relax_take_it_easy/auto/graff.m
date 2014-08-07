function [ ] = graff()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%La fonction squeeze les matrices 4D,
%je choisis des valeurs de g,d1,d2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=open('res_eric.mat')
for gamma=1:7
    for d1=1:7
        for d2=1:7
            hold;
            plot(res.x, squeeze(res.phi(gamma,d1,d2,:)),res.x,squeeze(res.chi(gamma,d1,d2,:)))
            pause(3);
            close;
            gamma
            d1
            d2
            pause(1);
            
        end
    end
end

end

