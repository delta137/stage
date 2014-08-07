function [res ] = action(phi,chi,phip,chip,t,p)
integrand=2*pi*t.*((1/2)*(phip.^2+chip.^2) + potentiel(phi,chi,t,p))
res=trapz(t,integrand)
end

