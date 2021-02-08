function [ dy ] = dfuncao( x )
% dy = pi*cos(pi*x);
% dy = - (x.^3 - 512)./(3*(x-1).^2);
% dy = 3*x.^2 - 1;
 
eps = 1e-3;
c2 = (exp(-1/sqrt(eps)) - 1)/(exp(1/sqrt(eps))- exp(-1/sqrt(eps)));
c1 = -1-c2;
dy = c1*(-1./sqrt(eps))*exp(-x./sqrt(eps)) + c2*(1./sqrt(eps))*exp(x/sqrt(eps));
end

