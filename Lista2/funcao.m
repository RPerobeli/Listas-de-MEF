function [ y ] = funcao( x )
% y = sin(pi*x);
% y = 514/3 - log(x-1) - 511./(3.*x-3) - x.^2./6 - 2.*x./3;
% y = x.^3 -x+1;


eps = 1e-3;
c2 = (exp(-1/sqrt(eps)) - 1)/(exp(1/sqrt(eps))- exp(-1/sqrt(eps)));
c1 = -1-c2;
y = c1*exp(-x./sqrt(eps)) + c2*exp(x/sqrt(eps)) + 1;
end

