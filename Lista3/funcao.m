function [ y ] = funcao( x )

eps = 1e-2;
kappa = 1;
num = 1-exp(kappa * x./eps);
den = 1-exp(kappa/eps);

y = 1/kappa * (x - num/den);
end

