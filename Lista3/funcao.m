function [ y ] = funcao( x )

% eps = 1e-2;
% kappa = 1;
% num = 1-exp(kappa * x./eps);
% den = 1-exp(kappa/eps);
% 
% y = 1/kappa * (x - num/den);

y = -1.90476*exp(-5.*x)+3.9604*exp(-x)+ 5.99498*10^-44*exp(100.*x) - 2.05563;
end

