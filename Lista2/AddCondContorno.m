function [ M, F ] = AddCondContorno( M ,F, ka, kb, ga, gb, qa, qb )
%Fun��o que adiciona de forma fraca a condi��o composta de Robin no
%problema de elementos finitos. Mexendo nas matrizes M e F globais

[m,n] = size(M);
M(1,1) = M(1,1) + ka;
M(m,n) = M(m,n) + kb;

F(1) = F(1) + (ka*ga - qa);
F(m) = F(m) + (kb*gb - qb);


end

