function [ M, F ] = AddMatrizMEF( M, Me, nel, F, Fe )
[r,c] = size(Me);
%encontra a posição adequada para inserir a matriz Me:
id = 1+(nel-1)*(r-1);
for i = 1:r
    F(id+(i-1)) = F(id+(i-1)) + Fe(i);
    for j = 1:c
        M(id+(i-1),id+(j-1)) = M(id+(i-1),id+(j-1)) + Me(i,j);
    end
end

end

