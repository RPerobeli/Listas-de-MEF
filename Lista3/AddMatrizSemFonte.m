function [M] = AddMatrizSemFonte(M, Me, nel)
[r,c] = size(Me);
%encontra a posição adequada para inserir a matriz Me:
id = 1+(nel-1)*(r-1);
for i = 1:r
    for j = 1:c
        M(id+(i-1),id+(j-1)) = M(id+(i-1),id+(j-1)) + Me(i,j);
    end
end

end

