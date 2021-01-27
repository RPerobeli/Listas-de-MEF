function [ F ] = AddVetorMEF( F, Fe, nel )
tam = length(Fe);
%encontra a posição adequada para inserir a matriz Me:
id = 1+(nel-1)*(tam-1); %tam-1 = grau
for i = 1:tam
    F(id+(i-1)) = F(id+(i-1)) + Fe(i);
end
end

