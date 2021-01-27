function [ w,pontos ] = MontaQuadraturaGaussiana( grau )
%Baseado no grau do polinomio, retorna um vetor com os valores dos pesos e
%um vetor com pontos a serem usados na função F da quadratura
switch grau
    case 1
        w = [1;1];
        pontos = [-sqrt(3)/3; sqrt(3)/3];
    case 2
        w = [0.5555556; 0.888889; 0.5555556];
        pontos = [-0.7746; 0.0; 0.7746];
    case 3
        w = [0.3478;0.65214515;0.3478;0.65214515];
        pontos = [-0.86113631;-0.33998;0.86113631;0.33998];
    case 4
        w = [0.23692689; 0.478628;0.5688889;0.236926; 0.478628];
        pontos = [-0.90617985;-0.53846931;0;0.90617985;0.53846931];
    otherwise
        display('ERRO: erro na montagem da quadratura')
end


end

