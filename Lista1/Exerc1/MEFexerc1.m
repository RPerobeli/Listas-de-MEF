clear all;
close all;
clc;

%{
Resolver o sistema M*alfa = b, onde alfa sao parametros em que a combinação
linear sum(alfa*phi) gera a função exata que se quer aproximar.
%}

%% Definições prévias

%Definição do domínio ômega
xi = -2.0;
xf = 2.0;

vetErro = [];
for k= 2:5
%Definição malha
nel = 4^k;
% nel = 4;
h = (xf-xi)/nel; %tamanho do elemento
x = xi:h:xf;
x_ex = xi:0.001:xf;

%função exata desejada:
exata = funcao(x);

%Definição do grau do polinomio a ser modelado, 1 - modelado por uma reta,
%2- modelado por uma parábola
%são necessários (grau + 1) nós no elemento para representar as funções
%shg.
grau = 4;
nen = grau + 1;

%número de nós global: grau*nel +1
numNos_g = grau*nel + 1;
h_el = h/grau;
x_global = xi:h_el:xf;
%% Transformaçoes para elemento de referencia

%Para a quadratura de Gauss, é possível determinar  integrais exatas de 
%polinomios de até grau 2n+1 usando n como grau do polinomio da base phi.
nint = nen; %nint é o numero de pontos necessários para realizar integral exata por QG
dx = h/2;
%transformação isoparamétrica:
%t = (2*x - xi - xf)/h;
%x_el = dx*t + (xi_el+xf_el)/2
%são importantes para passar a função shg, ou phi, para o elemento de
%referencia

%% Loop de montagem da matriz Global
M = zeros(numNos_g,numNos_g);
F = zeros(numNos_g,1);
[w, p] = MontaQuadraturaGaussiana(grau);
shg = MontaSHG(p,nint);

for n = 1:nel
    %matriz do elemento e vetor fonte do elemento
    Me = zeros(nen,nen);
    Fe = zeros(nen,1);
    %definir as posições de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posições dos nós do elemento no frame global e pontos secundários
    
    for l = 1:nint
        x_ref = 0;
        for i = 1:nen
            %necessário avaliar a função u no ponto de referencia, para o
            %vetor fonte, x_ref = t(x) transformação isoparametrica, se não se
            %sabe: faz-se como o mostrado abaixo
            x_ref = x_ref + shg(i,l)*xl(i);
        end
        for j = 1:nen
            Fe(j) = Fe(j) + funcao(x_ref) * shg(j,l)*w(l)*dx;
            for i = 1:nen
                Me(i,j) = Me(i,j)+ shg(i,l)*shg(j,l)*w(l)*dx;   
            end
        end
    end
    [M,F] = AddMatrizMEF(M, Me, n, F, Fe);
end

%% Resolve o sistema linear
alfa = M\F;

%% Calculo do erro
erroL2 = 0;

for n = i:numNos_g
    erro(n) = abs(alfa(n) - funcao(x_global(n)));
end
%% Plot função e aproximacao
figure;
plot(x_global,alfa, x_ex,funcao(x_ex));
title("Função aproximada com malha "+nel);
xlabel("x");
ylabel("f(x)")

%% Preparacao para o erro
erroL2 = max(erro);
passos(k-1) = h;
vetErro(k-1) = erroL2;
end



%% Plot Convergencia
figure;
plot(-log10(passos),log10(vetErro),'*-b');
xlabel("-log10(h)");
ylabel("log10(erro)")