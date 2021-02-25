clear all
close all
clc

%% Definicoes previas

eps = 1e-2;
kappa = 1;
gamma = 0;

%Definicao do dominio Omega
xi = 0.0;
xf = 1.0;

%% condicoes de contorno (gerais) para formulacao fraca
kap_a = 1e10;
kap_b = 1e10;
%Dirichlet
g_a = 0;
g_b = 1;
%Neumman
q_a = 0;
q_b = 0;

%% inicio do loop para diferentes refinamentos de malha
vetErro = [];
%Definição malha
nel = 10;
h = (xf-xi)/nel;

Pe_h = kappa*h/(2*eps);
% beta = coth(Pe_h(i))-1/Pe_h(i);
beta = [0.2850];
erroBeta=[];
for(k = 1:length(beta))

x = xi:h:xf;
x_ex = xi:0.001:xf;

%funcao exata desejada:
exata = funcao(x);

%sao necessarios (grau + 1) nos no elemento para representar as funcoes
%shg.
grau = 1;
nen = grau + 1;

%numero de nos global: grau*nel +1
numNos_g = grau*nel + 1;
h_el = h/grau;
x_global = xi:h_el:xf;
%% Transformacoes para elemento de referencia

%Para a quadratura de Gauss, e possivel determinar  integrais exatas de 
%polinomios de ate grau 2n+1 usando n como grau do polinomio da base phi.
nint = nen; %nint e o numero de pontos necesserios para realizar integral exata por QG
dx = h/2;
%transformacao isoparametrica:
%t = (2*x - xi - xf)/h;
%x_el = dx*t + (xi_el+xf_el)/2
%sao importantes para passar a funcao shg, ou phi, para o elemento de
%referencia

%% Loop de montagem da matriz Global via Bolha
M = zeros(numNos_g,numNos_g);
K = zeros(numNos_g,numNos_g);
C = zeros(numNos_g,numNos_g);
F = zeros(numNos_g,1);
[w, p] = MontaQuadraturaGaussiana(grau);
shg = MontaSHG(p,nint,h_el);
shgBolha = MontaSHGBolha(p,nint,h_el,beta(k),kappa);

for n = 1:nel
    %matriz do elemento e vetor fonte do elemento
    Me = zeros(nen,nen);
    Ke = zeros(nen,nen);
    Ce = zeros(nen,nen);
    Fe = zeros(nen,1);
    %definir as posicoes de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posicoes dos nos do elemento no frame global e pontos secundarios
    
    for l = 1:nint
        x_ref = 0;
        for i = 1:nen
            %necessario avaliar a funcao u no ponto de referencia, para o
            %vetor fonte, x_ref = t(x) transformacao isoparametrica, se nao se
            %sabe: faz-se como o mostrado abaixo
            x_ref = x_ref + shg(1,i,l)*xl(i);
        end
        for j = 1:nen
            Fe(j) = Fe(j) + fonte(x_ref)*shg(1,j,l)*w(l)*dx;
            for i = 1:nen
                Me(i,j) = Me(i,j)+ gamma*shg(1,i,l)*shgBolha(1,j,l)*w(l)*dx; 
                Ke(i,j) = Ke(i,j)+ (eps)*shg(2,i,l)*1/dx*shgBolha(2,j,l)*1/dx*w(l)*dx;
                Ce(i,j) = Ce(i,j)+ kappa*shg(2,i,l)*1/dx*shgBolha(1,j,l)*w(l)*dx;
            end
        end
    end
    [M,F] = AddMatrizMEF(M, Me, n, F, Fe);
    K = AddMatrizSemFonte(K, Ke, n);
    C = AddMatrizSemFonte(C, Ce', n);
end

%% Resolve o sistema linear
M_final = K+M+C;
[M_final,F] = AddCondContorno(M_final, F, kap_a, kap_b, g_a, g_b, q_a, q_b);
alfa = M_final\F;

%% Plot funcao e aproximacao
figure;
plot(x_global,alfa,'b', x_ex,funcao(x_ex), 'r');
nome = num2str(beta);
frase = 'Funcao aproximada por Bolha com beta ';
title(strcat(frase,nome));
% title("Funcao aproximada com malha "+ nel);
xlabel('x');
ylabel('f(x)');
end
%% Calcula o erro com a função exata:
%{
erroL2 = 0;

for n = 1:nel
    erro = 0;
    %definir as posições de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posições dos nós do elemento no frame global e pontos secundários
    for l = 1:nint
        uh = 0;
        x_ref=0;
        for i =1:nen
            uh = uh + shg(1,i,l)*alfa(i+(n-1)*grau);
            x_ref = x_ref + shg(1,i,l)*xl(i);
        end
        erro = erro + ((funcao(x_ref)- uh)^2)*w(l)*dx;
    end
    erroL2 = erroL2 +erro;
end
erroL2 = sqrt(erroL2);
erroBeta = [erroBeta, erroL2];
end
%% plot do erro
figure;
plot(beta, erroBeta);
title("Erro(L2) em funcao de Beta");
xlabel('\beta');
ylabel('erro');


%% Procura beta otimo
minErro = min(erroBeta);
betaOtimo = -1;
for g=1:length(erroBeta)
    if(erroBeta(g) == minErro)
        betaOtimo = beta(g);
        break;
    end
end
disp(betaOtimo);
%}
