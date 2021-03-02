clear all
close all
clc

%% Definicoes previas

eps = 1e-2;
kappa = 1;
gamma = 0;

%Definicao do dominio Omega
xi = 0.0;
xf = 2.0;

%% inicio do loop para diferentes refinamentos de malha
vetErro = [];

Pe_h = [1,5,10];
pe = 1;
h = Pe_h(pe)*2*eps/abs(kappa);

dt = h^2;
t = 0:dt:1.25;

beta = coth(Pe_h(pe))-1/Pe_h(pe);
tau = beta*h/(2*abs(kappa));
%Definicao malha
nel = (xf-xi)/h;

x = xi:h:xf;
x_ex = xi:0.001:xf;

%funcao exata desejada:
exata = funcao(x,1.25);

%sao necessarios (grau + 1) nos no elemento para representar as funcoes
%shg.
grau = 1;
nen = grau + 1;

%numero de nos global: grau*nel +1
numNos_g = grau*nel + 1;
h_el = h/grau;
x_global = xi:h_el:xf;

%% Set condicao inicial
alfa = funcao(x_global,0);

for k = 1:length(t)
    %% condicoes de contorno (gerais) para formulacao fraca
    kap_a = 1e10;
    kap_b = 1e10;
    %Dirichlet
    g_a = funcao(xi, t(k));
    g_b = funcao(xf, t(k));

    %Neumman
    q_a = 0;
    q_b = 0;
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

    %% Solucao via SUPG
    M = zeros(numNos_g,numNos_g);
    K = zeros(numNos_g,numNos_g);
    C = zeros(numNos_g,numNos_g);
    F = zeros(numNos_g,1);
    [w, p] = MontaQuadraturaGaussiana(grau);
    shg = MontaSHG(p,nint,h_el);
    shgSUPG = MontaSHGsupg(p,nint,h_el,beta,kappa);

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
                Fe(j) = Fe(j) + fonte(x_ref)*(shg(1,j,l)+tau*kappa*shg(2,j,l)*1/dx)*w(l)*dx;
                for i = 1:nen
                    Me(i,j) = Me(i,j)+ gamma*shg(1,i,l)*(shg(1,j,l)+tau*kappa*shg(2,j,l)*1/dx)*w(l)*dx; 
                    Ke(i,j) = Ke(i,j)+ eps*shg(2,i,l)*1/dx*shg(2,j,l)*1/dx*w(l)*dx;
                    Ce(i,j) = Ce(i,j)+ kappa*shg(2,i,l)*1/dx*(shg(1,j,l)+tau*kappa*shg(2,j,l)*1/dx)*w(l)*dx;
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
    solSUPG = M_final\F;
end

%% Plot funcao e aproximacao
    figure;
    plot(x_global,alfa,'b', x_ex,funcao(x_ex));
    nome = num2str(nel);
    frase = 'Funcao aproximada com malha ';
    title(strcat(frase,nome));
    % title("Funcao aproximada com malha "+ nel);
    xlabel('x');
    ylabel('f(x)');