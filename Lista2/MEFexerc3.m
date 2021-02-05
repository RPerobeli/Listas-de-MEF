clear all
close all
clc

%% Defini��es pr�vias

%Defini��o do dom�nio �mega
xi = 0.0;
xf = 1.0;

%% condicoes de contorno (gerais) para formulacao fraca
kap_a = 1e10;
kap_b = 1e10;
%Dirichlet
g_a = 0;
g_b = 0;
%Neumman
q_a = 0;
q_b = 0;

%% inicio do loop para diferentes refinamentos de malha
vetErro = [];
for g = 1:4
for k= 2:5
%Defini��o malha
nel = 2^k;
% nel = 42*(k-1);
h = (xf-xi)/nel; %tamanho do elemento
% h = 0.07; %tamanho fixo do elemento para alguns testes
x = xi:h:xf;
x_ex = xi:0.001:xf;

%fun��o exata desejada:
exata = funcao(x);

%Defini��o do grau do polinomio a ser modelado, 1 - modelado por uma reta,
%2- modelado por uma par�bola
%s�o necess�rios (grau + 1) n�s no elemento para representar as fun��es
%shg.
grau = g;
nen = grau + 1;

%n�mero de n�s global: grau*nel +1
numNos_g = grau*nel + 1;
h_el = h/grau;
x_global = xi:h_el:xf;
%% Transforma�oes para elemento de referencia

%Para a quadratura de Gauss, � poss�vel determinar  integrais exatas de 
%polinomios de at� grau 2n+1 usando n como grau do polinomio da base phi.
nint = nen; %nint � o numero de pontos necess�rios para realizar integral exata por QG
dx = h/2;
%transforma��o isoparam�trica:
%t = (2*x - xi - xf)/h;
%x_el = dx*t + (xi_el+xf_el)/2
%s�o importantes para passar a fun��o shg, ou phi, para o elemento de
%referencia

%% Loop de montagem da matriz Global
M = zeros(numNos_g,numNos_g);
K = zeros(numNos_g,numNos_g);
F = zeros(numNos_g,1);
[w, p] = MontaQuadraturaGaussiana(grau);
shg = MontaSHG(p,nint);

for n = 1:nel
    %matriz do elemento e vetor fonte do elemento
    Me = zeros(nen,nen);
    Ke = zeros(nen,nen);
    Fe = zeros(nen,1);
    %definir as posi��es de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posi��es dos n�s do elemento no frame global e pontos secund�rios
    
    for l = 1:nint
        x_ref = 0;
        for i = 1:nen
            %necess�rio avaliar a fun��o u no ponto de referencia, para o
            %vetor fonte, x_ref = t(x) transforma��o isoparametrica, se n�o se
            %sabe: faz-se como o mostrado abaixo
            x_ref = x_ref + shg(1,i,l)*xl(i);
        end
        for j = 1:nen
            Fe(j) = Fe(j) + fonte(x_ref)*shg(1,j,l)*w(l)*dx;
            for i = 1:nen
                Me(i,j) = Me(i,j)+ funcaoGamma(x_ref)*shg(1,i,l)*shg(1,j,l)*w(l)*dx; 
                Ke(i,j) = Ke(i,j)+ funcaoK(x_ref)*shg(2,i,l)*1/dx*shg(2,j,l)*1/dx*w(l)*dx;   
            end
        end
    end
    [M,F] = AddMatrizMEF(M, Me, n, F, Fe);
    K = AddMatrizSemFonte(K, Ke, n);
end

%% Resolve o sistema linear
M_final = K+M;
[M_final,F] = AddCondContorno(M_final, F, kap_a, kap_b, g_a, g_b, q_a, q_b);
alfa = M_final\F;

%% Calculo do erro
erroL2 = 0;

for n = 1:nel
    erro = 0;
    %definir as posi��es de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posi��es dos n�s do elemento no frame global e pontos secund�rios
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

%% Erro da derivada da fun��o

erroL2d = 0;

for n = 1:nel
    errod = 0;
    %definir as posi��es de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posi��es dos n�s do elemento no frame global e pontos secund�rios
    for l = 1:nint
        duh = 0;
        x_ref=0;
        for i =1:nen
            duh = duh + shg(2,i,l)*1/dx*alfa(i+(n-1)*grau);
            x_ref = x_ref + shg(1,i,l)*xl(i);
        end
        errod = errod + ((dfuncao(x_ref)- duh)^2)*w(l)*dx;
    end
    erroL2d = erroL2d +errod;
end
erroL2d = sqrt(erroL2d);
%% Plot fun��o e aproximacao
% figure;
% plot(x_global,alfa,'b', x_ex,funcao(x_ex), 'r' );
% nome = num2str(nel);
% frase = 'Fun��o aproximada com malha ';
% title(strcat(frase,nome));
% % title("Fun��o aproximada com malha "+ nel);
% xlabel('x');
% ylabel('f(x)')

%% Preparacao para o erro
passos(k-1) = h;
vetErro(k-1) = erroL2;
vetErrod(k-1) = erroL2d;
end

%% Plot Convergencia
figure(1);
plot(-log10(passos),log10(vetErro), '-*');
hold on
title('Erro da fun��o');
xlabel('-log10(h)');
ylabel('log10(erro)')
figure(2);
plot(-log10(passos),log10(vetErrod),'-*');
hold on
title('Erro da derivada da fun��o');
xlabel('-log10(h)');
ylabel('log10(erro_{derivada})')
end