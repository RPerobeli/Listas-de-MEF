clear all;
close all;
clc;

%{
Resolver o sistema M*alfa = b, onde alfa sao parametros em que a combina��o
linear sum(alfa*phi) gera a fun��o exata que se quer aproximar.
%}

%% Defini��es pr�vias

%Defini��o do dom�nio �mega
xi = -2.0;
xf = 2.0;

vetErro = [];
for k= 2:5
%Defini��o malha
nel = 4^k;
% nel = 4;
h = (xf-xi)/nel; %tamanho do elemento
x = xi:h:xf;
x_ex = xi:0.001:xf;

%fun��o exata desejada:
exata = funcao(x);

%Defini��o do grau do polinomio a ser modelado, 1 - modelado por uma reta,
%2- modelado por uma par�bola
%s�o necess�rios (grau + 1) n�s no elemento para representar as fun��es
%shg.
grau = 1;
nen = grau + 1;

%n�mero de n�s global: grau*nel +1
numNos_g = grau*nel + 1;

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
alfaFinal=[];
[w, p] = MontaQuadraturaGaussiana(grau);
shg = MontaSHG(p,nint);
x_global = [];
for n = 1:nel
    %matriz do elemento e vetor fonte do elemento
    Me = zeros(nen,nen);
    Fe = zeros(nen,1);
    %definir as posi��es de cada elemento no frame global
    xi_el = xi+(n-1)*h;
    xf_el = xi+n*h;
    h_sec = (xf_el - xi_el)/grau;
    xl = xi_el:h_sec:xf_el; %vetor com as posi��es dos n�s do elemento no frame global e pontos secund�rios
    x_global = [x_global, xl];
    for l = 1:nint
        x_ref = 0;
        for i = 1:nen
            %necess�rio avaliar a fun��o u no ponto de referencia, para o
            %vetor fonte, x_ref = t(x) transforma��o isoparametrica, se n�o se
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
    %Resolve o sistema linear
    alfa = Me\Fe;
    alfaFinal = [alfaFinal; alfa];
end



%% Calculo do erro
erroL2 = 0;

for n = i:length(x_global)
    erro(n) = abs(alfaFinal(n) - funcao(x_global(n)));
end
erroL2 = max(erro);
passos(k-1) = h;
vetErro(k-1) = erroL2;

%% Plot fun��o e aproximacao
figure;
plot(x_global,alfaFinal, x_ex,funcao(x_ex));
title("Fun��o aproximada com malha "+nel);
xlabel("x");
ylabel("f(x)")
end

%% Plot Convergencia
figure;
plot(-log10(passos),log10(vetErro), '*-b');