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
for k= 2:6
%Defini��o malha
nel = 2^k;
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
F = zeros(numNos_g,1);
[w, p] = MontaQuadraturaGaussiana(grau);
shg = MontaSHG(p,nint);

for n = 1:nel
    %matriz do elemento e vetor fonte do elemento
    Me = zeros(nen,nen);
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
% erroL2 = 0;
% for n = 1:nel
%     erro = 0;
%     %definir as posi��es de cada elemento no frame global
%     xi_el = xi+(n-1)*h;
%     xf_el = xi+n*h;
%     h_sec = (xf_el - xi_el)/grau;
%     xl = xi_el:h_sec:xf_el; %vetor com as posi��es dos n�s do elemento no frame global e pontos secund�rios
%     for l = 1:nint
%         uh = 0;
%         x_ref=0;
%         for i = 1:nen
%             %Os alfas devem ser pegos nas posi�oes globais, n+i-1
%             uh = uh + shg(i,l)*alfa(i+(n-1)*grau); 
%             x_ref = x_ref + shg(i,l)*xl(i);
%         end
%         erro = erro + ((funcao(x_ref)- uh)^2)*w(l)*dx;
%     end
%     erroL2 = erroL2 +erro;    
% end
% erroL2 = sqrt(erroL2);

% erro na norma do m�ximo
for n = i:numNos_g
    erro(n) = abs(alfa(n) - funcao(x_global(n)));
end
erroL2 = max(erro);
passos(k-1) = h;
vetErro(k-1) = erroL2;
end

%% Plot fun��o e aproximacao
figure;
plot(x_global,alfa, x_ex,funcao(x_ex));

%% Plot Convergencia
figure(2);
plot(-log10(passos),log10(vetErro), '*-b');