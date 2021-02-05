function [ w,pontos ] = MontaQuadraturaGaussiana( grau )
%Baseado no grau do polinomio, retorna um vetor com os valores dos pesos e
%um vetor com pontos a serem usados na fun��o F da quadratura
switch grau
    case 1
        w = [1;1];
        pontos = [-sqrt(3.)/3.; sqrt(3.)/3.];
    case 2
        w = [5./9.; 8./9.; 5./9.];
        pontos = [-sqrt(3./5.); 0.0; sqrt(3./5.)];
    case 3
         w = [
             (18.-sqrt(30.))/36.;
             (18.+sqrt(30.))/36.;
             (18.+sqrt(30.))/36.;
             (18.-sqrt(30.))/36.];
        pontos = [
            -sqrt(3./7. +(2./7.)*sqrt(6./5.));
            -sqrt(3./7. -(2./7.)*sqrt(6./5.));
            sqrt(3./7. -(2./7.)*sqrt(6./5.));
            sqrt(3./7. +(2./7.)*sqrt(6./5.))];
    case 4
        %w = [0.23692689; 0.478628;0.5688889;0.236926; 0.478628];
        %pontos = [-0.90617985;-0.53846931;0;0.90617985;0.53846931];
        w = [(322.-13.*sqrt(70.))/900;
            (322.+13.*sqrt(70.))/900;
            128/225;
            (322.+13.*sqrt(70.))/900;
            (322.-13.*sqrt(70.))/900];
        pontos = [-1./3.*sqrt(5.+2.*sqrt(10./7.));
            -1./3.*sqrt(5.-2.*sqrt(10./7.));
            0;
            1./3.*sqrt(5.-2.*sqrt(10./7.));
            1./3.*sqrt(5.+2.*sqrt(10./7.));];
    otherwise
        display('ERRO: erro na montagem da quadratura')
end


end

