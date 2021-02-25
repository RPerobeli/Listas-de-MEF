function [shg] = MontaSHGBolha(p,nint,h,beta,kappa)
%formação de uma matriz shg, que contém os valores das funções base
%baseado na quadratura gaussiana
grau = nint-1;
switch grau
    case 1
        disp('shgBolha: 1');
        for i = 1:nint
            t = p(i);
            sum = 0.75*beta*(1-t^2);
            shg(1,1,i) = 0.5*(1-t)-sign(kappa)*sum;
            shg(1,2,i) = 0.5*(1+t)+sign(kappa)*sum;
            shg(2,1,i) = -0.5+3/2*beta*t;
            shg(2,2,i) = 0.5-3/2*beta*t;
        end
end

