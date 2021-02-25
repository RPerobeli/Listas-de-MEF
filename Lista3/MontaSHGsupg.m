function [shg] = MontaSHGsupg(p,nint,h,beta,kappa)
%formação de uma matriz shg, que contém os valores das funções base
%baseado na quadratura gaussiana
grau = nint-1;
tau = beta*h/(2*abs(kappa));
switch grau
    case 1
        disp('shgSUPG: 1');
        for i = 1:nint
            t = p(i);
            shg(1,1,i) = 0.5*(1-t)+tau*kappa*(-0.5);
            shg(1,2,i) = 0.5*(1+t)+tau*kappa*(0.5);
            shg(2,1,i) = -0.5;
            shg(2,2,i) = 0.5;
        end
end

