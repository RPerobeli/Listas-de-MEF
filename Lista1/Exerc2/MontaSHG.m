function [ shg ] = MontaSHG( p,nint )
%formação de uma matriz shg, que contém os valores das funções base
%baseado na quadratura gaussiana
grau = nint-1;
switch grau
    case 1
        disp('shg: 1');
        for i = 1:nint
            t = p(i);
            shg(1,i) = 0.5*(1-t);
            shg(2,i) = 0.5*(1+t); 
        end
    case 2
        display('shg: 2');
        for i = 1:nint
            t = p(i);
            shg(1,i) = 0.5*t*(t-1);
            shg(2,i) = -(t-1)*(t+1);
            shg(3,i) = 0.5*t*(t+1);
        end
        
    case 3
        display('shg: 3');
        for i = 1:nint
            t = p(i);
            shg(1,i) = -9/16*(t+1/3)*(t-1/3)*(t-1);
            shg(2,i) = 27/16*(t+1)*(t-1/3)*(t-1);
            shg(3,i) = -27/16*(t+1)*(t+1/3)*(t-1);
            shg(4,i) = 9/16*(t+1)*(t+1/3)*(t-1/3);
        end
    case 4
        display('shg: 4');
        for i = 1:nint
            t = p(i);
            shg(1,i) = 2/3*(t+0.5)*t*(t-0.5)*(t-1);
            shg(2,i) = -8/3*(t+1)*t*(t-0.5)*(t-1);
            shg(3,i) = 4*(t+1)*(t+0.5)*(t-0.5)*(t-1);
            shg(4,i) = -8/3*(t+1)*(t+0.5)*t*(t-1);
            shg(5,i) = 2/3*(t+1)*(t+0.5)*t*(t-0.5);
        end
        
    otherwise
        display('ERRO: Erro na montagem de shg');
end


end

