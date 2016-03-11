function [DTQ] = BuildDTQ(fw,gw,alpha,t)

global SYLVESTER_BUILD_METHOD

[r,~] = size(fw);
m = r - 1;

[r,~] = size(gw);
n = r - 1;


switch SYLVESTER_BUILD_METHOD
    case 'Standard'
        % Build Matrices D,T and Q.
        D = BuildD(m,n,t);
        T = BuildT(fw,gw,alpha,t);
        Q = BuildQ(m,n,t);
        DTQ = D*T*Q;
        
    case 'Rearranged'
        
        DTQ = BuildDTQ_ElementWise(fw,gw,alpha,t);
    otherwise
        error('Error : Build Method is either (Standard) or (Rearranged)')
end

end