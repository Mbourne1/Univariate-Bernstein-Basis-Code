function [DTQ] = BuildDTQ(fw,gw,alpha,t)
% Build the matrix DTQ 
%
% fw : Coefficients of polynomial f(\omega,\theta)
%
% gw : Coefficients of polynomial g(\omega,\theta)
%
% alpha : 
%
% t : Degree of GCD d(x).

% Global Variables.
global SYLVESTER_BUILD_METHOD

% Get degree of polynomial f(w)
m = size(fw,1) - 1;

% Get degree of polynomial g(w)
n = size(gw,1) - 1;
 

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