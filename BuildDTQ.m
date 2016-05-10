function [DTQ] = BuildDTQ(fx,gx,t)
% BuildDTQ(fw,gw,alpha,t)
%
% Build the matrix DTQ = D^{-1}T(f,g)*Q.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% t : Degree of GCD d(x).
%
% Note: If you wish to build the Sylvester matrix for preprocessed and
% scaled polynomials f(\omega)g(\omega), the preprocessed forms must be the
% inputs to this function.

% Global Variables.
global SETTINGS

% Get degree of polynomial f(w)
m = GetDegree(fx);

% Get degree of polynomial g(w)
n = GetDegree(gx);
 

switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'Standard'
        % Build matrices D,T and Q.
        D = BuildD(m,n-t);
        T = BuildT(fx,gx,t);
        Q = BuildQ(m,n,t);
        DTQ = D*T*Q;
        
    case 'Rearranged'
        
        DTQ = BuildDTQ_ElementWise(fx,gx,t);
    otherwise
        error('Error : SETTINGS.SYLVESTER_BUILD_METHOD is either (Standard) or (Rearranged)')
end

end