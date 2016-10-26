function [DTQ] = BuildDTQ(fx,gx,k)
% BuildDTQ(fx,gx,t)
%
% Build the matrix DTQ = D^{-1}T(f,g)*Q.
%
% Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% k : index of subresultant
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
        
        % Build matrix D^{-1}
        D = BuildD(m,n-k);
        
        % Build matrix T(f,g) = T_{n-k}(f) T_{m-k}*(g)
        T = BuildT(fx,gx,k);
        
        % Build matrix Q = [Q_{n-k} Q_{m-k}]
        Q = BuildQ(m,n,k);
        
        % Get D^{-1} * T(f,g) * Q
        DTQ = D*T*Q;
        
    case 'Rearranged'
        
        DTQ = BuildDTQ_ElementWise(fx,gx,k);
        
    otherwise
        error('Error : SETTINGS.SYLVESTER_BUILD_METHOD is either (Standard) or (Rearranged)')
end

end