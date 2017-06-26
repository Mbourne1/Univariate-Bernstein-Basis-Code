function [ux, vx] = GetQuotients_2Polys(fx, gx, k)
% Given polynomials f(x) and g(x), get the quotient polynomials u(x) and
% v(x) such that f(x)*v(x) = g(x)*u(x).
%
% % Inputs
%
% [fx, gx] : Coefficients of the polynomial f(x) and g(x) in the Bernstein 
% basis
% 
% k : Degree of common divisor
%
% % Outputs
%
% [ux, vx] : Coefficients of the polynomial u(x) and v(x) in the Bernstein
% basis.

global SETTINGS

% Get degree of input polynomial g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the t^th subresultant
St = BuildSubresultant_2Polys(fx, gx, k);

% Get the optimal column for removal
[~,idx_col] = GetMinDistance(St);

% Remove optimal column
At = St;
At(:,idx_col) = [];

% Get the optimal column c_{t} removed from S_{k}
ct = St(:,idx_col);

% Obtain the solution vector x = [-v;u]
x_ls = SolveAx_b(At,ct);

% Insert a zero into the position corresponding to the index of the optimal
% column so that S(f,g)*vec_x = 0.
vec_x =[
    x_ls(1:idx_col-1);
    -1;
    x_ls(idx_col:end);
    ]  ;

% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.

nCoefficients_vx = n - k + 1;

switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        vx_bi = vec_x(1:nCoefficients_vx);
        ux_bi = -vec_x(nCoefficients_vx + 1:end);
        
        ux = GetWithoutBinomials(ux_bi);
        vx = GetWithoutBinomials(vx_bi);
        
    case 'DT'
        
        vx_bi = vec_x(1:nCoefficients_vx);
        ux_bi = -vec_x(nCoefficients_vx + 1:end);
        ux = GetWithoutBinomials(ux_bi);
        vx = GetWithoutBinomials(vx_bi);
        
    case 'DTQ'
        
        vx = vec_x(1:nCoefficients_vx);
        ux = -vec_x(nCoefficients_vx +1:end);
        
    case 'TQ'
        
        
        
        vx = vec_x(1:nCoefficients_vx);
        ux = -vec_x(nCoefficients_vx + 1:end);
        
        
    case 'DTQ Denominator Removed'
        
        vx = vec_x(1:nCoefficients_vx);
        ux = -vec_x(nCoefficients_vx + 1:end);
        
    case 'DTQ Rearranged'
        
        vx = vec_x(1:nCoefficients_vx);
        ux = -vec_x(nCoefficients_vx +1 :end);
        
    otherwise 
        error('err')
end



end