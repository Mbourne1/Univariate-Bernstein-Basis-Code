function [dx] = GetGCDCoefficients_2Polys(ux, vx, fx, gx, t)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
% % Inputs
%
% [ux, vx] : Coefficients of cofactor polynomials u(x) and v(x) in the
% Bernstein basis.
%
% [fx, gx] : Coefficients of polynomial f(x) and g(x) in the Bernstein basis.
%
% k : Degree of common divisor.
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta


% Global variables
global SETTINGS


switch SETTINGS.GCD_COEFFICIENT_METHOD
    case 'ux and vx'
        % Build solution vector bk = [f;g]
        bk = [fx ; gx];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG_2Polys(ux, vx, t);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dx = SolveAx_b(HCG,bk);
        
    case 'ux'
        bk = fx;
        
        H1C1G = BuildH1C1G(ux,t);
        
        dx = SolveAx_b(H1C1G,bk);
        
        
    otherwise
        error('GCD_COEFFICIENT_METHOD is either (ux) or (ux and vx)')
end




