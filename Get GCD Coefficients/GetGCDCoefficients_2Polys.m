function [dx] = GetGCDCoefficients_2Polys(ux, vx, fx, gx, t, alpha, theta)
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

% Get f(w) and g(w)
fw = GetWithThetas(fx, theta);
gw = GetWithThetas(gx, theta);

% Get u(w) and v(w)
uw = GetWithThetas(ux, theta);
vw = GetWithThetas(vx, theta);


switch SETTINGS.GCD_COEFFICIENT_METHOD
    case 'ux and vx'
        % Build solution vector bk = [f;g]
        bk = [fw ; alpha .* gw];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG_2Polys(uw, vw, t);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dw = SolveAx_b(HCG,bk);
        
        % Get d(x) without thetas
        dx = GetWithoutThetas(dw,theta);
        
    case 'ux'
        bk = fw;
        
        H1C1G = BuildH1C1G(uw,t);
        
        dw = SolveAx_b(H1C1G,bk);
        
        dx = GetWithoutThetas(dw,theta);
    otherwise
        error('GCD_COEFFICIENT_METHOD is either (ux) or (ux and vx)')
end




