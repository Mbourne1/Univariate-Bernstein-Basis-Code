function [dx] = GetGCDCoefficients_3Polys(ux, vx, wx, fx, gx, hx, k, alpha, beta, gamma, theta)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
% % Inputs
%
% ux : (Vector) Coefficients of the polynomial u(x)
%
% vx : (Vector) Coefficients of the polynomial v(x)
%
% wx : (Vector) Coefficients of the polynomial w(x)
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% k : (Int) Degree of common divisor.
%
% alpha : (Float) Optimal value of \alpha
%
% beta : (Float) Optimal value of \beta
%
% gamma : (Float) Optimal value of \gamma
%
% theta : (Float) Optimal value of \theta
%
% % Outputs
%
% dx : (Vector) Coefficients of the polynomial d(x)

if (nargin ~= 11)
   error('Not enough input arguments'); 
end

% Global variables
global SETTINGS

% Get f(\omega) and g(\omega) and h(\omega)
fw = GetWithThetas(fx, theta);
gw = GetWithThetas(gx, theta);
hw = GetWithThetas(hx, theta); 

% Get u(\omega), v(\omega) and w(\omega)
uw = GetWithThetas(ux, theta);
vw = GetWithThetas(vx, theta);
ww = GetWithThetas(wx, theta);

switch SETTINGS.GCD_COEFFICIENT_METHOD
    
    case 'ux and vx'
        % Build solution vector bk = [f;g]
        bk = [alpha.*fw ; beta.*gw; gamma.*hw];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG_3Polys(uw, vw, ww, k);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dw = SolveAx_b(HCG, bk);
        
        % Get d(x) without thetas
        dx = GetWithoutThetas(dw, theta);
        
    case 'ux'
        bk = fw;
        
        H1C1G = BuildH1C1G(uw, k);
        
        dw = SolveAx_b(H1C1G, bk);
        
        dx = GetWithoutThetas(dw, theta);
    otherwise
        error('GCD_COEFFICIENT_METHOD is either (ux) or (ux and vx)')
end




