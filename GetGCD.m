function [dx] = GetGCD(ux,vx,fx_n,gx_n,t,alpha,theta)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
% Inputs
%
%
% uw    : Quotient of f where uw is in the form u_{i}\theta^{i}.
%
% vw    : Quotient of g where vw is in the form v_{i}\theta^{i}.
%
% fw_n  : Coefficients of polynomial f in modified bernstein basis.
%
% gw_n  : Coefficients of polynomial g in modified bernstein basis.
%
% t     : Degree of AGCD.
%

% Get degree of polynomials f(x)
m = GetDegree(fx_n);

% Get degree of polynomial g(x)
n = GetDegree(gx_n);

% Get fw and gw
fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);

% Get uw and vw
uw = GetWithThetas(ux,theta);
vw = GetWithThetas(vx,theta);



%
global SETTINGS

switch SETTINGS.GCD_COEFFICIENT_METHOD
    case 'ux and vx'
        % Build solution vector bk = [f;g]
        bk = [fw ; alpha .* gw];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG(uw,vw,m,n,t);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dw = SolveAx_b(HCG,bk);
        
        dx = GetWithoutThetas(dw,theta);
    case 'ux'
        bk = fw;
        
        H1C1G = BuildH1C1G(uw,t);
        
        dw = SolveAx_b(H1C1G,bk);
        
        dx = GetWithoutThetas(dw,theta);
    otherwise
        error('GCD_COEFFICIENT_METHOD is either ux or ux and vx')
end




