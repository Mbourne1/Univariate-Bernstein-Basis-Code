function [dx] = GetGCD(ux,vx,fx_n,gx_n,t,alpha,theta)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
%                       Inputs
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
m = size(fx_n,1) - 1;

% Get degree of polynomial g(x)
n = size(gx_n,1) - 1;


% Get fw and gw
fw_n = fx_n .* theta.^(0:1:m)';
gw_n = gx_n .* theta.^(0:1:n)';

% Get uw and vw
uw = ux .* theta.^(0:1:m-t)';
vw = vx .* theta.^(0:1:n-t)';

% Build solution vector bk = [f;g]
bk = [fw_n ; alpha .* gw_n];

% Build the coefficient vector HCG
HCG = BuildHCG(uw,vw,m,n,t);

% Get the vector d(w), which is the solution of a problem of the form Ax=b
dw = SolveAx_b(HCG,bk);

dx = dw./(theta.^(0:1:t)');


end




