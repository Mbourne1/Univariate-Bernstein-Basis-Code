function [fx_o, gx_o, hx_o, dx_o, ux_o, vx_o, wx_o, alpha_o, theta_o, t ] = ...
    o_gcd_3Polys_mymethod(fx, gx, hx, deg_limits)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
% Inputs:
%
%
% [fx, gx, hx] : Coefficients of the polynomial f(x), g(x) and h(x)
%
% deg_limits : Upper and lower limits for GCD degree may be defined here
% otherwise set to [0,min(m,n)]
%
%
%
% Outputs:
%
%
% [fx, gx, hx] : Coefficients of polynomials    
%                   f(x) + \delta f(x)
%                   g(x) + \delta g(x)
%                   h(x) + \delta h(x)
%
% dx : The GCD of f(x) + \delta f(x) and g(x) + \delta g(x)
%
% [ux, vx, wx] : Coefficients of polynomials
%                   u(x) where u(x) = f(x)/d(x)
%                   v(x) where v(x) = g(x)/d(x)
%                   w(x) where w(x) = h(x)/d(x)
%
% alpha : Optimal \alpha
%
% theta : Optimal \theta




% % Get the degree of the GCD of f(x) g(x) and h(x)
[t, alpha, theta, gm_fx, gm_gx, gm_hx] = ...
    Get_GCD_Degree_3Polys(fx, gx, hx, deg_limits);


LineBreakLarge();


if t == 0 % If degree of GCD is 0, polynomials are coprime
    
    fprintf([mfilename ' : ' sprintf('f(x) and g(x) are coprime \n')])
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    alpha_o = 1;
    theta_o = 1;
    return
    
end

% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end

% Normalise f(x) and g(x) by Geometric mean to obtain fx_n and gx_n.
% Normalise by geometric mean obtained by entries of f(x) and g(x) in the
% subresultant S_{t}
fx_n = fx./gm_fx;
gx_n = gx./gm_gx;
hx_n = hx./gm_hx;

% % Get the optimal column of the sylvester matrix to be removed. Where
% % removal of the optimal column gives the minmal residual in (Ak x = ck)

% Get f(\omega) and \alpha.*g(\omega)
fw = GetWithThetas(fx_n,theta);
a_gw = alpha.*GetWithThetas(gx_n,theta);
hw = GetWithThetas(hx_n,theta);


% Build S_{t}(f,g)
St_preproc = BuildSubresultant_3Polys(fw, a_gw, hw, t);

% Get index of optimal column for removal
[~,idx_col] = GetMinDistance(St_preproc);

% % Get Low rank approximation of the Sylvester matrix S_{t}
[fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
    LowRankApproximation_3Polys(fx_n, gx_n, hx_n, alpha, theta, t, idx_col);

% Get the coefficients of the GCD by APF or other method.
[fx_alr, gx_alr, hx_alr, dx_alr, ux_alr, vx_alr, wx_alr, alpha_alr, theta_alr] = ...
    APF_3Polys(fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr, t);

% Get outputs
fx_o = fx_alr;
gx_o = gx_alr;
hx_o = hx_alr;

dx_o = dx_alr;

ux_o = ux_alr;
vx_o = vx_alr;
wx_o = wx_alr;

alpha_o = alpha_alr;
theta_o = theta_alr;




end








