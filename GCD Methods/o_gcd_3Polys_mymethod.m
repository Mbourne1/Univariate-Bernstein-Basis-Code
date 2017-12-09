function [fx_o, gx_o, hx_o, dx_o, ux_o, vx_o, wx_o,  lambda_o, mu_o, rho_o, t ] = ...
    o_gcd_3Polys_mymethod(fx, gx, hx, limits_t, rank_range)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
% Inputs:
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% limits_t : [Int Int] Upper and lower limits for GCD degree may be defined here
% otherwise set to [0,min(m,n)]
%
%
% % Outputs:
%
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% dx : (Vector) The GCD of f(x) + \delta f(x) and g(x) + \delta g(x)
%
% ux : (Vector) Coefficients of polynomial u(x) = f(x)/d(x)
%
% vx : (Vector) Coefficients of polynomial v(x) = g(x)/d(x)
%
% wx : (Vector) Coefficients of polynomial w(x) = h(x)/d(x)
%
% lambda : (Float) Optimal \lambda
%
% mu : (Float) Optimal \mu
%
% rho : (Float) Optimal \rho
%
% theta : (Float) Optimal \theta






% % Get the degree of the GCD of f(x) g(x) and h(x)
[t, lambda, mu, rho, theta, gm_fx, gm_gx, gm_hx] = ...
    Get_GCD_Degree_3Polys(fx, gx, hx, limits_t, rank_range);

LineBreakLarge();


if t == 0 % If degree of GCD is 0, polynomials are coprime
    
    fprintf([mfilename ' : ' sprintf('f(x) and g(x) are coprime \n')])
    
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    
    lambda_o = 1;
    mu_o = 1; 
    rho_o = 1;
    
    return
    
end

% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end

% Normalise f(x) and g(x) by Geometric mean to obtain fx_n and gx_n.
% Normalise by geometric mean obtained by entries of f(x) and g(x) in the
% subresultant S_{t}
fx_n = fx ./ gm_fx;
gx_n = gx ./ gm_gx;
hx_n = hx ./ gm_hx;

% % Get the optimal column of the sylvester matrix to be removed. Where
% % removal of the optimal column gives the minmal residual in (Ak x = ck)

% Get f(\omega) and \alpha.*g(\omega)
lambda_fw = lambda .* GetWithThetas(fx_n, theta);
mu_gw  = mu .* GetWithThetas(gx_n, theta);
rho_hw = rho .* GetWithThetas(hx_n, theta);


% Build S_{t}(f,g,h)
St_preproc = BuildSubresultant_3Polys(lambda_fw, mu_gw, rho_hw, t);

% Get index of optimal column for removal
[~, idx_col] = GetMinDistance(St_preproc);

% % Get Low rank approximation of the Sylvester matrix S_{t}
[fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, lambda_lr, mu_lr, ...
    rho_lr, theta_lr] = LowRankApproximation_3Polys(fx_n, ...
    gx_n, hx_n, lambda, mu, rho, theta, t, idx_col);

% Get the coefficients of the GCD by APF or other method.
[fx_alr, gx_alr, hx_alr, dx_alr, ux_alr, vx_alr, wx_alr, alpha_alr, beta_alr, theta_alr] = ...
    APF_3Polys(fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, lambda_lr, mu_lr, rho_lr, theta_lr, t);

% Get outputs
fx_o = fx_alr;
gx_o = gx_alr;
hx_o = hx_alr;

dx_o = dx_alr;

ux_o = ux_alr;
vx_o = vx_alr;
wx_o = wx_alr;

lambda_o = alpha_alr;
mu_o = beta_alr;
rho_o = theta_alr;




end








