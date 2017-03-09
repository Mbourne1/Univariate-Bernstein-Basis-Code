function [fx_o, gx_o, dx_o, ux_o, vx_o, alpha_o, theta_o, t ] = ...
    o_gcd_2Polys_mymethod(fx, gx, deg_limits)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
% Inputs:
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% deg_limits : [(Int) (Int)] Upper and lower limits for GCD degree may be defined here
% otherwise set to [0,min(m,n)]
%
%
%
% Outputs:
%
%
% fx : (Vector) f(x) + \delta f(x)
%
% gx : (Vector) g(x) + \delta g(x
%
% dx : (Vector) The GCD of f(x) + \delta f(x) and g(x) + \delta g(x)
%
% ux : (Vector) Coefficients of polynomial u(x) where u(x) = f(x)/d(x)
%
% vx : (Vector) Coefficeints of polynomial v(x) where v(x) = g(x)/d(x)
%
% alpha : (Float) Optimal \alpha
%
% theta : (Float) Optimal \theta




% Use this to compute the degree of the GCD only using one Sylvester matrix
% onesub = 'y'. This code is included as a test to show that computing GCD degree by
% all subresultant matrices is better than just the first.
onesub = 'n';
if onesub == 'y'
    Get_GCD_Degree_One_Subresultant(fx,gx);
    fx_o = [];
    gx_o = [];
    dx_o = [];
    ux_o = [];
    vx_o = [];
    alpha_o = [];
    theta_o = [];
    t = [];
    return;
end

% % Get the degree of the GCD

[t, alpha, theta, GM_fx, GM_gx] = Get_GCD_Degree_2Polys(fx, gx, deg_limits);
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

% Normalise f(x) and g(x) by Geometric mean to obtain fx_n and gx_n.
% Normalise by geometric mean obtained by entries of f(x) and g(x) in the
% subresultant S_{t}
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;

% % Get the optimal column of the sylvester matrix to be removed. Where
% % removal of the optimal column gives the minmal residual in (Ak x = ck)

% Get f(\omega) and \alpha.*g(\omega)
fw = GetWithThetas(fx_n, theta);
a_gw = alpha.* GetWithThetas(gx_n, theta);

% Build S_{t}(f,g)
St_preproc = BuildSubresultant_2Polys(fw, a_gw, t);

% Get index of optimal column for removal
[~,idx_col] = GetMinDistance(St_preproc);

% % Get Low rank approximation of the Sylvester matrix S_{t}(f,g)
% suffix 'lr' stands for 'low rank'
[fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    LowRankApproximation_2Polys(fx_n, gx_n, alpha, theta, t, idx_col);

% Get the coefficients of the GCD by APF or other method.
% suffix alr stands for 'approx low rank' 
[fx_alr, gx_alr, dx_alr, ux_alr, vx_alr, alpha_alr, theta_alr] = ...
    APF_2Polys(fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr, t);

% Get outputs
fx_o = fx_alr;
gx_o = gx_alr;
dx_o = dx_alr;
ux_o = ux_alr;
vx_o = vx_alr;
alpha_o = alpha_alr;
theta_o = theta_alr;




end








