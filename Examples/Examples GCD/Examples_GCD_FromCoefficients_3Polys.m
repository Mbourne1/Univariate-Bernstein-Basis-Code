function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_FromCoefficients_3Polys(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
% % Outputs
%
% [fx, gx, hx] : Coefficients of the polynomials 
%                   f(x)
%                   g(x)
%                   d(x)
%
% dx : Coefficients of the GCD d(x)
%
% [ux, vx, wx] : Coefficients of the polynomials 
%                   u(x) given by f(x)/d(x)
%                   v(x) given by g(x)/d(x)
%                   w(x) given by h(x)/d(x)

addpath('../Examples')

% Get roots and multiplicities from example file
[f_root_mult_arr, g_root_mult_arr, h_root_mult_arr, ...
    d_root_mult_arr,...
    u_root_mult_arr, v_root_mult_arr, w_root_mult_arr] = Univariate_GCD_Examples_3Polys(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and h(x).
fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);
gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);
hx = BuildPolyFromRootsSymbolic(h_root_mult_arr);

% Get the coefficients of the polynomial d(x)
dx = BuildPolyFromRootsSymbolic(d_root_mult_arr);

% Get the coefficients of the polynomials u(x) v(x) and w(x)
ux = BuildPolyFromRootsSymbolic(u_root_mult_arr);
vx = BuildPolyFromRootsSymbolic(v_root_mult_arr);
wx = BuildPolyFromRootsSymbolic(w_root_mult_arr);

% Get the symbolic polynomials
fx_sym = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);
hx_sym = GetSymbolicPolyFromSymbolicRoots(h_root_mult_arr);

dx_sym = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);

ux_sym = GetSymbolicPolyFromSymbolicRoots(u_root_mult_arr);
vx_sym = GetSymbolicPolyFromSymbolicRoots(v_root_mult_arr);
wx_sym = GetSymbolicPolyFromSymbolicRoots(w_root_mult_arr);

% display polynomials f(x) g(x) and h(x)
display(fx_sym)
display(gx_sym)
display(hx_sym)

% Display common divisor d(x)
display(dx_sym)

% Display cofactor polynomials u(x) v(x) and w(x)
display(ux_sym)
display(vx_sym)
display(wx_sym)
end


