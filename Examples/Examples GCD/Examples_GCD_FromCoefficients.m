function [fx,gx,dx,ux,vx] = Examples_GCD_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
% % Outputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% dx : (Vector) Coefficients of the GCD d(x)
%
% ux : (Vector) Coefficients of the polynomial u(x) given by f(x)/d(x)
%
% vx : (Vector) Coefficients of the polynomial v(x) given by g(x)/d(x)

addpath(genpath('../Examples'));

% Get roots and multiplicities from example file
[f_root_mult_arr, g_root_mult_arr, d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ...
    GCD_Examples_Univariate_2Polys(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);
gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);
dx = BuildPolyFromRootsSymbolic(d_root_mult_arr);
ux = BuildPolyFromRootsSymbolic(u_root_mult_arr);
vx = BuildPolyFromRootsSymbolic(v_root_mult_arr);

% Get the symbolic polynomials
fx_sym = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);
dx_sym = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);
ux_sym = GetSymbolicPolyFromSymbolicRoots(u_root_mult_arr);
vx_sym = GetSymbolicPolyFromSymbolicRoots(v_root_mult_arr);

display(fx_sym)
display(gx_sym)
display(dx_sym)
display(ux_sym)
display(vx_sym)

end


