function [fx,gx,dx] = Examples_GCD_FromCoefficients(ex_num)


addpath('../Examples')

% Get roots and multiplicities from example file
[f_root_mult_arr,g_root_mult_arr,d_root_mult_arr,u_root_mult_arr,v_root_mult_arr] = Univariate_GCD_Examples(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);

gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);

dx = BuildPolyFromRootsSymbolic(d_root_mult_arr);



end


