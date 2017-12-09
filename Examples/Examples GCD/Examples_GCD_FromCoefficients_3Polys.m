function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_FromCoefficients_3Polys(ex_num_var)
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

addpath(genpath('../Examples'))

ex_num = ex_num_var(1 : end - 1);
ex_num_variant = ex_num_var(end);

% Get roots and multiplicities from example file
[PolyA_root_mult_arr, PolyB_root_mult_arr, PolyC_root_mult_arr, ...
    d_root_mult_arr,...
    Poly1_root_mult_arr, Poly2_root_mult_arr, Poly3_root_mult_arr] = ...
    GCD_Examples_Univariate_3Polys(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and h(x).
PolyA_sym = GetSymbolicPolyFromSymbolicRoots(PolyA_root_mult_arr);
PolyB_sym = GetSymbolicPolyFromSymbolicRoots(PolyB_root_mult_arr);
PolyC_sym = GetSymbolicPolyFromSymbolicRoots(PolyC_root_mult_arr);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
PolyA = BuildPolyFromRootsSymbolic(PolyA_root_mult_arr);
PolyB = BuildPolyFromRootsSymbolic(PolyB_root_mult_arr);
PolyC = BuildPolyFromRootsSymbolic(PolyC_root_mult_arr);

Poly1 = BuildPolyFromRootsSymbolic(Poly1_root_mult_arr);
Poly2 = BuildPolyFromRootsSymbolic(Poly2_root_mult_arr);
Poly3 = BuildPolyFromRootsSymbolic(Poly3_root_mult_arr);

% Get the coefficients of the polynomial d(x)
dx = BuildPolyFromRootsSymbolic(d_root_mult_arr);
dx_sym = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);

bad_scaling = false;

if bad_scaling == true
   
    PolyA = PolyA * 10^(5);
    PolyB = PolyB * 10^(5);
    PolyC = PolyC * 10^(-5);
end

switch ex_num_variant
    
    case 'a'
        
        fx = PolyA;
        fx_sym = PolyA_sym;
        
        gx = PolyB;
        gx_sym = PolyB_sym;
        
        hx = PolyC;
        hx_sym = PolyC_sym;
        
        ux = Poly1;
        vx = Poly2;
        wx = Poly3;
        
        
    case 'b'
        
        gx = PolyA;
        gx_sym = PolyA_sym;
        
        fx = PolyB;
        fx_sym = PolyB_sym;
        hx = PolyC;
        hx_sym = PolyC_sym;
        vx = Poly1;
        ux = Poly2;
        wx = Poly3;
        
    case 'c'
        
        hx = PolyA;
        hx_sym = PolyA_sym;
        
        gx = PolyB;
        gx_sym = PolyB_sym;
        
        fx = PolyC;
        fx_sym = PolyC_sym;
        
        wx = Poly1;
        vx = Poly2;
        ux = Poly3;
        
    otherwise
        error('Not valid ordering')
end


% Display polynomials f(x) g(x) and h(x)
display(fx_sym)
display(gx_sym)
display(hx_sym)

% Display common divisor d(x)
display(dx_sym)


end


