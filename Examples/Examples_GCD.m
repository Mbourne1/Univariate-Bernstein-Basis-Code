function [f_exact,g_exact,d_exact]=Examples_GCD(ex_num)
% Inputs.
%
% n - Index of example to be used
%
% Outputs.
%
% a - Roots and multiplicities of polynomial f.
%
% b - Roots and multiplicities of polynomial g.
%
% c - Roots and multiplicities of polynomial d, the GCD of f and g.
%
% u - Roots and multiplicities of quotient polynomial f/u = d.
%
% v - Roots and multiplicities of quotient polynomial g/v = d.

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Coefficients'
        
        [f_exact,g_exact,d_exact] = Examples_GCD_FromCoefficients(ex_num);
        
        
    case 'From Roots'
        [f_root_mult_arr,g_root_mult_arr,d_root_mult_arr] = Examples_GCD_FromRoots(ex_num);
        
        f_exact_bi = BuildPolyFromRoots(f_root_mult_arr);
        g_exact_bi = BuildPolyFromRoots(g_root_mult_arr);
        d_exact_bi = BuildPolyFromRoots(d_root_mult_arr);

        % Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
        % polynomials f, g, u, v and d in standard bernstein form.
        f_exact = GetWithoutBinomials(f_exact_bi);
        g_exact = GetWithoutBinomials(g_exact_bi);
        d_exact = GetWithoutBinomials(d_exact_bi);

    otherwise
        error('err')
        
end

end



