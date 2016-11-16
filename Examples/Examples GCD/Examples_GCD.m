function [fx,gx,dx,ux,vx] = Examples_GCD(ex_num)
% Inputs. 
%
% ex_num : Example number
%
% Outputs.
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% dx : Coefficients of the polynomial d(x), the GCD of f(x) and g(x)
%
% ux : Coefficients of the polynomial u(x), given by f(x)/d(x)
%
% vx : Coefficients of the polynomial v(x), given by g(x)/d(x)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Coefficients'
        
        [fx,gx,dx,ux,vx] = Examples_GCD_FromCoefficients(ex_num);
        
        
    case 'From Roots'
        [f_root_mult_arr,g_root_mult_arr,d_root_mult_arr] = Examples_GCD_FromRoots(ex_num);
        
        f_exact_bi = BuildPolyFromRoots(f_root_mult_arr);
        g_exact_bi = BuildPolyFromRoots(g_root_mult_arr);
        d_exact_bi = BuildPolyFromRoots(d_root_mult_arr);

        % Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
        % polynomials f, g, u, v and d in standard bernstein form.
        fx = GetWithoutBinomials(f_exact_bi);
        gx = GetWithoutBinomials(g_exact_bi);
        dx = GetWithoutBinomials(d_exact_bi);

    otherwise
        error('err')
        
end

end



