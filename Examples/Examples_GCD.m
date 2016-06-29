
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

EXAMPLE_TYPE = 'From Roots';

switch EXAMPLE_TYPE
    case 'From Coefficients'
        [f_exact,g_exact,d_exact] = FromCoefficients(ex_num);
        
        
    case 'From Roots'
        [f_roots,g_roots,d_roots] = FromRoots(ex_num);
        
        f_exact_bi = B_poly(f_roots);
        g_exact_bi = B_poly(g_roots);
        d_exact_bi = B_poly(d_roots);

        % Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
        % polynomials f, g, u, v and d in standard bernstein form.
        f_exact = GetWithoutBinomials(f_exact_bi);
        g_exact = GetWithoutBinomials(g_exact_bi);
        d_exact = GetWithoutBinomials(d_exact_bi);

    otherwise
        error('err')
        
end

end

function [f,g,d] = FromCoefficients(ex_num)

    switch ex_num
        case '1'
            x = sym('x');
            
            f = (x+1)^2 * (x+2)^2;
            g = (x+1)^2 * (x+3)^2;
            d = (x+1)^2;
            
        f = double(fliplr(coeffs(f,x,'All')))';
        g = double(fliplr(coeffs(g,x,'All')))';
        d = double(fliplr(coeffs(d,x,'All')))';
            
        f = PowerToBernstein(f);
        g = PowerToBernstein(g);
        d = PowerToBernstein(d);
    end

end

function [root_mult_array_fx, root_mult_array_gx , root_mult_array_dx] = FromRoots(ex_num)


pattern = 'Custom:m=(\d+).n=(\d+).t=(\d+).low=(-?\d+).high=(-?\d+)';

if ~isempty(regexp(ex_num,pattern,'start'))
    
    str = ex_num;
    
    expression_m = regexp(str,'m=(\d+)','tokens');
    m_str = expression_m{1};
    
    expression_n = regexp(str,'n=(\d+)','tokens');
    n_str = expression_n{1};
    
    expression_t = regexp(str,'t=(\d+)','tokens');
    t_str = expression_t{1};
    
    expression_low = regexp(str, 'low=(-?\d+)', 'tokens');
    low_str = expression_low{1};
    
    expression_high = regexp(str, 'high=(-?\d+)','tokens');
    high_str = expression_high{1};
    


    m = str2double(m_str);
    n = str2double(n_str);
    t = str2double(t_str);
    intvl_low = str2double(low_str);
    intvl_high = str2double(high_str);
    
    [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
    
    
    
    
else
    switch ex_num
        
        
        case '1'
            root_mult_array_fx = [
                0.2 2
                0.5 1
                0.7 1
                0.9 1
                1.1 1
                ];
            root_mult_array_gx = [
                0.5 1
                0.1 1
                0.9 1
                0.3 1
                ];
            
            
        case '2'
            root_mult_array_fx = [
                0.2 1
                0.4 1
                0.6 1
                0.8 1
                ];
            
            root_mult_array_gx = [0.2 1
                0.3 1];
            
        case '3'
            root_mult_array_fx = [
                0.56 20
                0.75 3
                0.82 3
                0.37 3];
            
            root_mult_array_gx = [
                0.56    20
                0.75    3
                0.99    4
                0.37    3
                0.12    3
                0.20    3
                ];
            
        case '4'
            
            root_mult_array_fx = [0.1    20
                0.5    2];
            
            root_mult_array_gx = [0.1    20
                0.9    1];
            
        case '5'
            
            root_mult_array_fx = [0.1    2
                0.3     2
                0.5    2];
            
            root_mult_array_gx = [0.1    2];
            
            % From Winkler Paper - Methods for the computation of the degree of an
            % approximate greatest common divisor of two inexact bernstein basis
            % polynomials.
        case '6'
            root_mult_array_fx = [
                0.10    3
                0.56    4
                0.75    3
                0.82    3
                1.37    3
                -0.27   3
                1.46    2
                ];
            
            root_mult_array_gx = [
                0.10    2
                0.56    4
                0.75    3
                0.99    4
                1.37    3
                2.12    3
                1.20    3
                ];
            
        case '7'
            root_mult_array_fx = [
                0.23   4
                0.43   3
                0.57   3
                0.92   3
                1.70   3
                ];
            root_mult_array_gx = [
                0.23   4
                0.30   2
                0.77   5
                0.92   2
                1.20   5
                ];
            
        case '8'
            root_mult_array_fx = [0.1  5
                0.56 4
                0.75 3
                0.82 3
                1.37 3];
            root_mult_array_gx = [0.1  5
                0.56 4
                0.75 3
                0.99 4
                1.37 3
                2.12 3
                1.20 3];
            
            % Example 6.2
        case '9'
            root_mult_array_fx = [
                0.14    3
                0.56    3
                0.89    4
                1.45    4
                2.37    3
                -3.61   5
                ];
            root_mult_array_gx = [
                0.14    4
                0.99    1
                2.37    3
                -0.76   2
                -1.24   2
                -3.61   7
                ];
            
        case '10'
            root_mult_array_fx = [
                0.14    3
                0.56    3
                0.89    4
                0.37    3
                ];
            
            root_mult_array_gx = [
                0.14    3
                0.37    3
                0.76   2
                0.24   2
                ];
            
            
            % Example 6.2
        case '11'
            root_mult_array_fx = [
                0.10    3
                0.56    5
                1.40    4
                1.79    3
                2.69    2
                -2.68   3
                ];
            root_mult_array_gx = [
                0.10    4
                0.56    4
                1.79    3
                2.68    2
                2.69    3
                -1.40   2
                ];
        case '12'
            root_mult_array_fx = [
                0.10    7;
                0.50    12;
                ];
            
            root_mult_array_gx = [
                0.10    6;
                0.50    11;
                0.99    1;
                ];
        case '13'
            root_mult_array_fx = ...
                [
                0.5     4;
                ];
            root_mult_array_gx = ...
                [
                0.5     3;
                ]
            
        case '13'
            
            intvl_low = -1;
            intvl_high = 1;
            
            t = 5;
            m = 10;
            ex_num = 7;
            [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,ex_num,t,intvl_low, intvl_high);
            
        case 'Custom'
            intvl_low = -1;
            intvl_high = 1;
            
            prompt = 'Enter the degree of Polynomial f(x) :';
            m = input(prompt);
            
            prompt = 'Enter the degree of Polynomial g(x) :';
            ex_num = input(prompt);
            
            prompt = 'Enter the degree of Polynomial d(x) :';
            t = input(prompt);
            
            
            %         t = 5;
            %         m = 10;
            %         n = 7;
            [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,ex_num,t,intvl_low, intvl_high)
            
            
        otherwise
            error('error: Example Number not valid')
            
    end
    
end
root_mult_array_dx = GetGCDRoots(root_mult_array_fx,root_mult_array_gx)
root_mult_array_ux = GetQuotientRoots(root_mult_array_fx,root_mult_array_dx);
root_mult_array_vx = GetQuotientRoots(root_mult_array_gx,root_mult_array_dx);

m = sum(root_mult_array_fx(:,2));
ex_num = sum(root_mult_array_gx(:,2));
d = sum(root_mult_array_dx(:,2));
t = d;

end




