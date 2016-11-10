
function fx = Examples_Roots_FromCoefficients(ex_num)

x = sym('x');

% Get the factors and multiplicities of f(x)
addpath('../Examples');
f_root_mult_array = Univariate_Roots_Examples(ex_num);


% Get the coefficients of f(x) in Bernstein form
fx = GetCoefficients(f_root_mult_array);




end


function [fx] = GetCoefficients(f_roots_mult_array)
% Get the coefficients of the polynomial f(x) in Bernstein form, when given
% the factors and the multiplicity of each factor in f(x).
%
% Inputs.
%
% f_factor_mult_array : A matrix whose rows contain [factor, mult] pairs.
% Each factor is given as a symbolic expression.


addpath 'Basis Conversion'
addpath 'Bernstein Methods'

syms x

% Get the number of factors
nFactors = size(f_roots_mult_array);

% For each factor
for i = 1:1:nFactors
    
    % Get the factor
    sym_factor = f_roots_mult_array(i,1);
    
    % Get the multiplicity
    mult = f_roots_mult_array(i,2);
    
    % Get a temp symbolic polynomial of (x-r_{i})^(m_{i})
    temp_poly = sym_factor.^mult;
    
    % Get coefficients of the symbolic polynomial
    temp_poly_coeffs = flipud(double(coeffs(temp_poly,x,'All'))');
    
    % Store the coefficients of the temp poly in an array
    arr_factors_brn{i} = PowerToBernstein(temp_poly_coeffs);
    
end

% Get the product of all polynomials (x-r_{i})^{m_{i}}
poly = arr_factors_brn{1};


for i = 2:1:nFactors
    
    poly = Bernstein_Multiply(poly,arr_factors_brn{i});
    
end

    fx = poly;

end