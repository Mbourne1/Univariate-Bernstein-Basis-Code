function [fx,gx,dx] = Examples_GCD_FromCoefficients(ex_num)

x = sym('x');

addpath('../Examples')
[f,g,d,u,v] = Univariate_GCD_Examples(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
fx = GetCoefficients(f);
gx = GetCoefficients(g);
dx = GetCoefficients(d);



end


function fx = GetCoefficients(f_factor_mult_array)

addpath 'Basis Conversion'
addpath 'Bernstein Methods'


% This function takes the factors of the polynomial f(x) in a symbolic
% form, gets the coefficients of the factors in Bernstein form, then
% multiplies these to form polynomial f(x) in bernstein form.

syms x

nFactors = size(f_factor_mult_array,1);

for i = 1 : 1 : nFactors
   
   % Get the factor
   sym_factor = f_factor_mult_array(i,1);
   
   % Get the multiplicity
   mult = f_factor_mult_array(i,2);
   
   % Get the expression (factor).^{mult}
   sym_factor = sym_factor^mult;
   
   % Get coefficients in power basis
   pwr_poly = double(fliplr(coeffs(sym_factor,x,'All')))';
   
   
   % Convert the factor to a polynomial in Bernstein form
   arr_factors_f{i} = PowerToBernstein(pwr_poly); 
   
end

% Multiply all matrices of coefficients of factors to form coefficients of
% polynomial f(x)

fx = arr_factors_f{1};

for i = 2:1:nFactors
   fx = Bernstein_Multiply(fx,arr_factors_f{i});
   
end

end
