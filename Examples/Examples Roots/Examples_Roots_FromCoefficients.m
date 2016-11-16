
function fx = Examples_Roots_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
%
% % Outputs
%
% fx : Coefficients of the polynomial f(x) in the Bernstein form.

addpath('../Examples');

% Get the factors and corresponding multiplicities of f(x)
f_root_mult_array = Univariate_Roots_Examples(ex_num);


% Get the coefficients of f(x) in Bernstein form
fx = BuildPolyFromRootsSymbolic(f_root_mult_array);
fx_sym = BuildSymbolicPolyFromSymbolicRoots(f_root_mult_array);
disp(fx_sym);




end

