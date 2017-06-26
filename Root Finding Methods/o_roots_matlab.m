function [root_mult_mat] = o_roots_Matlab(fx)
% Get roots of a polynomial whose coefficients are given in Bernstein Basis
% Polynomial by an adaptation of the MATLAB Roots() function
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x) whose roots are to be computd
% 
% % Outputs
%
% root_mult_mat : (Matrix) Matrix containing the roots of f(x) and 
% corresponding multiplicities

% Get the Binomial coefficients corresponding to the coefficients of
% polynomial f.
fx_bi = GetWithBinomials(fx);

% Get roots wrt to y
roots_calc = roots(fx_bi);

% Convert the roots to the bernstein basis.
% r_{Bb} :- Roots Bernstein Basis.
% r_{calc} :- Roots Calculated by Matlab Method.
% r_{Bb} = 1-r_{calc} / (1+r_{calc})

root_mult_mat = [1- roots_calc./(1 + roots_calc) ones(length(roots_calc),1)];


end