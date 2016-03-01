function [roots_Bb] = o_roots_matlab(fx)
% Get roots of a polynomial whose coefficients are given in Bernstein Basis
% Polynomial by an adaptation of the MATLAB Roots() function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs.
%
% fx : Coefficients of polynomial f(x) whose roots are to be computd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the degree of polynomial f
[r,~] = size(fx);
m = r - 1;

% Get the Binomial coefficients corresponding to the coefficients of
% polynomial f.
Bi_m = GetBinomials(m);

% Get roots wrt to y
roots_calc = roots(fx.*Bi_m);

% Convert the roots to the bernstein basis.
% r_{Bb} :- Roots Bernstein Basis.
% r_{calc} :- Roots Calculated by Matlab Method.
% r_{Bb} = 1-r_{calc} / (1+r_{calc})

roots_Bb = [1- roots_calc./(1+roots_calc) ones(length(roots_calc),1)];

% Printout calculated roots by matlab method
PrintoutRoots('MATLAB', roots_Bb);



end