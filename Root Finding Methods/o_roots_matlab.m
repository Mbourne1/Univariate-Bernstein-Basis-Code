function [roots_Bb] = o_roots_matlab(ex,emin,emax,seed)
%% Get roots of Bernstein Basis Polynomial by matlab function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex - (Int) Example Number
% emin - Noise/Signal maximum threshold (minimum)
% emax - Noise/Signal maximum threshold (maximum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'Examples'


% Get the polynomial roots from the example number provided.
[f_roots_exact,~] = Root_Examples(ex,seed);

% Get coefficients of f_exact_bi (Polynomial in scaled Bernstein Basis)
f_exact_bi = B_poly(f_roots_exact);

% Get degree of the polynomial f.
m = length(f_exact_bi) - 1;

% Get the Binomial coefficients corresponding to the coefficients of
% polynomial f.
Bi_m = zeros(m+1,1);
for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

% Divide f_exact_bi (Polynomial in Scaled Bernstein Basis) by binomial
% coefficients to obtain f_exact (Polynomial in Bernstein Basis).
f_exact = f_exact_bi./Bi_m;

% Add variable noise to the coefficients of f_exact and obtain fx in
% the bernstein basis.
fx = VariableNoise(f_exact,emin,emax,seed);

% Get roots wrt to y
roots_calc = roots(fx.*Bi_m);

% Convert the roots to the bernstein basis.
% r_{Bb} :- Roots Bernstein Basis.
% r_{calc} :- Roots Calculated by Matlab Method.
% r_{Bb} = 1-r_{calc} / (1+r_{calc})

fprintf('ROOTS CALCULATED BY MATLAB ROOTS FUNCTION\n');
roots_Bb = [1- roots_calc./(1+roots_calc) ones(length(roots_calc),1)];
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f + %22.15f i  \t \t %3g \n',[real(roots_Bb(:,1)),imag(roots_Bb(:,1)),...
    roots_Bb(:,2)]');
fprintf('\n');



end