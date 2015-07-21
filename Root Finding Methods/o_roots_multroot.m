function [roots_Bb] = o_roots_multroot(example_number,e_low,e_up,seed)

addpath 'Examples'
addpath 'multroot/multroot'


% Get the exact roots
[f_roots_exact,~] = Root_Examples(example_number,seed);

% Get polynomial coefficients in scaled bernstein basis
f_exact_bi = B_poly(f_roots_exact);

% Get degree of polynomial f
m = length(f_exact_bi) - 1;

% Get binomial coefficients corresponding to f
Bi_m = zeros(m+1,1);
for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

% Get polynomial coefficients in standard bernstein basis.
f_exact = f_exact_bi./Bi_m ;

% Add variable noise to the coefficients
fx = VariableNoise(f_exact,e_low,e_up,seed);

% Get roots wrt to y
roots_calc = multroot(fx.*Bi_m);

% convert roots wrt t, Bernstein basis
roots_Bb = [1- roots_calc(:,1)./(1+roots_calc(:,1)) roots_calc(:,2)];

fprintf('\nROOTS CALCULATED BY MULTROOTS FUNCTION \n');
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f + %22.15f i  \t \t %3g \n',[real(roots_Bb(:,1)),imag(roots_Bb(:,1)),...
    roots_Bb(:,2)]');



end

