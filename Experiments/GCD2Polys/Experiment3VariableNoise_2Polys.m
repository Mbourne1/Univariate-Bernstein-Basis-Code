function [] = Experiment3VariableNoise_2Polys(ex_num)
%
% % Experiment : Shows how the degree of the GCD can be computed from the
% singular values, but as noise increases, ability to compute degree of GCD
% decreases.
%


% Variable - Noise level 1e-12 -> 1e-4
% Variable - Include/Exclude Preprocessing

close all;
clc;

% -------------------------------------------------------------------------

% Without Preprocessing
bool_alpha_theta = false;
mean_method = 'None';

el_arr = {1e-12, 1e-10, 1e-8, 1e-6, 1e-4};

for i = 1 : 1 : length(el_arr)
    
    el = el_arr{i};
    o_gcd_Univariate_2Polys(ex_num, el, el, mean_method, bool_alpha_theta, 'None', 'None', 'DTQ')
    
end


% -------------------------------------------------------------------------

% With Preprocessing
bool_alpha_theta = true;
mean_method = 'Geometric Mean Matlab Method';

for i = 1 : 1 : length(el_arr)
    
    el = el_arr{i};
    o_gcd_Univariate_2Polys(ex_num, el, el, mean_method, bool_alpha_theta, 'None', 'None', 'DTQ')
    
end

end