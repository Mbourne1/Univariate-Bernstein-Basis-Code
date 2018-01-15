function [] = Experiment3VariableNoise_2Polys(ex_num, bool_preproc)
%
% % Experiment : Shows how the degree of the GCD can be computed from the
% singular values, but as noise increases, ability to compute degree of GCD
% decreases.
%


% Variable - Noise level 1e-12 -> 1e-4
% Variable - Include/Exclude Preprocessing

close all;
clc;


rank_revealing_metric = 'Minimum Singular Values';
sylvester_build_method = 'DTQ';
low_rank_approx_method = 'None';
APF_method = 'None';


% -------------------------------------------------------------------------

% Without Preprocessing
switch bool_preproc
    case true
        bool_alpha_theta = true;
        mean_method = 'Geometric Mean Matlab Method';
    case false
        bool_alpha_theta = false;
        mean_method = 'None';
end


el_arr = {1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4};
%eu_arr = {1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2};
for i = 1 : 1 : length(el_arr)
    
    el = el_arr{i};
    eu = el_arr{i};
    
    o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, APF_method, sylvester_build_method, ...
        rank_revealing_metric)
    
end



end