function [] = Experiment2Preprocessing_3Polys(ex_num)

close all; clc;


% Good Examples

% 1
% 2
% 3
% 4
% 5
% 6


% 8a -
% 12 - too big
% 13 - too big
% 14 - too big
% 15 - too big



% Great Example
% Example 8
% Note - Coefficients of g(x) span many more orders of magnitude than f or
% h, and g(x) has much higher degree.
% 8a 1e-5, 1e-7 - Preprocessing recovers gcd degree,
% 8b

emin = 1e-7;
emax = 1e-4;

% Constants
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Normalised Minimum Singular Values';



global SETTINGS

SETTINGS.SCALING_METHOD = 'NONE';

nEquations = '2';
mean_method = 'None';
bool_alpha_theta = false;
o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, nEquations, rank_revealing_metric)


bool_alpha_theta = true;
mean_method = 'Geometric Mean Matlab Method';
%arrScaling_method = {'lambda_mu_rho', 'mu_rho', 'lambda_mu','lambda_rho'};
arrScaling_method = {'lambda_rho'};

for i = 1 : 1 : length(arrScaling_method)
    
    SETTINGS.SCALING_METHOD = arrScaling_method{i};
    o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, nEquations, rank_revealing_metric)
    
end

end