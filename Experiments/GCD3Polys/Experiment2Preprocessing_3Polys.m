function [] = Experiment2Preprocessing_3Polys(ex_num)

close all; clc;


% Constants for this experiment -------------------------------------------

emin = 1e-7;
emax = 1e-4;

% Constants
% low_rank_approx_method 
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';


apf_method = 'None';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

% rank_revealing_metric
%   'Normalised Minimum Singular Values'
%   'Minimum Singular Values'
rank_revealing_metric = 'Normalised Minimum Singular Values';



global SETTINGS

SETTINGS.SCALING_METHOD = 'NONE';

nEquations = '2';
mean_method = 'None';
bool_alpha_theta = false;
o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, nEquations, rank_revealing_metric)


bool_alpha_theta = true;
mean_method = 'Geometric Mean Matlab Method';
%arrScaling_method = {'lambda_mu_rho', 'mu_rho', 'lambda_mu','lambda_rho'};
arrScaling_method = {'lambda_rho'};

for i = 1 : 1 : length(arrScaling_method)
    
    SETTINGS.SCALING_METHOD = arrScaling_method{i};
    o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, nEquations, rank_revealing_metric)
    
end

end