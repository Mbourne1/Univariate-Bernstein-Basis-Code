function [] = Experiment2Preprocessing_2Polys(ex_num)
%
% % Input
%
% ex_num : (String) Example number


% Variable : Preprocessing 


close all; clc;

% Set upper and lower noise limit
emin = 1e-6;
emax = 1e-8;


% Constants
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';


% Perform GCD computation without preprocessing
mean_method = 'None';
bool_alpha_theta = false;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;

% Perform GCD computation with preprocessing
mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;



end