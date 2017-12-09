function [] = Experiment2Preprocessing_2Polys(ex_num)

% Variable : Preprocessing 

% Good Examples

% Example 11 - 1e-10 , 1e-12


close all; clc;

emin = 1e-6;
emax = 1e-8;


% Constants
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';



mean_method = 'None';
bool_alpha_theta = false;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;


mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;



end