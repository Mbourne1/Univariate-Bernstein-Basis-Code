function [] = Experiment6_LowRankApprox_3Polys(ex_num)

close all; clc; 

emin = 1e-10;
emax = 1e-8;

nEquations = '2';


apf_method = 'None';

sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';



global SETTINGS
 
SETTINGS.SCALING_METHOD = 'lambda_mu_rho';

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;
%mean_method = 'None';
%bool_alpha_theta = false;

low_rank_approx_method = 'None';

o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_build_method, nEquations, ...
    rank_revealing_metric)



low_rank_approx_method = 'Standard STLN';

o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_build_method, nEquations, ...
    rank_revealing_metric)


end


