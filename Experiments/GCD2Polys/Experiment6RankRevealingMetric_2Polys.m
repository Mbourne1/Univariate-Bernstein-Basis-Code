close all; 
clc; 

% % Good Examples


ex_num = '10';
el = 1e-8;
eu = 1e-6;

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;

%mean_method = 'None';
%bool_preproc = false;

sylvester_build_method = 'DTQ';
low_rank_approx_method = 'None';
apf_method = 'None';


rank_revealing_metric = 'Minimum Singular Values';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;

rank_revealing_metric = 'Normalised Minimum Singular Values';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;

rank_revealing_metric = 'R1 Row Diagonals';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;

rank_revealing_metric = 'R1 Row Norms';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;