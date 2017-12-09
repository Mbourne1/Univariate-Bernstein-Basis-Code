close all; clc; 


ex_num = '11';
el = 0;
eu = 0;

%mean_method = 'Geometric Mean Matlab Method';
%bool_preproc = true;

mean_method = 'None';
bool_preproc = false;

sylvester_build_method = 'T';
rank_revealing_metric = 'Minimum Singular Values';

low_rank_approx_method = 'None';

apf_method = 'Last Non-zero Row';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;



