close all; clc; 


ex_num = '2';
el = 1e-10;
eu = 1e-8;

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;
sylvester_build_method = 'DTQ';

low_rank_approx_method = 'None';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, 'None', sylvester_build_method, 'Minimum Singular Values') ;

low_rank_approx_method = 'Standard STLN';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, 'None', sylvester_build_method, 'Minimum Singular Values') ;

low_rank_approx_method = 'Standard SNTLN';
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, 'None', sylvester_build_method, 'Minimum Singular Values') ;
