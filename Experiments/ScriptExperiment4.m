close all; clc;

ex_num = '7';

el = 1e-10;

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;

o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ', 'Minimum Singular Values') ;



mean_method = 'None';
bool_preproc = false;

o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ', 'Minimum Singular Values') ;
