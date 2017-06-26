% Comparing the types of Sylvester subresultant matrix in the computation
% of the degree of the GCD.
%
% Variable : Preprocessing
% Variable : Sylvester Build Type



close all; clc;

ex_num = '3';

el = 1e-10;

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;

o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DT') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'TQ') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ Denominator Removed')





mean_method = 'Geometric Mean Matlab Method';
bool_preproc = false;


o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DT') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'TQ') ;
o_gcd_Univariate_2Polys(ex_num, el, 1e-12, mean_method, bool_preproc, 'None', 'None', 'DTQ Denominator Removed')
