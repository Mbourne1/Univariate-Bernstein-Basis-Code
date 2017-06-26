% This example, we want to see the computation of the degree of the GCD
% with different forms of the sylvester subresultant matrix. 
%
% Variable : Sylvester Build Type



close all; clc; 
 
 ex_num = '8'; 
 
 el = 1e-8; 
 eu = 1e-12;

 
 mean_method = 'None';
 bool_preproc = false;
 
 %mean_method = 'Geometric Mean Matlab Method';
 %bool_preproc = true;
 
 o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, 'None', 'None', 'DTQ');
 o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, 'None', 'None', 'DT');
 o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, 'None', 'None', 'TQ');
 o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, 'None', 'None', 'DTQ Denominator Removed');
 
