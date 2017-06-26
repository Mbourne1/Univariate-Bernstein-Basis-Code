
% Testing the different forms of the Sylvester subresultant matrices for
% the computation of the degree and coefficients of the gcd 

close all; 
clc;

ex_num = '8';
emax = 1e-8;
boolPreproc = false;
meanMethod = 'None';

o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'DTQ');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'TQ');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'DT');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'T');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'DTQ Rearranged');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc, 'None', 'None', 'DTQ Denominator Removed');


% With Preprocessing
boolPreproc = true;
meanMethod = 'Geometric Mean Matlab Method';

o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','DTQ');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','TQ');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','DT');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','T');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','DTQ Rearranged');
o_gcd_Univariate_2Polys(ex_num, emax, 1e-12, meanMethod, boolPreproc,'None','None','DTQ Denominator Removed');



LineBreakLarge()
fprintf('End of test \n')
LineBreakLarge()