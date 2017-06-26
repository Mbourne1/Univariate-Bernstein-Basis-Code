
% % Experiment : Shows how the degree of the GCD can be computed from the
% singular values, but as noise increases, ability to compute degree of GCD
% decreases.

% Variable - Noise level 1e-12 -> 1e-4
% Variable - Include/Exclude Preprocessing

close all; 
clc;

ex_num = '9';


% Without Preprocessing --> Noise increases
boolPreproc = false;
meanMethod = 'None';

el = 1e-12;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-10;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-8;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-6;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-4;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

% With Preprocessing

boolPreproc = true;
meanMethod = 'Geometric Mean Matlab Method';
el = 1e-12;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-10;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-8;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-6;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')

el = 1e-4;
o_gcd_Univariate_2Polys(ex_num, el, el, meanMethod, boolPreproc, 'None', 'None', 'DTQ')