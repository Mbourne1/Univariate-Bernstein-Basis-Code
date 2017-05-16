
close all; 
clc;

ex_num = '8';
emax = 1e-8;
preproc = 'n';

o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','DTQ');
o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','TQ');
o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','DT');
o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','T');
o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','DTQ Rearranged');
o_gcd_2Polys(ex_num,emax,1e-12,'None',preproc,'None','None','DTQ Rearranged Denom Removed');

LineBreakLarge()
fprintf('End of test \n')
LineBreakLarge()