close all; clc;

ex_num = '9';
el = 1e-10;
eu = 1e-12;


mean_method = 'None';
bool_preproc = false;



rank_revealing_metric = 'Minimum Singular Values';

sylvesterFormat = 'DTQ';


o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, 'None', 'None', sylvesterFormat, rank_revealing_metric);