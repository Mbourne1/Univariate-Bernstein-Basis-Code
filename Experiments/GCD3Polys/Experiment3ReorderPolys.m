close all; clc;

% ex_num (String) : This string has the form a a number followed by the
% letter 'a' 'b' or 'c', for the three different polynomial orderings.
ex_num = '11b';

% emin : (Float
% emax : (Float)
emin = 1e-12;
emax = 1e-10;

% arrEx_num = {'11'}; emin = 1e-9; emax = 1e-9;

% Good Examples
% 1 2 3 4 5 6 7 8 9 10
% 11a - 
% 11b

% Notes
% example 1

%mean_method = 'Geometric Mean Matlab Method';
%bool_alpha_theta = true;

mean_method = 'None';
bool_alpha_theta = false;

low_rank_approx_method = 'None';
apf_method = 'None';

nEquations = '2';


rank_revealing_metric = 'Minimum Singular Values';
sylvester_build_method = 'DTQ';


o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_build_method, nEquations, rank_revealing_metric)

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;


o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_build_method, nEquations, rank_revealing_metric)
