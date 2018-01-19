function [] = Experiment3ReorderPolys(ex_num, bool_preproc)
%
% % Inputs
%
% ex_num : (String) Example Number 
% 
% bool_preproc : (Boolean) 
%

close all; clc;


emin = 1e-12;
emax = 1e-10;

% arrEx_num = {'11'}; emin = 1e-9; emax = 1e-9;



switch bool_preproc
    case 1
    mean_method = 'Geometric Mean Matlab Method';
    bool_alpha_theta = true;
    case 0
        
    mean_method = 'None';
    bool_alpha_theta = false;
end

% low_rank_approx_method
low_rank_approx_method = 'None';
apf_method = 'None';

nEquations = '2';


rank_revealing_metric = 'Minimum Singular Values';
sylvester_matrix_variant = 'DTQ';


o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;


o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)
