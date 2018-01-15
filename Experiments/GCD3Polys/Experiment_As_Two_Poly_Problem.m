function [] = Experiment_As_Two_Poly_Problem(ex_number, bool_preproc)

% Variable : Preprocessing

% Good Examples

% Example 11 - 1e-10 , 1e-12



close all; clc;

% Set max and min noise level
emin = 1e-9;
emax = 1e-9;

% set preprocessing related variables
switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta= false;
end

% Constants
apf_method = 'None';

% % Low Rank Approximation Method
% 'None'
% 'Standard STLN'
% 'Standard SNTLN'
low_rank_approx_method = 'None';

% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_build_method = 'DTQ';

arr_ex_variant = {'a','b','c'};

for i = 1 : 1  : length(arr_ex_variant)
    
    ex_variant = arr_ex_variant{i};
    ex_num = strcat(ex_number, ex_variant);
    o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;
    
    
end






three_poly_problem = false;

if three_poly_problem == true
    
    % Three Polynomials - 2 x 3 subresutlant - f appears twice S(f,g,h)
    
    % nEquations refers to the choice of a (2x3) or a (3x3) partitioned
    % subresultant matrix.
    nEquations = '2';
    
    for i = 1 : 1 : length(arr_ex_variant)
        ex_variant = arr_ex_variant{i};
        ex_num = strcat(ex_number, ex_variant);
        o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, nEquations, rank_revealing_metric)
    end
    
    
    
    % Three polynomials - 3 x 3 subresultant
    
    nEquations = '3';
    o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, nEquations, rank_revealing_metric)
    
end




end