function [] = Experiment4LowRankApprox_2Polys(ex_num)
% Experiment with optimal low rank approximation method for the t-th
% subresultant matrix. Where 't' is the degree of the GCD.
%
% % Inputs
%
% ex_num : (String) Example number

close all;
clc;


% Set upper and lower noise level
el = 1e-8;
eu = 1e-7;


% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% Set preprocessing related variables
bool_preprocessing = true;
if (bool_preprocessing == true)
    mean_method = 'Geometric Mean Matlab Method';
    bool_preproc = true;
    
else
    mean_method = 'None';
    bool_preproc = false;
    
end

% APF method
apf_method = 'None';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_build_method = 'DTQ';


% Get array of low rank approximation methods
arr_low_rank_approx_method = {'None', 'Standard STLN', 'Standard SNTLN'};

% For each low rank approximation method compute the GCD
for i = 1 : 1 : length(arr_low_rank_approx_method)
    
    low_rank_approx_method = arr_low_rank_approx_method{i};
    
    o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;
    
end




end

