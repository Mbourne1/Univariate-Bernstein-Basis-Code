function [] = Experiment4LowRankApprox_2Polys(ex_num)

close all; clc;

% % Good Examples
%
% ex_num = '2'; el = 1e-8; eu = 1e-7;
%
% ex_num = '12' ; el = 1e-6 ; eu = 1e-5;
% ex_num = '15'; el = 1e-12; eu = 1e-8;

% 3 - STLN gives worse results
% 6 - worse results

el = 1e-8;
eu = 1e-7;



rank_revealing_metric = 'Minimum Singular Values';
%rank_revealing_metric = 'R1 Row Norms' ;

% R1 Row Norms, R1 Row Diagonals, Minimum Singular Values, Normalised Minimum Singular Values


%mean_method = 'None';
%bool_preproc = false;

apf_method = 'None';

sylvester_build_method = 'DTQ';

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;


arr_low_rank_approx_method = {'None', 'Standard STLN', 'Standard SNTLN'};


for i = 1 : 1 : length(arr_low_rank_approx_method)
    
    low_rank_approx_method = arr_low_rank_approx_method{i};
    o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_build_method, rank_revealing_metric) ;

end




end

