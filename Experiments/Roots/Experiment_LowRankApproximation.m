function [] = Experiment_LowRankApproximation(ex_num, bool_preproc)
%
%
%
% >> Experiment_LowRankApproximation('1')
% >> Experiment_LowRankApproximation('15')
close all; clc;


%
% -------------------------------------------------------------------------
% Constants

emin = 1e-10;
emax = 1e-8;

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end


apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';




%deconvolution_method_hx = 'Batch';
deconvolution_method_hx = 'Batch Constrained';


deconvolution_method_wx = 'Batch';

bool_deconvolution_preproc = true;









low_rank_approx_method = 'None';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


low_rank_approx_method = 'Standard STLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)



low_rank_approx_method = 'Standard SNTLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


end


