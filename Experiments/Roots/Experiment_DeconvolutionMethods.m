function [] = Experiment_DeconvolutionMethods(ex_num)
%
%
%
% >> Experiment_DeconvolutionMethods('1')

close all; clc;


%
% -------------------------------------------------------------------------
% Constants

emin = 1e-10;
emax = 1e-10;

%mean_method = 'None';
%bool_alpha_theta = false;

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';




%low_rank_approx_method = 'Standard SNTLN';
low_rank_approx_method = 'None';
%low_rank_approx_method = 'Standard STLN';


deconvolution_method_wx = 'Batch';

%
% -------------------------------------------------------------------------
% Variables 



bool_deconvolution_preproc = false;



deconvolution_method_hx = 'Separate';
        
o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


deconvolution_method_hx = 'Batch';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


deconvolution_method_hx = 'Batch With STLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


deconvolution_method_hx = 'Batch Constrained';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)

deconvolution_method_hx = 'Batch Constrained With STLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)











bool_deconvolution_preproc = true;



deconvolution_method_hx = 'Batch';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


deconvolution_method_hx = 'Batch With STLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


deconvolution_method_hx = 'Batch Constrained';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)

deconvolution_method_hx = 'Batch Constrained With STLN';

o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, bool_deconvolution_preproc)


sameaxes()

end


