function [] = Experiment_LowRankApprox_DeconvolutionMethods(ex_num, bool_preproc)
%
%
%
% >> Experiment_LowRankApprox_DeconvolutionMethods('10')

close all; clc;


%
% -------------------------------------------------------------------------
% Constants

emin = 1e-4;
emax = 1e-4;


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

%
% -------------------------------------------------------------------------
% Irrelevant - only interested in what happens before deconvolution stages





%low_rank_approx_method = 'Standard SNTLN';

%
% -------------------------------------------------------------------------
% Variables


deconvolution_method_wx = 'Batch';
bool_deconvolution_preproc = true;

%arr_deconvolution_method_hx = {'Separate', 'Batch', 'Batch With STLN', 'Batch Constrained', 'Batch Constrained With STLN'};
arr_deconvolution_method_hx = {'Separate', 'Batch'};

arr_low_rank_approx_method = {'None', 'Standard SNTLN'};



for i1 = 1 : 1 : length(arr_deconvolution_method_hx)
    
    deconvolution_method_hx = arr_deconvolution_method_hx{i1};
    
    for i2 = 1 : 1 : length(arr_low_rank_approx_method)
        
        low_rank_approx_method = arr_low_rank_approx_method{i2};
        
        
        try
        [epsilon_fi(i1,i2), epsilon_hi(i1,i2), epsilon_wi(i1,i2)] = o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
            low_rank_approx_method, apf_method, sylvester_build_method, ...
            rank_revealing_metric, deconvolution_method_hx, ...
            deconvolution_method_wx, bool_deconvolution_preproc)
        
        
        catch
        
        end
        
    end
    
end

display(epsilon_fi)

display(epsilon_hi)

display(epsilon_wi)

sameaxes()


end


