
close all; clc;

% Good examples
% 3 1e-7 1e-12

% Bad Examples
%
% 4 - Fails at final couple of gcd computations
% 5 - Correct multiplicity structure, roots are not precise, large backward
% error


% Queries
%
% 


ex_num = '10';
emin = 1e-10;
emax = 1e-12;
mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';


arrDeconvolution_method = {'Separate', 'Batch', 'Batch With STLN', 'Batch Constrained', 'Batch Constrained With STLN'};
arrDeconvolution_preproc = {true, false};

%deconvolution_method_hx = 'Separate';
deconvolution_method_wx = 'Batch';


for i = 1:1:length(arrDeconvolution_method)
    
    for i2 = 1:1: length(arrDeconvolution_preproc)
        
        deconvolution_method_hx = arrDeconvolution_method{i};
        deconvolution_preproc = arrDeconvolution_preproc{i2};
        
        o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
            low_rank_approx_method, apf_method, sylvester_build_method, ...
            rank_revealing_metric, deconvolution_method_hx, ...
            deconvolution_method_wx, deconvolution_preproc)
        
    
        
    end
end



