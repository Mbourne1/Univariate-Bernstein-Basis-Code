

close all; clc;




ex_num = '1';
arrEmax = {1e-10,1e-8,1e-6};
emin = 1e-12;

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

low_rank_approx_method = 'None';
apf_method = 'None';

sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';


arrDeconvolution_method = {'Separate', 'Batch', 'Batch With STLN', 'Batch Constrained', 'Batch Constrained With STLN'};
arrDeconvolution_preproc = {true, false};

arrDeconvolution_method_wx = {'Batch', 'Separate'};


parfor i = 1:1:length(arrDeconvolution_method)
    for i2 = 1:1: length(arrDeconvolution_preproc)
        for i3 = 1: 1: length(arrDeconvolution_method_wx)
            for i4 = 1 : 1 : length(arrEmax)
                
                emax = arrEmax{i4};
                deconvolution_method_hx = arrDeconvolution_method{i};
                deconvolution_preproc = arrDeconvolution_preproc{i2};
                deconvolution_method_wx = arrDeconvolution_method_wx{i3};
                
                o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
                    low_rank_approx_method, apf_method, sylvester_build_method, ...
                    rank_revealing_metric, deconvolution_method_hx, ...
                    deconvolution_method_wx, deconvolution_preproc)
                
            end
        end
    end
end