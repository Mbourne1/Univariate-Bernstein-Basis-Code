function [] = o_roots_batch

ex_num_arr = {'Example Zeng'};
emin_arr = {1e-12};
mean_method_arr = {'None', 'Geometric Mean Matlab Method'};

emax = 1e-12;

bool_alpha_theta_arr = {'y','n'};
low_rank_approx_method_arr = {'None','Standard STLN'};
apf_method_arr = {'None'};

bool_q_arr = {'y','n'};
bool_log_arr = {'y','n'};

% Deconvolution method

DECONVOLVE_METHOD_arr = {'Separate', 'Batch'};

% the sequence h_{x} is obtained from
roots_hx_arr = {'From ux','From Deconvolutions'};

global SETTINGS


for i1 = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i1};
    
    for i2 = 1:1:length(emin_arr)
        
        emin = emin_arr{i2};
        
        for i3 = 1:1:length(mean_method_arr)
            
            mean_method = mean_method_arr{i3};
            
            for i4 = 1:1:length(bool_alpha_theta_arr)
                
                bool_alpha_theta = bool_alpha_theta_arr{i4};
                
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    
                    low_rank_approx_method = low_rank_approx_method_arr{i5};
                    
                    for i6 = 1:1:length(apf_method_arr)
                        
                        apf_method = apf_method_arr{i6};
                        
                        for i7 = 1:1:length(bool_q_arr)
                            
                            SETTINGS.BOOL_Q = bool_q_arr{i7};
                            
                            for i8 = 1:1:length(bool_log_arr)
                                
                                SETTINGS.BOOL_LOG = bool_log_arr{i8};
                                
                                for i9 = 1:1:length(roots_hx_arr)
                                    
                                    
                                    SETTINGS.ROOTS_HX = roots_hx_arr{i9};
                                    
                                    for i10 = 1:1:length(DECONVOLVE_METHOD_arr)
                                        
                                        SETTINGS.DECONVOLVE_METHOD = DECONVOLVE_METHOD_arr{i10};
                                        
                                        o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
                                        
                                    end
                                end
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
end