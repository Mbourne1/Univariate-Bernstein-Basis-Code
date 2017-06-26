function [] = o_gcd_Univariate_3Polys_batch()
% Perform a batch of polynomial GCD computations
%
% >> o_gcd_3Polys_batch
%
%

ex_num_arr = {'1'};
%ex_num_arr = {  '2'};
bool_alpha_theta_arr = {true, false};
emin_arr = {1e-08,1e-10,1e-12};
%emin_arr = {1e-08,1e-10};
low_rank_approx_method_arr = {'Standard SNTLN', 'Standard STLN', 'None'};
apf_method_arr = {'Standard APF Nonlinear', 'None'};
mean_method_arr = {'Geometric Mean Matlab Method','None'};
bool_log_arr = {false};
gcd_coefficient_method_arr = {'ux and vx'};
%Sylvester_Build_Method_arr = {'T','DT','DTQ','TQ','DTQ Denominator Removed','DTQ Rearranged'};
%Sylvester_Build_Method_arr = {'DTQ'};

Sylvester_Build_Method_arr = {'T','DT','DTQ','TQ'};


global SETTINGS



for i1 = 1:1:length(bool_log_arr)
    
    SETTINGS.BOOL_LOG = bool_log_arr{i1};
    
    for i2 = 1:1:length(gcd_coefficient_method_arr)
        
        SETTINGS.GCD_COEFFICIENT_METHOD = gcd_coefficient_method_arr{i2};
        
        % Changing example number
        parfor i3 = 1:1:length(ex_num_arr)
            
            ex_num = ex_num_arr{i3};
            
            % Changing lower noise boundary
            for i4 = 1:1:length(emin_arr)
                
                emin = emin_arr{i4};
                emax = 1e-12;
                
                % Changing low rank approx method
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    
                    low_rank_approx_method = low_rank_approx_method_arr{i5};
                    
                    % Changing alpha theta boolean
                    for i6 = 1:1:length(bool_alpha_theta_arr)
                        
                        bool_alpha_theta = bool_alpha_theta_arr{i6};
                        
                        for i7 = 1:1:length(mean_method_arr)
                            mean_method = mean_method_arr{i7};
                            
                            for i8 = 1:1:length(Sylvester_Build_Method_arr)
                                
                                sylvester_build_method = Sylvester_Build_Method_arr{i8};
                                
                                for i9 = 1:1:length(apf_method_arr)
                                    
                                    apf_method = apf_method_arr{i9};
                                    try
                                        
                                        close all;
                                        clc;
                                        o_gcd_Univariate_3Polys(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,sylvester_build_method)
                                        fileId = fopen('log_GCD_3Polys.txt','a')
                                        fprintf(fileId,'%s \n',datetime('now') , 'Success');
                                        fclose(fileId);
                                        
                                    catch err
                                        
                                        fileId = fopen('log_GCD_3Polys.txt','a')
                                        fprintf(fileId,'%s %s \n\n\n',datetime('now'), getReport(err));
                                        fclose(fileId);
                                        
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



