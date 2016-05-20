function [] = o_gcd_batch()


ex_num_arr = {'1'};

bool_alpha_theta_arr = {'y','n'};

emin_arr = {1e-10,1e-11,1e-12};

low_rank_approx_method_arr = {'Standard STLN','None'};

apf_method_arr = {'None', 'Standard APF'};

mean_method_arr = {'None','Geometric Mean Matlab Method'};

bool_q_arr = {'y','n'};

bool_log_arr = {'y','n'};

global SETTINGS

% Changing example number
for i1 = 1:1:length(ex_num_arr)
    ex_num = ex_num_arr{i1};
    
    % Changing lower noise boundary
    for i2 = 1:1:length(emin_arr)
        
        emin = emin_arr{i2};
        emax = 1e-12;
        
        % Changing low rank approx method
        for i3 = 1:1:length(low_rank_approx_method_arr)
            
            low_rank_approx_method = low_rank_approx_method_arr{i3};
            
            % Changing alpha theta boolean
            for i4 = 1:1:length(bool_alpha_theta_arr)
                
                bool_alpha_theta = bool_alpha_theta_arr{i4};
                
                for i5 = 1:1:length(mean_method_arr)
                    
                    mean_method = mean_method_arr{i5};
                    
                    for i6 = 1:1:length(bool_q_arr)
                        SETTINGS.BOOL_Q = bool_q_arr{i6};
                        
                        for i7 = 1:1:length(bool_log_arr)
                            SETTINGS.BOOL_LOG = bool_log_arr{i7};
                            
                            apf_method = 'None';
                            
                            try
                                o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
                            catch
                                fprintf([mfilename ' : ' 'Error'])
                            end
                        end
                        
                        
                    end
                    
                    
                end
                
                
            end
        end
    end
    
end



