function [] = o_roots_batch

% count = 1;
% for i = 3:1:10
%     ex_num_arr{count} = sprintf('Custom:m=%i low=-10 high=10',i);
%     count = count + 1;
% end

ex_num_arr = ...
    {
    '1','2','3'
    %      '4'
    %      '5'
    %      '6'
    %      '7'
    %      '8'
    %      '9'
    %      '10'
    %      '11'
    %      '12'
    %      '13'
    %      'Example Zeng'
    };


emin_arr = ...
    {
    1e-12,...
    %1e-11,...
    % 1e-10,...
    % 1e-9,...
    };

mean_method_arr = ...
    {...
    'None',...
    'Geometric Mean Matlab Method'...
    };

emax = 1e-12;

bool_alpha_theta_arr = ...
    {...
    'y',...
    'n'...
    };

low_rank_approx_method_arr = ...
    {...
    'None',...
    'Standard STLN'...
    'Standard SNTLN'
    };

apf_method_arr = {'None'};

arr_sylvester_build_method =...
    {
        'DTQ',...
        'TQ',...
        'DT',...
        'T'
    };

bool_log_arr = ...
    {
    true,...
    false
    };

% Deconvolution method

DECONVOLVE_METHOD_arr = ...
    {
    'Separate',...
    'Batch'
    };

% the sequence h_{x} is obtained from
roots_hx_arr = ...
    {
    'From ux',...
    'From Deconvolutions'
    };


var_arr = [];
count = 1;
for i1 = 1:1:length(ex_num_arr)
    for i2 = 1:1:length(emin_arr)
        for i3 = 1:1:length(mean_method_arr)
            for i4 = 1:1:length(bool_alpha_theta_arr)
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    for i6 = 1:1:length(apf_method_arr)
                        var_arr(count,:) = [i1,i2,i3,i4,i5,i6];
                        count = count + 1;
                    end
                end
            end
        end
    end
end


global SETTINGS

SETTINGS.GCD_COEFFICIENT_METHOD = 'ux';





for i8 = 1:1:length(bool_log_arr)
    
    SETTINGS.BOOL_LOG = bool_log_arr{i8};
    
    for i9 = 1:1:length(roots_hx_arr)
        
        SETTINGS.ROOTS_HX = roots_hx_arr{i9};
        
        for i10 = 1:1:length(DECONVOLVE_METHOD_arr)
            
            SETTINGS.DECONVOLVE_METHOD = DECONVOLVE_METHOD_arr{i10};
            
            for i=1:1:size(var_arr,1)
                
                var = var_arr(i,:) ;
                ex_num = ex_num_arr{var(1)};
                emin = emin_arr{var(2)};
                mean_method = mean_method_arr{var(3)};
                bool_alpha_theta = bool_alpha_theta_arr{var(4)};
                low_rank_approx_method = low_rank_approx_method_arr{var(5)};
                apf_method = apf_method_arr{var(6)};
                
                try
                    o_roots(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
                catch err
                    fprintf(err.message);
                end
                
            end
            
            
            
        end
    end
end





end