function [] = o_gcd_batch()

ex_num_arr = ...
    {
    '1',...
    '2',...
    %'Custom:m=15 n=9.t=7 low=-1 high=1'
    };

% 
% % Build a custom array of examples
% count = 1;
% for i = 3:1:14
%     for j = 3:1:7
%         for k = 1:1:min(i,j)
%             %
%             m = i;
%             n = j;
%             t = k;
%             
%             ex_num_arr{count} = sprintf('Custom:m=%i n=%i t=%i low=-1 high=1',m,n,t);
%             
%             count = count + 1;
%             
%         end
%         
%         
%     end
% end


bool_alpha_theta_arr = {'y','n'};

emin_arr = ...
    {
    1e-10,...
    1e-12...
    };


low_rank_approx_method_arr = ...
    {
    'Standard STLN',...
    'None'...
    };

apf_method_arr = ...
    {
    'Standard APF',...
    'None'
    };

mean_method_arr = ...
    {
    'Geometric Mean Matlab Method',...
    'None'...
    };

bool_q_arr = ...
    {
    'y',...
    'n'...
    };

bool_log_arr = ...
    {
    'y',...
    'n'
    };

gcd_coefficient_method_arr = ...
    {
    'ux',...
    'ux and vx'
    };





global SETTINGS

var_arr = [];

% Changing example number
for i1 = 1:1:length(ex_num_arr)
    for i2 = 1:1:length(emin_arr)
        for i3 = 1:1:length(low_rank_approx_method_arr)
            for i4 = 1:1:length(bool_alpha_theta_arr)
                for i5 = 1:1:length(mean_method_arr)
                    var_arr(count,:) = [i1,i2,i3,i4,i5];
                    count = count + 1;
                end
            end
        end
    end
end

for i6 = 1:1:length(bool_q_arr)
    SETTINGS.BOOL_Q = bool_q_arr{i6};
    for i7 = 1:1:length(bool_log_arr)
        SETTINGS.BOOL_LOG = bool_log_arr{i7};
        apf_method = 'None';
        for i8 = 1:1:length(gcd_coefficient_method_arr)
            SETTINGS.GCD_COEFFICIENT_METHOD = gcd_coefficient_method_arr{i8};
            
            parfor i = 1:1:length(var_arr)
                
                var = var_arr(i,:);
                
                ex_num = ex_num_arr{var(1)};
                emin = emin_arr{var(2)};
                emax = 1e-12;
                low_rank_approx_method = low_rank_approx_method_arr{var(3)};
                bool_alpha_theta = bool_alpha_theta_arr{var(4)};
                mean_method = mean_method_arr{var(5)};

                try
                    o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
                catch err
                    fprintf([mfilename ' : ' 'Error'])
                    display(err.message)
                end

                
                
            end
            
        end
    end
    
end



end



