function [] = o_Deconvolution_Batch

arr_ex_num = {'1', '2', '3', '4', '5', '6','7','8'};
arr_noise = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12};
arr_bool_preproc = {true, false};

log_filename = 'log-deconvolution.txt';


parfor i1 = 1:1:length(arr_ex_num)
    ex_num = arr_ex_num{i1};
    
    for i2 = 1:1:length(arr_noise)
        emin = arr_noise{i2};
        
        for i3 = 1:1:length(arr_bool_preproc)
            bool_preproc = arr_bool_preproc{i3};
            
            
            try
                close all;
                clc;
                o_Deconvolution(ex_num, emin, bool_preproc)
                
                % Print success to log file
                fileId = fopen(log_filename, 'a');
                
                fprintf(fileId, '%s %s \n',...
                    datetime('now'),...
                    'Success'...
                );
                
                fclose(fileId);
                
            catch err
                % Print failure to log file
                fileId = fopen(log_filename, 'a');
                
               
                
                fprintf(fileId, '%s %s \n', ...
                    datetime('now'),...
                    getReport(err)...
                );
                
                fclose(fileId);
           end
        end
    end
end

    


end