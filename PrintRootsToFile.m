
function [] = PrintRootsToFile(arr_RootFindingMethod, arr_ForwardErrors, arr_BackwardErrors)
% Print results of root finding computation to a text file
%
% % Inputs
%
% arr_RootFindingMethod : (Array of Strings) Array of the names of the root
% finding methods used
%
% arr_errors : (Array of Floats)


global SETTINGS


nMethods = length(arr_RootFindingMethod);

fullFileName = sprintf('Results/Results_o_roots.txt');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1 : 1 : nMethods
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        
        % Only print results if my method
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, method_ForwardError);
        end
        
        
    end
    
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    for i = 1 : 1 : nMethods
        
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        
        % Only print results if my method
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, method_ForwardError)
        end
        
    end
    fclose(fileID);
    
end



    function WriteNewLine(method_name, myBackwardError, myForwardError)
        % Write a new line of the text file
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            num2str(SETTINGS.BOOL_LOG),...
            SETTINGS.SYLVESTER_BUILD_METHOD,...
            SETTINGS.GCD_COEFFICIENT_METHOD,...
            method_name,...
            SETTINGS.DECONVOLUTION_METHOD_HX,...
            SETTINGS.DECONVOLUTION_METHOD_WX,...
            myForwardError,...
            myBackwardError...
            );
    end


    function WriteHeader()
        % If the file doesnt already exist, write a header to the text file
        fprintf(fileID,'DATE, EX_NUM, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD, LOW_RANK_ITE, APF_METHOD, APF_ITE, BOOL_LOG, SYLVESTER_BUILD_METHOD, GCD_METHODM, METHOD_NAME, DECONVOLUTION_METHOD HX, DECONVOLUTION_METHOD WX, FORWARD_ERROR, BACKWARD_ERROR \n');
    end


end
