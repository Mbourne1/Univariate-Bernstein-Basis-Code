
function [] = PrintGCDToFile(m, n, t_exact, t_calc, error)
% Print results of gcd computation to a text file
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% t_exact : (Int) Computed degree of the GCD d(x)
%
% t_calc : (Int)
%
% error : [Float Float Float] Array of errors e
%   error.dx
%   error.ux
%   error.vx


global SETTINGS

%myDate = datetime('today');
%[year,month,day] = ymd(myDate);
%fullFileName = sprintf('Results/Results_o_gcd-%s-%s-%s.txt',num2str(year),num2str(month),num2str(day));

fullFileName = sprintf('Results/Results_o_gcd.txt');


% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end



    function WriteNewLine()
        % Write a new line of the text file
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t_exact),...
            int2str(t_calc),...
            num2str(error.ux),...
            num2str(error.vx),...
            num2str(error.dx),...
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
            SETTINGS.RANK_REVEALING_METRIC...
            );
    end


    function WriteHeader()
        % If the file doesnt already exist, write a header to the text file
        fprintf(fileID,'DATE, EX_NUM, m, n, t_exact, t_calc, ERROR_UX, ERROR_VX, ERROR_DX, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,BOOL_LOG,SYLVESTER_BUILD_METHOD,GCD_METHOD, RANK_REVEALING_METRIC\n');
    end


end