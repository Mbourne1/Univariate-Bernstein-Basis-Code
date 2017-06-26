function [] = o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    Sylvester_Build_Method, rank_revealing_metric)
% o_gcd_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method, Sylvester_Build_Method)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% emin: (Float) Signal to noise ratio (minimum)
%
% emax: (Float) Signal to noise ratio (maximum)
%
% mean_method : (String) Method for taking mean of entires in S_{k}
%   * 'Geometric Mean Matlab Method'
%   * 'Geometric Mean My Method'
%
% bool_alpha_theta : (Bool) true or false if preprocessing is performed
%   * true
%   * false
%
% low_rank_approx_method : (String)
%   * 'Standard STLN'
%   * 'Standard SNTLN'
%   * 'Root Specific SNTLN'
%   * 'None'
%
% apf_method : (String)
%   * 'Standard APF NonLinear'
%   * 'Standard APF Linear'
%   * 'None'
%
% Sylvester_Build_Method : (String)
%   * 'T'
%   * 'DT'
%   * 'DTQ'
%   * 'TQ'
%   * 'DTQ Denominator Removed'
%   * 'DTQ Rearranged'
%
% rank_revealing_metric : (String)
%   * Singular Values :   
%   * Max/Min Singular Values :
%   * R1 Row Norms : 
%   * R1 Row Diagonals :
%   * Residuals :
%
% % Example
% >> o_gcd_Univariate_2Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'Minimum Singular Values')
% >> o_gcd_Univariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'DTQ', 'Minimum Singular Values')

global SETTINGS

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% % Set global variables.
SetGlobalVariables_GCD(...
    ex_num,...
    emin,...
    emax,...
    mean_method,...
    bool_alpha_theta,...
    low_rank_approx_method,...
    apf_method,...
    Sylvester_Build_Method, ...
    rank_revealing_metric);

% Print the parameters.
LineBreakLarge()
fprintf('INPUT VARIABLES\n')
fprintf('\t EXAMPLE NUMBER : %s \n',ex_num)
fprintf('\t EMIN : %s \n' , num2str(SETTINGS.EMIN))
fprintf('\t EMAX : %s \n' , num2str(SETTINGS.EMAX))
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA))
fprintf('\t RANK REVEALING METRIC : %s \n', SETTINGS.RANK_REVEALING_METRIC);
fprintf('\t LOW RANK APPROX METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD)
fprintf('\t LOG: %s \n', num2str(SETTINGS.BOOL_LOG))
fprintf('\t SYLVESTER BUILD METHOD: %s \n', SETTINGS.SYLVESTER_BUILD_METHOD)
LineBreakLarge()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get roots from example file
[fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);

% Get degree of f(x) and g(x)
m = GetDegree(fx_exact);
n = GetDegree(gx_exact);

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);
gx = AddVariableNoiseToPoly(gx_exact, emin, emax);

% Set upper and lower limit of the degree of the GCD
lowerLimit_t = 1;
upperLimit_t = min(m, n);
limits_t = [lowerLimit_t upperLimit_t];

% Set range of rank revealing metric
rank_range = [0 0];

% Obtain the coefficients of the GCD d and quotient polynomials u(x) and v(x).
[~, ~, dx_calc, ux_calc, vx_calc] = o_gcd_2Polys_mymethod(fx, gx, limits_t, rank_range);

t_calc = GetDegree(dx_calc);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.

LineBreakMedium();
try
    
    error.dx = GetError('d', dx_exact, dx_calc);
    error.ux = GetError('u', ux_exact, ux_calc);
    error.vx = GetError('v', vx_exact, vx_calc);
    
catch
    
    fprintf('Error : Can not compare computed GCD with Exact GCD. Check degree of GCD is computed correctly \n')
    error.dx = 1000;
    error.ux = 1000;
    error.vx = 1000;
    
end

% Print results to results file
PrintToFile(m, n, t_calc, error)
LineBreakMedium();

end





function [error] = GetError(name, fx_calc,fx_exact)
% Get distance between f(x) and the calulated f(x)

% Get the angle between the two vectors
% angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
% angle_error = 1 - angle;
% fprintf('\tCalculated angle error : %8.2e \n', angle_error)

fx_calc  = NormaliseVector(fx_calc);
fx_exact = NormaliseVector(fx_exact);

% Calculate relative errors in outputs
rel_error_fx = norm(abs(fx_calc - fx_exact) ./ fx_exact);

fprintf('\tCalculated relative error %s : %8.2e \n ',name, rel_error_fx);

error = norm(abs(fx_calc - fx_exact) );

fprintf('\tCalculated error %s : %8.2e \n', name,error);



end



function [] = PrintToFile(m, n, t, error)
% Print results of gcd computation to a text file
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% t : (Int) Computed degree of the GCD d(x)
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
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t),...
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
        fprintf(fileID,'DATE, EX_NUM, m, n, t, ERROR_UX, ERROR_VX, ERROR_DX, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,BOOL_LOG,SYLVESTER_BUILD_METHOD,GCD_METHOD, RANK_REVEALING_METRIC\n');
    end


end


