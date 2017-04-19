function [] = o_gcd_Univariate_2Polys_DegreeElevationTest(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, Sylvester_Build_Method, p, q)
% o_gcd_Univariate_2Polys_DegreeElevationTest(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, p , q)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex: (String) Example Number
%
% emin: (Float) Signal to noise ratio (minimum)
%
% emax: (Float) Signal to noise ratio (maximum)
%
% mean_method : (String) Method for taking mean of entires in S_{k}
%           'Geometric Mean Matlab Method'
%           'Geometric Mean My Method'
%
% bool_alpha_theta : (Boolean) true or false if preprocessing is performed
%
% low_rank_approx_method : (String)
%           'Standard STLN'
%           'Standard SNTLN'
%           'Root Specific SNTLN'
%           'None'
%
% apf_method : (String)
%           'Standard APF NonLinear'
%           'Standard APF Linear'
%           'None'
%
% Sylvester_Build_Method : (String)
%           'T'
%           'DT'
%           'DTQ'
%           'TQ'
%           'DTQ Rearranged Denom Removed'
%           'DTQ Rearranged'
%
% p : (Int) Degree elevation of polynomial f(x)
%
% q : (Int) Degree elevation of polynomial g(x)
%
% % Example
% >> o_gcd_Univariate_2Polys_DegreeElevationTest('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 2, 5)
% >> o_gcd_Univariate_2Polys_DegreeElevationTest('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear','DTQ', 4, 3)
%
% % Custom Example
%
% ex_num = 'Custom:m=10 n=10 t=5 low=0 high=1'
% >> o_gcd_Univariate_2Polys_DegreeElevationTest(ex_num,1e-10,1e-12,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ', 2, 5)

global SETTINGS

% Set the problem type to a GCD problem
problemType = 'GCD';

% Add subfolders
restoredefaultpath
addpath(...
    'Build Matrices',...
    'Formatting',...
    'GCD Methods',...
    'Get Cofactor Coefficients',...
    'Get GCD Coefficients',...
    'Measures',...
    'Plotting',...
    'Preprocessing',...
    'Results',...
    'Sylvester Matrix')

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Low Rank Approximation'));


% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% % Set global variables.
SetGlobalVariables(problemType,...
    ex_num,...
    emin,...
    emax,...
    mean_method,...
    bool_alpha_theta,...
    low_rank_approx_method,...
    apf_method,...
    Sylvester_Build_Method);

% Print the parameters.
LineBreakLarge()
fprintf('INPUT VARIABLES\n')
fprintf('\t EXAMPLE NUMBER : %s \n',ex_num)
fprintf('\t EMIN : %s \n' , num2str(SETTINGS.EMIN))
fprintf('\t EMAX : %s \n' , num2str(SETTINGS.EMAX))
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n', SETTINGS.BOOL_ALPHA_THETA)
fprintf('\t LOW RANK APPROX METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD)
fprintf('\t LOG: %s \n', SETTINGS.BOOL_LOG)
fprintf('\t SYLVESTER BUILD METHOD: %s \n', SETTINGS.SYLVESTER_BUILD_METHOD)
LineBreakLarge()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get roots from example file
[fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);

m = GetDegree(fx_exact);
n = GetDegree(gx_exact);

max_r = min(p, q);

fx_star_exact = Bernstein_DegreeElevate_Univariate(fx_exact, p);
gx_star_exact = Bernstein_DegreeElevate_Univariate(gx_exact, q);

m_star = m + p;
n_star = n + q;

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx_star = AddVariableNoiseToPoly(fx_star_exact, emin, emax);
gx_star = AddVariableNoiseToPoly(gx_star_exact, emin, emax);

% set upper and lower limit of the degree of the GCD
lower_lim = 1;
upper_lim = min(m_star, n_star);
limits = [lower_lim upper_lim];

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~, ~, dx_calc, ux_calc, vx_calc] = o_gcd_2Polys_mymethod(fx_star, gx_star, limits);

t_star_calc = GetDegree(dx_calc);
t = t_star_calc - max_r;

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
PrintToFile(m, n, t_star_calc, error)
LineBreakMedium();

end

function [error] = GetError(name,f_calc,f_exact)
% Get distance between f(x) and the calulated f(x)

% Get the angle between the two vectors
% angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
% angle_error = 1 - angle;
% fprintf('\tCalculated angle error : %8.2e \n', angle_error)

f_calc  = Normalise(f_calc);
f_exact = Normalise(f_exact);

% Calculate relative errors in outputs
rel_error_f = norm(abs(f_calc - f_exact) ./ f_exact);

fprintf('\tCalculated relative error %s : %8.2e \n ',name,rel_error_f);

error = norm(abs(f_calc - f_exact) );

fprintf('\tCalculated error %s : %8.2e \n', name,error);



end



function [] = PrintToFile(m,n,t,error)
% Print results of gcd computation to a text file
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Computed degree of the GCD d(x)
%
% error : array of errors e
%   error.dx
%   error.ux
%   error.vx


global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd%s.txt',datetime('today'));

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
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
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
            SETTINGS.GCD_COEFFICIENT_METHOD...
            );
    end


    function WriteHeader()
        % If the file doesnt already exist, write a header to the text file
        fprintf(fileID,'DATE,EX_NUM,m,n,t,ERROR_UX,ERROR_VX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,BOOL_LOG,SYLVESTER_BUILD_METHOD,GCD_METHOD\n');
    end


end


