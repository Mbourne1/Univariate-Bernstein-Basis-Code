function [] = o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,Sylvester_Build_Method)
% o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex:   Example Number
%
% emin: Signal to noise ratio (minimum)
%
% emax: Signal to noise ratio (maximum)
%
% mean_method : Method for taking mean of entires in S_{k}
%
%           'Geometric Mean Matlab Method'
%           'Geometric Mean My Method'
%
% bool_alpha_theta : 'y' or 'n' if preprocessing is performed
%
% low_rank_approx_method :
%           'Standard STLN'
%           'Standard SNTLN'
%           'Root Specific SNTLN'
%           'None'
%
% apf_method :
%           'Standard APF NonLinear'
%           'Standard APF Linear'
%           'None'
%
% Sylvester_Build_Method : 
%           'T'
%           'DT'
%           'DTQ'
%           'TQ'
%           'DTQ Rearranged Denom Removed'
%           'DTQ Rearranged'
%
% % Example
% >> o_gcd('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard STLN','Standard APF NonLinear','DTQ')
% 
% % Custom Example
%
% ex_num = 'Custom:m=10 n=10 t=5 low=0 high=1'
% >> o_gcd(ex_num,1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ')


% Set the problem type to a GCD problem
problemType = 'GCD';

restoredefaultpath
addpath(...
    'Build Matrices',...
    'Formatting',...
    'Measures',...
    'Plotting',...
    'Preprocessing',...
    'Results',...
    'Sylvester Matrix')

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approx'));

% Consistency of input parameters.

% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
% if emin > emax
%     fprintf('minimum noise greater than maximum noise \n swapping values...\n')
%     [emin,emax] = swap(emin,emax);
% end

% Set the global variables
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
PrintGlobalVariables()
global SETTINGS

fprintf('PARAMETERS:\n\n')
fprintf('\tmin noise : %2.4e \n\tmax noise : %2.4e \n',emin,emax)
fprintf('INPUT VARIABLES\n')

LineBreakLarge()
fprintf('\t EXAMPLE NUMBER : %s \n',ex_num)
fprintf('\t MEAN METHOD : %s \n',SETTINGS.MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n',SETTINGS.BOOL_ALPHA_THETA)
fprintf('\t Low Rank Approximation Method : %s \n',SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF Method : %s \n ',SETTINGS.APF_METHOD)
fprintf('\t LOG: %s \n',SETTINGS.BOOL_LOG)
fprintf('\t Sylvester Build Method : %s \n',SETTINGS.SYLVESTER_BUILD_METHOD)
fprintf('')
LineBreakLarge()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get roots from example file
[f_exact, g_exact,d_exact] = Examples_GCD(ex_num);

% Print the coefficients of f(x) and g(x)
%PrintCoefficients_Bivariate_Bernstein(f_exact,'f')
%PrintCoefficients_Bivariate_Bernstein(g_exact,'g')

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = AddVariableNoiseToPoly(f_exact,emin,emax);
gx = AddVariableNoiseToPoly(g_exact,emin,emax);

% set upper and lower limit of the degree of the GCD
lower_lim = 1;
upper_lim = min(GetDegree(fx),GetDegree(gx));

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~,~,dx_calc,ux_calc,vx_calc] = o_gcd_mymethod(fx,gx,[lower_lim,upper_lim]);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.

LineBreakMedium();
try
error_dx = GetError('d',d_exact,dx_calc);
catch
    error_dx = 1000;
end
PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error_dx)
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



function [] = PrintToFile(m,n,t,error_dx)


global SETTINGS

fullFileName = 'Results/Results_o_gcd.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        num2str(error_dx),...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.APF_METHOD,...
        SETTINGS.BOOL_LOG,...
        SETTINGS.SYLVESTER_BUILD_METHOD,...
        SETTINGS.GCD_COEFFICIENT_METHOD...
    );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end


