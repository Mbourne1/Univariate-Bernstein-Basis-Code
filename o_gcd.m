function [] = o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
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
%           'Standard APF'
%           'None'
%
% % Example
% >> o_gcd('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard STLN','None')
% 
% % Custom Example
%
% ex_num = 'Custom:m=10 n=10 t=5 low=0 high=1'
% >> o_gcd(ex_num,1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard STLN','None')


% Set the problem type to a GCD problem
problemType = 'GCD';
addpath(...
    'BasisConversion',...
    'BernsteinMethods',...
    'Examples',...
    'Formatting',...
    'GCDDegree',...
    'GetCofactorCoefficients',...
    'GetGCDCoefficients',...
    'Low Rank Approx',...
    'MatrixBuilding',...
    'Plotting',...
    'Preprocessing',...
    'Results',...
    'APF')

% Consistency of input parameters.

% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
% if emin > emax
%     fprintf('minimum noise greater than maximum noise \n swapping values...\n')
%     [emin,emax] = swap(emin,emax);
% end

% Set the global variables
global SETTINGS
if isempty(SETTINGS)
    fprintf('Set Q and log')
    SETTINGS.BOOL_Q = 'y';
    SETTINGS.BOOL_LOG = 'n';
    SETTINGS.GCD_COEFFICIENT_METHOD = 'ux';
end

SetGlobalVariables(problemType,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)

% Print the parameters.
PrintGlobalVariables()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails

% Add neccesary paths.
addpath 'Measures'

% Get roots from example file
[f_exact, g_exact,d_exact] = Examples_GCD(ex_num);

display(f_exact)
display(g_exact)
display(d_exact)

% Print the coefficients of f(x) and g(x)
%PrintCoefficients_Bivariate_Bernstein(f_exact,'f')
%PrintCoefficients_Bivariate_Bernstein(g_exact,'g')

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = VariableNoise(f_exact,emin,emax);
gx = VariableNoise(g_exact,emin,emax);

% set upper and lower limit of the degree of the GCD
lower_lim = 1;
upper_lim = min(GetDegree(fx),GetDegree(gx));

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~,~,dx_calc,ux_calc,vx_calc] = o_gcd_mymethod(fx,gx,[lower_lim,upper_lim]);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.
%PrintCoefficients('u',u_exact, ux_calc);
%PrintCoefficients('v',v_exact, vx_calc);
PrintCoefficients('d',d_exact, dx_calc);

error_dx = GetError('d',d_exact,dx_calc);

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error_dx)


end


function [] = PrintRoots(f,f_roots)

% print out the exact roots of f,g and d
fprintf('\nRoots of %s \n',f);
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[f_roots(:,1),f_roots(:,2)]');
fprintf('\n');

end

function [] = PrintCoefficients(name,u_exact,u_calc)

% Normalise quotient polynomial u
u_calc  = Normalise(u_calc);
u_exact = Normalise(u_exact);

fprintf('\nCoefficients of %s \n\n',name);
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(u_exact(:,1))';  real(u_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n', mat);
fprintf('\n');
GetError(name,u_calc,u_exact);

end

function [error] = GetError(name,f_calc,f_exact)
% Get distance between f(x) and the calulated f(x)

% Get the angle between the two vectors
angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
angle_error = 1 - angle;
fprintf('Calculated angle error : %8.2e \n', angle_error)

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
        SETTINGS.BOOL_Q,...
        SETTINGS.GCD_COEFFICIENT_METHOD...
    );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end


