function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method)
% O_ROOTS(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities
% m_{i}.
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Noise/Signal maximum threshold (minimum)
%
% emax : (Float) Noise/Signal maximum threshold (maximum)
%
% mean_method : (string) method used to compute the mean of entries in C_{n-k}(f)
%               and C_{m-k}(g)
%           'None' - No mean
%           'Geometric Mean Matlab Method'
%           'Geometric Mean My Method'
%
%
% bool_alpha_theta : (Boolean)
%           true
%           false
%
% low_rank_approx_method : (string)
%           'None'
%           'Standard STLN'
%           'Standard SNTLN'
%           'Root Specific SNTLN'
%
% apf_method ('string')
%           'None'
%           'Standard APF NonLinear'
%           'Standard APF Linear'
%
% Sylvester_Build_Method
%           'T'
%           'DT'
%           'DTQ'
%           'TQ'
%           'DTQ Rearranged Denom Removed'
%           'DTQ Rearranged'
%
% % Example
% >> O_ROOTS_UNIVARIATE('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ')
% >> O_ROOTS_UNIVARIATE('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard SNTLN', 'Standard APF NonLinear','DTQ')
%
% >> O_ROOTS_UNIVARIATE('Custom:m=10 low=-1 high=1',1e-12,1e-10,'Geometric Mean Matlab Method',true,'Standard SNTLN','Standard APF','DTQ')


% add paths
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
addpath(genpath('Deconvolution'));
addpath(genpath('Examples'));
addpath(genpath('GCD Methods'));
addpath(genpath('Get Cofactor Coefficients'));
addpath(genpath('Get GCD Coefficients'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Low Rank Approx'));
addpath(genpath('Root Finding Methods'));

% Set the problem type to a roots type problem.
ProblemType = 'Roots';

% Set the global variables
global SETTINGS
if isempty(SETTINGS)
    
    fprintf('Set Q and log \n');
    SETTINGS.BOOL_LOG = false;
    SETTINGS.ROOTS_HX = 'From Deconvolutions';
    SETTINGS.GCD_COEFFICIENT_METHOD = 'ux';
    
end

SetGlobalVariables(ProblemType, ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method,sylvester_build_method)

% Use the examples from a set of given roots
global Example_Type
Example_Type = 'From Roots'; % fromRoots/fromCoefficients

% Print Settings to console
PrintGlobalVariables();

% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    fprintf('minimum noise greater than maximum noise \n swapping values...\n')
    emin_wrong = emin;
    emax_wrong = emax;
    emin = emax_wrong;
    emax = emin_wrong;
end




% Get the polynomial f(x) as a column vector of coefficients.
[fx_exact] = Examples_Roots(ex_num);

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = AddVariableNoiseToPoly(fx_exact,emin,emax);

% %
% %
% %
% Calculate roots by mymethod.
%try
% Start timer
myMethodStart = tic;

% Get roots by my method
arr_root_mult_MyMethod = o_roots_mymethod(fx);

% End timer
time.MyMethod = toc(myMethodStart);

% Get error
errors.MyMethod = GetErrorMeasure(arr_root_mult_MyMethod,fx_exact);
LineBreakLarge()
comp_roots.MyMethod = mat2str(arr_root_mult_MyMethod(:,1));

%catch err
%   fprintf([mfilename ' : ' 'Error computing Roots by My Method \n' ])
%
%   fprintf(err.message);
%  errors.MyMethod = 999999999;
% time.MyMethod   = 999999999;
%  comp_roots.MyMethod = mat2str([0 0 0 0 0 0]);
%end



% %
% %
% %
% Calculate roots by Musser Method
try
    MusserMethodStart = tic;
    arr_root_mult_MusserMethod = o_roots_Musser(fx);
    time.MusserMethod = toc(MusserMethodStart);
    errors.MusserMethod = GetErrorMeasure(arr_root_mult_MusserMethod,fx_exact);
    LineBreakLarge()
    
catch err
    
    fprintf([mfilename ' : ' 'Error computing Roots by Musser Method \n' ])
    fprintf(err.message);
    errors.MusserMethod = 9999999;
    time.MusserMethod = 9999999;
end


% %
% %
% %
% Calculate roots by matlab 'roots' function.
arr_root_mult_MatlabMethod = o_roots_Matlab(fx);
errors.MatlabMethod = GetErrorMeasure(arr_root_mult_MatlabMethod,fx_exact);
LineBreakLarge()


% %
% %
% %
% Calculate roots by 'multroot' function.
try
    arr_root_mult_Multroot = o_roots_multroot(fx);
    errors.MultrootMethod = GetErrorMeasure(arr_root_mult_Multroot,fx_exact);
    LineBreakLarge()
catch
    errors.MultrootMethod = 9999999;
end


% Calculate roots by 'Interval Bisection' function
% clc_roots_intervalBisection = o_roots_bisection(fx);

% Calculate roots by 'Subdivisiton' Method
%clc_roots_subdivision = o_roots_subdivision(fx);

% Calculate roots by bezier clipping
%clc_roots_clipping = o_roots_BezierClipping(fx);




% %


% Get vector of roots calculated by my method, and extract multiplicities.
% eg: if r_{1} has multiplicity 5, then the vector would be given by
% X1 = [r1 r1 r1 r1 r1 ...].


%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if(SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
    figure('name',figure_name)
    hold on;
    scatter( real(arr_root_mult_MyMethod(:,1)), imag(arr_root_mult_MyMethod(:,1)),'yellow','*','DisplayName','My Method');
    scatter( real(arr_root_mult_MatlabMethod(:,1)), imag(arr_root_mult_MatlabMethod(:,1)),'red','DisplayName','Matlab Roots');
%    scatter( real(arr_root_mult_Multroot(:,1)), imag(arr_root_mult_Multroot(:,1)),'green','s','filled','DisplayName','MultRoots');
    grid on
    
    xlabel('Real');
    ylabel('Imaginary');
    legend(gca,'show')
    ylim()
    str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
    title(str);
    hold off
    
end

% %
% %
PrintToFile(comp_roots,errors,time)


end



function []= PrintToFile(comp_roots,error,time)

global SETTINGS

fullFileName = 'Results/Results_o_roots.txt';

if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    
    str_errors = ' %s, \t %s, \t %s, \t %s, \t %s, \t';
    str_globals = '%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \t';
    str_types = strcat(str_errors,str_globals, '\n');
    
    fprintf(fileID,str_types,...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        num2str(error.MyMethod),...
        num2str(time.MyMethod),...
        num2str(error.MusserMethod),...
        num2str(time.MusserMethod),...
        num2str(error.MatlabMethod),...
        num2str(error.MultrootMethod),...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA, ...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.APF_METHOD, ...
        SETTINGS.SYLVESTER_BUILD_METHOD,...
        SETTINGS.DECONVOLVE_METHOD_HX_FX,...
        SETTINGS.BOOL_LOG,...
        SETTINGS.DECONVOLVE_METHOD_WX_HX,...
        comp_roots.MyMethod...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end


end


function [my_error] = GetErrorMeasure(root_mult_arr, fx_exact)
% Get the distance between the polynomial f_{comp} and f_{exact}
%
% Inputs.
%
% root_mult_arr : Matrix whose rows contain a computed root of f(x) and
%                 its corresponding multiplicitiy.
%
% fx_exact : Coefficients of exact Polynomial f(x)
%
% % Outputs
%
% Error : Measure of error between the exact polynomial f(x) and the
% polynomial f(x) computed by the root finding method.

% Get coefficients of computed polynomial f(x)
fx_comp = GetWithoutBinomials( BuildPolyFromRoots( root_mult_arr));

% Get distance between f_{comp}(x) and f_{exact}(x)
my_error  = norm(Normalise(fx_comp) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(my_error);

end