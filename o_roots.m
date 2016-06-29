function [] = o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
% o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities 
% m_{i}. 
%
% % Inputs
%
% ex_num : (string) Example Number
%
% emin : (float) Noise/Signal maximum threshold (minimum)
%
% emax : (float) Noise/Signal maximum threshold (maximum)
%
% mean_method : (string) method used to compute the mean of entries in C_{n-k}(f)
%               and C_{m-k}(g)
%               'None' - No mean
%               'Geometric Mean Matlab Method'
%               'Geometric Mean My Method'
%              
%
% bool_alpha_theta : (string) {'y,'n''}
%
% low_rank_approx_method : (string) {'None','Standard STLN', 'Standard SNTLN'}
%
% apf_method ('string') {'None', 'Standard APF', 'Root Specific APF'}
%
% % Example
%
% >> o_roots('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard SNTLN','Standard APF')
%
% >> o_roots('Custom:m=10 low=-1 high=1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard SNTLN','Standard APF')


% Set the problem type to a roots type problem.
ProblemType = 'Roots';

% Set the global variables
global SETTINGS
if isempty(SETTINGS)
    fprintf('Set Q and log \n');
    SETTINGS.BOOL_Q = 'y';
    SETTINGS.BOOL_LOG = 'n';
    SETTINGS.ROOTS_HX = 'From Deconvolutions';
    SETTINGS.DECONVOLVE_METHOD = 'Batch';
    SETTINGS.GCD_COEFFICIENT_METHOD = 'ux';
end

SetGlobalVariables(ProblemType, ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)

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



% Add subdirectories
addpath 'Root Finding Methods'

% Get the polynomial f(x) as a column vector of coefficients.
[fx_exact] = Examples_Roots(ex_num);

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = VariableNoise(fx_exact,emin,emax);

% %
% %
% %
% Calculate roots by mymethod.
try
    myMethodStart = tic;
    arr_root_mult_MyMethod = o_roots_mymethod(fx);
    time.MyMethod = toc(myMethodStart);
    errors.MyMethod = GetErrorMeasure(arr_root_mult_MyMethod,fx_exact);
    LineBreakLarge()
    comp_roots.MyMethod = mat2str(arr_root_mult_MyMethod(:,1));
catch err
    fprintf([mfilename ' : ' 'Error computing Roots by My Method \n' ])
    fprintf(err.message);
    errors.MyMethod = 999999999;
    time.MyMethod   = 999999999;
    comp_roots.MyMethod = mat2str([0 0 0 0 0 0]);
end

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
arr_root_mult_MatlabMethod = o_roots_matlab(fx);
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
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
        figure('name',figure_name)
        hold on;
        scatter( real(arr_root_mult_MyMethod(:,1)), imag(arr_root_mult_MyMethod(:,2)),'yellow','*','DisplayName','My Method');
        scatter( real(arr_root_mult_MatlabMethod(:,1)), imag(arr_root_mult_MatlabMethod(:,2)),'red','DisplayName','Matlab Roots');
        scatter( real(arr_root_mult_Multroot(:,1)), imag(arr_root_mult_Multroot(:,2)),'green','s','filled','DisplayName','MultRoots');
        xlabel('Real');
        ylabel('Imaginary');
        legend(gca,'show')
        str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
        title(str);
        hold off
    case 'n'
        % Dont plot graph
    otherwise
        error('error: plot_graphs is either y or n')
end

% %
% %
PrintToFile(comp_roots,errors,time)


end


function [nondistinctRoots_mymthd] = GetRepeatedRoots(mat_Root_Mult) 
% Given the matrix whose columns are a [root,multiplicity] pair, get a
% vector which contains each root r_{i} m_{i} times, where m_{i} is the
% multiplicity of r_{i}.


% Let sum_rt_mult_mymthd be the sum of all of the multiplicities of all of the
% roots obtained by my method
sum_rt_mult_mymthd = sum(mat_Root_Mult(:,2));

% Initialise a vector to store the nondistinct roots
nondistinctRoots_mymthd = zeros(sum_rt_mult_mymthd,1);

% Initialise a count
count = 1;

% for each unique root i
for i = 1:1:size(mat_Root_Mult,1)
    
    % Get multiplicty of root i
    m = mat_Root_Mult(i,2);
    
    % for each of the m roots at r_{i}
    for j = 1:1:m
        
        % Add the root to a vector of nondistinct roots
        nondistinctRoots_mymthd(count,1) = mat_Root_Mult(i,1);
        
        % Increment the counter
        count = count + 1;
    end
end
end


function []= PrintToFile(comp_roots,error,time)

global SETTINGS

fullFileName = 'Results_o_roots.txt';

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
        SETTINGS.BOOL_Q,...
        SETTINGS.DECONVOLVE_METHOD,...
        SETTINGS.BOOL_LOG,...
        SETTINGS.ROOTS_HX,...
        comp_roots.MyMethod...
    );
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end


end


function [error] = GetErrorMeasure(arr_root_mult,fx_exact)
% Get the distance between the polynomial f_{comp} and f_{exact}
%
% Inputs.
%
% arr_root_mult : Matrix whose rows contain a computed root of f(x) and 
%                 its corresponding multiplicitiy.
%
% fx_exact : Coefficients of exact Polynomial f(x)
%

% Get coefficients of computed polynomial f(x)
fx_comp = GetWithoutBinomials(B_poly(arr_root_mult));

% Get distance between f_{comp}(x) and f_{exact}(x)
error  = norm(Normalise(fx_comp) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(error);

end