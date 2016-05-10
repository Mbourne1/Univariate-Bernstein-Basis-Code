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



SetGlobalVariables(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)

global problemType 
problemType = 'fromRoots'; % fromRoots/fromCoefficients


% Validate Inputs.
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
arr_root_mult_MyMethod = o_roots_mymethod(fx);
fx_mymethod = GetWithoutBinomials(B_poly(arr_root_mult_MyMethod));
error.MyMethod = norm(Normalise(fx_mymethod) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(error.MyMethod);
LineBreakLarge()

% %
% %
% %
% Calculate roots by Musser Method
arr_root_mult_MuzzerMethod = o_roots_Musser(fx);
fx_MuzzerMethod = GetWithoutBinomials(B_poly(arr_root_mult_MuzzerMethod));
error.MuzzerMethod = norm(Normalise(fx_MuzzerMethod) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(error.MuzzerMethod);
LineBreakLarge()

% %
% %
% %
% Calculate roots by matlab 'roots' function.
arr_root_mult_Matlab = o_roots_matlab(fx);
fx_matlab   = GetWithoutBinomials(B_poly(arr_root_mult_Matlab));
error.MatlabMethod   = norm(Normalise(fx_matlab) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(error.MatlabMethod);
LineBreakLarge()


% %
% %
% %
% Calculate roots by 'multroot' function.
arr_root_mult_Multroot = o_roots_multroot(fx);
fx_multroot = GetWithoutBinomials(B_poly(arr_root_mult_Multroot));
error.MultrootMethod = norm(Normalise(fx_multroot) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
display(error.MultrootMethod);
LineBreakLarge()

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


nondistinctRoots_mymethod = GetRepeatedRoots(arr_root_mult_MyMethod);
nondistinctRoots_matlab   = GetRepeatedRoots(arr_root_mult_Matlab);
nondistinctRoots_multroot  = GetRepeatedRoots(arr_root_mult_Multroot);
%nondistinctRoots_bisection  = GetRepeatedRoots(clc_roots_intervalBisection);


%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
        figure('name',figure_name)
        hold on;
    
        scatter((real(nondistinctRoots_mymethod)),imag(nondistinctRoots_mymethod),'yellow','*','DisplayName','My Method');
        scatter((real(nondistinctRoots_matlab)),imag(nondistinctRoots_matlab),'red','DisplayName','Matlab Roots');
        scatter((real(nondistinctRoots_multroot)),imag(nondistinctRoots_multroot),'green','s','filled','DisplayName','MultRoots');
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
PrintToFile(error)


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


function []= PrintToFile(error)

global SETTINGS

fullFileName = 'o_roots_results.txt';

if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    
    str_errors = '%s \t %s \t %s \t %s \t';
    str_globals = '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t';
    str_types = strcat(str_errors,str_globals, '\n');
    
    fprintf(fileID,str_types,...
        error.MyMethod,...
        error.MuzzerMethod,...
        error.MatlabMethod,...
        error.MultrootMethod,...
        SETTINGS.EX_NUM,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA, ...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.APF_METHOD, ...
        SETTINGS.BOOL_Q,...
        SETTINGS.DECONVOLVE_METHOD,...
        SETTINGS.BOOL_LOG);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end


end

