function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method)
% O_ROOTS(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method)
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


% Add subfolders
restoredefaultpath
addpath(genpath(pwd));

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
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);



% Root Finding Methods
%   My Method
%   Musser Method
%   Yun Method
%   Matlab Method

arr_RootFindingMethod = {'My Method', 'Musser Method', 'Matlab Method'};
%arr_RootFindingMethod = {'Yun Method'};


% Initialise some arrays
nMethods = length(arr_RootFindingMethod);
arr_StartTime = cell(nMethods, 1);
arr_root_mult = cell(nMethods,1);
arr_errors = cell(nMethods,1);
arr_Duration = cell(nMethods,1);


% For each method, compute the roots and multiplicities

for i = 1 : 1 : nMethods
    
    LineBreakLarge();
    fprintf('');
    fprintf(['Method : ' arr_RootFindingMethod{i} '\n']);
    fprintf('');
    LineBreakLarge();
    LineBreakLarge();
    
    % Start Timer
    arr_StartTime{i} = tic;
    
    % Get roots and multiplicities
    arr_root_mult{i} = GetRootsAndMultiplicities(fx, arr_RootFindingMethod{i});
    
    % Get Error Measure
    arr_errors{i} = GetErrorMeasure(arr_root_mult{i}, fx_exact);    
    
    % Get Duration 
    arr_Duration{i} = toc(arr_StartTime{i});
    fprintf('Elapsed Time : %2.4f \n',arr_Duration{i});
    
    LineBreakLarge()
    
end



% Calculate roots by 'Interval Bisection' function
clc_roots_intervalBisection = o_roots_Bisection(fx);

% Calculate roots by 'Subdivisiton' Method
clc_roots_subdivision = o_roots_Subdivision(fx);

% Calculate roots by bezier clipping
clc_roots_clipping = o_roots_BezierClipping(fx);




%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if (SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
    figure('name',figure_name)
    hold on;
    
    
    for i = 1:1:length(arr_RootFindingMethod)
        
        % Get method name
        methodName = arr_RootFindingMethod{i};
        
        % Get matrix of roots and corresponding multiplicities
        mat_Root_Mult = arr_root_mult{i};
        
        % Get Vector of real part of each root
        vRealPart = real(mat_Root_Mult(:,1));
        
        % Get Vector of imag part of each root
        vImagPart = imag(mat_Root_Mult(:,1));
        
        % Get scatter data
        scatter( vRealPart, vImagPart, 'red','DisplayName', methodName);
        
    end
    
    
    
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


end


function root_mult_matrix = GetRootsAndMultiplicities(fx, str_method)
% 
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% str_method : (String)

switch str_method
    case 'My Method'
        root_mult_matrix = o_roots_mymethod(fx);
        
    case 'Musser Method'
        root_mult_matrix = o_roots_Musser(fx);
        
    case 'Yun Method'
        root_mult_matrix = o_roots_Yun(fx);
        
    case 'Matlab Method'
        root_mult_matrix = o_roots_Matlab(fx);
        
    case 'Zeng Method'
      
        
    otherwise
        error('err')
end

end