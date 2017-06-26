function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_build_method, rank_revealing_metric, deconvolution_method_hx,...
    deconvolution_method_wx)
% O_ROOTS(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method, sylvester_build_method, ...
%   rank_revealing_metric, deconvolution_method_hx, deconvolution_method_wx)
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
%   * 'T'
%   * 'DT'
%   * 'DTQ'
%   * 'TQ'
%   * 'DTQ Denominator Removed'
%   * 'DTQ Rearranged'
%
% rank_revealing_metric
%   * Singular Values
%   * Max/Min Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals
%
% deconvolution_method_hx
%   * Separate
%   * Batch
%   * Batch With STLN
%   * Batch Constrained
%   * Batch Constrained With STLN
%
%
% deconvolution_method_wx
%   * Separate
%   * Batch
%
%
% % Example
% >> O_ROOTS_UNIVARIATE('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'Minimum Singular Values', 'Batch Constrained', 'Batch')



% Add subfolders
restoredefaultpath
addpath(genpath(pwd));




% Set the global variables
global SETTINGS

SetGlobalVariables_Roots(ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, ...
    apf_method, sylvester_build_method, rank_revealing_metric, ...
    deconvolution_method_hx, deconvolution_method_wx)



% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    
    fprintf('Minimum noise greater than maximum noise \n Swapping values...\n')
    emin_wrong = emin;
    emax_wrong = emax;
    emin = emax_wrong;
    emax = emin_wrong;
    fprintf('')
    
end


% Print Settings to console
PrintGlobalVariables();



% Get the polynomial f(x) as a column vector of coefficients.
[fx_exact, fx_factor_multiplicity_array] = Examples_Roots(ex_num);

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);

% Root Finding Methods
%   My Method
%   Musser Method
%   Yun Method
%   Matlab Method

arr_RootFindingMethod = {...
    %'My Method', ...
    %'Musser Method', ...
    %'Matlab Method', ...
    %'Zeng Method',...
    %'Bisection Method', ...
    %'Subdivision Method', ...
    'Bezier Clipping Method'...
    };
%arr_RootFindingMethod = {'Yun Method'};


% Initialise some arrays
nMethods = length(arr_RootFindingMethod);
arrStartTime = cell(nMethods, 1);
arrRootMultiplicity = cell(nMethods,1);
arrForwardErrors = cell(nMethods,1);
arrBackwardErrors = cell(nMethods,1);
arrDuration = cell(nMethods,1);


% For each method, compute the roots and multiplicities

for i = 1 : 1 : nMethods
    
    LineBreakLarge();
    fprintf(['Method : ' arr_RootFindingMethod{i} '\n']);
    LineBreakLarge();
    LineBreakLarge();
    
    % Get name of root finding Method
    method_name = arr_RootFindingMethod{i};
    
    % Start Timer
    arrStartTime{i} = tic;
    
    % Get roots and multiplicities
    arrRootMultiplicity{i} = GetRootsAndMultiplicities(fx, method_name);
    
    % Get Error Measure
    try
        arrBackwardErrors{i} = GetBackwardErrorMeasure(arrRootMultiplicity{i}, fx_exact);
        arrForwardErrors{i} = GetForwardErrorMeasure(arrRootMultiplicity{i}, fx_factor_multiplicity_array);
        
    catch
        arrForwardErrors{i} = 1000;
    end
    
    % Get Duration
    arrDuration{i} = toc(arrStartTime{i});
    fprintf('Elapsed Time : %2.4f \n', arrDuration{i});
    
    
    LineBreakLarge()
    
end






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
        mat_Root_Mult = arrRootMultiplicity{i};
        
        try
            % Get Vector of real part of each root
            vRealPart = real(mat_Root_Mult(:,1));
            
            % Get Vector of imag part of each root
            vImagPart = imag(mat_Root_Mult(:,1));
            
            % Get scatter data
            scatter( vRealPart, vImagPart, 'DisplayName', methodName);
            
        catch
            fprintf('Roots found by %s method can not be plotted. \n', methodName);
        end
        
    end
    
    
    
    grid on
    xlabel('Real');
    ylabel('Imaginary');
    legend(gca,'show');
    ylim();
    str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
    title(str);
    hold off
    
end


PrintRootsToFile(arr_RootFindingMethod, arrForwardErrors, arrBackwardErrors);


end











end