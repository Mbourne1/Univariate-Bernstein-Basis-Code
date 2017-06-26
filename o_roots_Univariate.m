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


PrintToFile(arr_RootFindingMethod, arrForwardErrors, arrBackwardErrors);


end


function [] = PrintToFile(arr_RootFindingMethod, arr_ForwardErrors, arr_BackwardErrors)
% Print results of gcd computation to a text file
%
% % Inputs
%
% arr_RootFindingMethod : (Array of Strings) Array of the names of the root
% finding methods used
%
% arr_errors : (Array of Floats)


global SETTINGS


nMethods = length(arr_RootFindingMethod);

fullFileName = sprintf('Results/Results_o_roots.txt');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1 : 1 : nMethods
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, method_ForwardError);
        end
    end
    
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    for i = 1 : 1 : nMethods
        
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, method_ForwardError)
        end
    end
    fclose(fileID);
    
end



    function WriteNewLine(method_name, myBackwardError, myForwardError)
        % Write a new line of the text file
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
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
            method_name,...
            SETTINGS.DECONVOLUTION_METHOD_HX,...
            SETTINGS.DECONVOLUTION_METHOD_WX,...
            myForwardError,...
            myBackwardError...
            );
    end


    function WriteHeader()
        % If the file doesnt already exist, write a header to the text file
        fprintf(fileID,'DATE, EX_NUM, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD, LOW_RANK_ITE, APF_METHOD, APF_ITE, BOOL_LOG, SYLVESTER_BUILD_METHOD, GCD_METHODM, METHOD_NAME, DECONVOLUTION_METHOD HX, DECONVOLUTION_METHOD WX, FORWARD_ERROR, BACKWARD_ERROR \n');
    end


end





function [my_error] = GetBackwardErrorMeasure(root_mult_arr, fx_exact)
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
my_error  = norm(NormaliseVector(fx_comp) - NormaliseVector(fx_exact)) ./ norm(NormaliseVector(fx_exact));


end

function [my_error] = GetForwardErrorMeasure(root_mult_arr_comp, root_mult_arr_exact)
%
% % Inputs
%
% root_mult_arr_comp : (Matrix)
%
% root_mult_arr_comp : (Matrix)
%
% % Outputs
%
% my_error : (Float)


syms x

% Sort both matrices based on multiplicity

% Sort the exact matrix
[values, order] = sort(root_mult_arr_exact(:,2));
root_mult_arr_exact = root_mult_arr_exact(order,:);


% Sort the computed matrix
[values, order] = sort(root_mult_arr_comp(:,2));
root_mult_arr_comp = root_mult_arr_comp(order,:);


nFactors = size(root_mult_arr_exact, 1);

for i = 1 : 1 : nFactors
    
    % Get the factor
    sym_factor = root_mult_arr_exact(i,1);
    
    % Get coefficients in power basis
    try
        pwr_poly = double((coeffs(sym_factor,x,'All')))';
    catch
        pwr_poly = double((coeffs(sym_factor,x)))';
    end
    
    root = - pwr_poly(2) ./ pwr_poly(1);
    
    root_mult_arr_exact(i,1) = root;
    
end



% Now compare the computed roots
try
    vRoots = double(root_mult_arr_exact(:,1));
    vRoots2 = double(root_mult_arr_comp(:,1));
    
    myDistance = abs(vRoots - vRoots2);
    
    total_error = sum( myDistance );
    
    my_error = total_error;
    
catch
    my_error = 1000;
end
end


function root_mult_matrix = GetRootsAndMultiplicities(fx, str_method)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% str_method : (String) Factorisation method name
%
% % Outputs
%
% root_mult_matrix : (Matrix) where the columns contain the roots and
% corresponding multiplicity of the root.

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
        root_mult_matrix = o_roots_multroot(fx);
        
        
    case 'Bisection Method'
        root_mult_matrix = o_roots_Bisection(fx);
        
    case 'Subdivision Method'
        root_mult_matrix = o_roots_Subdivision(fx);
        
    case 'Bezier Clipping Method'
        root_mult_matrix = o_roots_BezierClipping(fx);
        
        
    otherwise
        error('err')
end

% Print the calculated roots and the corresponding multiplicities.
try
    PrintoutRoots(str_method, root_mult_matrix);
catch
    fprintf('No Roots Found in %s \n', str_method)
end


end