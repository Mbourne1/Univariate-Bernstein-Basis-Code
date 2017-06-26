function [] = o_Deconvolution(ex_num, emin, bool_preproc)
% Test the different methods of deconvolving the polynomials f_{i}(x), to
% form the set of polynomials h_{i} where h_{i} = f{i}/f_{i+1}
%
% % Inputs
%
% ex_num : (String) Example number (String)
%
% noise : (Float) noise level
%
% bool_preproc : (Boolean) Bool determining whether to include preprocessing
%
% % Outputs
%
% Results are printed to a file
%
% % Example
%
% >> o_Deconvolution('1', 1e-12, true)


% Set global settings
global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;
SETTINGS.BOOL_LOG = false;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;
SETTINGS.SEED = 1024;


restoredefaultpath();
addpath(genpath('../Examples'));
% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

syms x y;


% Get the array containing symbolic factors and their corresponding
% multiplicity in f_{0}(x)
[arr_factor_mult] = Deconvolution_Examples_Univariate(ex_num);

% Get a vector of the factors of f(x)
vFactors = arr_factor_mult(:,1);

% Get a vector of the multiplicity of the factors of f(x)
vMult = double(arr_factor_mult(:,2));

% Get highest power of any factor
highest_pwr = max(vMult);

% %
% %
% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).



% Initialise the matrix to store multiplicity of each root across its
% columns, for each f_{i}(x).
mult_mat_arr_fx = zeros(highest_pwr+1, length(vMult));

nPolynomials_arr_fx = highest_pwr +1;

arr_sym_fx = cell(nPolynomials_arr_fx, 1);
vDeg_arr_fx = zeros(nPolynomials_arr_fx, 1);

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the factors of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Add the multiplicity structure to the matrix of multiplicities.
    mult_mat_arr_fx(i+1,:) = mults';
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fx{i+1,1} = prod(vFactors.^(mults));
    
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_arr_fx(i+1) = double(feval(symengine, 'degree', (arr_sym_fx{i+1})));
end


% Get the number of polynomials f_{i}(x)
nPolynomials_arr_fx = size(arr_sym_fx, 1);

% Get the number of polynomials h_{i}(x)
nPolynomials_arr_hx = nPolynomials_arr_fx - 1;


% Get multiplicity of the factors of h_{i}(x) as a matrix
% Each of the n rows store multiplicity of each factor for each polynomial
% h_{i}(x).

mult_matrix_arr_hx = abs(diff(mult_mat_arr_fx));


% % Build the set of polynomials h_{i}(x)

% Initialise a cell array to store polynomials h_{i}(x)
arr_symbolic_hx = cell(nPolynomials_arr_hx, 1);

for i = 1:1:nPolynomials_arr_hx
    
    % Get multiplicity structure of h_{i}(x)
    mults = mult_matrix_arr_hx(i, :)';
    
    % Get symbolic polynomial h_{i}(x)
    arr_symbolic_hx{i} = prod(vFactors.^(mults));
    
end



% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDegree_arr_hx = diff(vDeg_arr_fx);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDegree_arr_wx = diff([vDegree_arr_hx; 0]);

% Get the multiplicities of the factors of f(x)
vMultiplicities = find(vDegree_arr_wx~=0);


arr_fx = cell(nPolynomials_arr_fx, 1);
arr_hx_exact = cell(nPolynomials_arr_hx, 1);

for i = 1:1:nPolynomials_arr_fx
    if i <= nPolynomials_arr_hx
        
        arr_fx{i,1} = sym2poly(arr_sym_fx{i})';
        arr_hx_exact{i,1} = sym2poly(arr_symbolic_hx{i})';
        
    else
        
        arr_fx{i,1} = 1;
        
    end
    
end


% %
% %
% Convert the polynomials f_{i}(x) to Bernstein Basis
for i = 1 : 1 : nPolynomials_arr_fx
    
    arr_fx{i,1} = PowerToBernstein(arr_fx{i,1});
    
end

% Convert the polynomials h_{i}(x) to Bernstein Basis
for i = 1 : 1 : nPolynomials_arr_hx
    
    arr_hx_exact{i,1} = PowerToBernstein(arr_hx_exact{i,1});
    
end



% %
% %
% %
% Add noise to the coefficients of f_{i}(x)

% Initialise a cell array to store noisy polynomials f_{i}(x)
arr_fx_noisy = cell(nPolynomials_arr_fx, 1);

% Get noisy polynomials f_{i}(x)
for i = 1 : 1 : nPolynomials_arr_fx
    
    arr_fx_noisy{i,1} = AddNoiseToPoly(arr_fx{i},emin);
    
end


% Define an array of deconvolution methods to be used
arr_DeconvolutionMethod = {...
    'Separate' ...
    'Batch', ...
    'Batch With STLN',...
    'Batch Constrained',...
    'Batch Constrained With STLN'...
    };


nMethods = length(arr_DeconvolutionMethod);

% Testing deconvolution
LineBreakLarge();
arr_hx = cell(nMethods,1);
arr_Error = cell(nMethods,1);

for i = 1 : 1 : nMethods
    
    % Get deconvolution method
    method_name = arr_DeconvolutionMethod{i};
    
    switch method_name
        
        case 'Separate'
            arr_hx{i,1} = Deconvolve_Separate(arr_fx_noisy);
            
        case 'Batch'
            arr_hx{i,1} = Deconvolve_Batch(arr_fx_noisy);
            
        case 'Batch With STLN'
            arr_hx{i,1} = Deconvolve_Batch_With_STLN(arr_fx_noisy);
            
        case 'Batch Constrained'
            arr_hx{i,1} = Deconvolve_Batch_Constrained(arr_fx_noisy, vMultiplicities);
            
        case 'Batch Constrained With STLN'
            arr_hx{i,1} = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy, vMultiplicities);
            
        otherwise
            error('err')
            
    end
    
    fprintf([mfilename ' : ' sprintf('%s \n', method_name )]);
    
    arr_Error{i,1} = GetErrors(arr_hx{i,1}, arr_hx_exact);
    
end








PlotErrors(arr_Error, arr_DeconvolutionMethod)



%--------------------------------------------------------------------------
% Console writing

for i = 1:1:nMethods
    
    methodName = arr_DeconvolutionMethod{i};
    vError = arr_Error{i};
    display([mfilename ' : ' sprintf('Error %s : %e', methodName, mean(vError))]);
    
end



% Initialise array to store error for each method
arr_ErrorNorm = cell(nMethods,1);

for i = 1 : 1 : nMethods
    arr_ErrorNorm{i} = norm(arr_Error{i});
end


PrintToResultsFile(ex_num, bool_preproc, emin, arr_DeconvolutionMethod, arr_ErrorNorm);

end

function [] = PrintToResultsFile(ex_num, bool_preproc, noise, arr_DeconvolutionMethod, arr_ErrorNorm)
%
% % Inputs
%
% ex_num : (String)
%
% bool_preproc : (Boolean)
%
% noise : (float)
%
% arr_DeconvolutionMethod : (Array of Strings)
%
% arr_ErrorNorm : (Array of floats)


fullFileName = sprintf('Results/Results_o_deconvolutions.txt');
nMethods = length(arr_DeconvolutionMethod);

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    
    
    for i = 1:1:nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        
        
        WriteNewLine(method_name, error_norm)
        
        
    end
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    
    WriteHeader()
    for i = 1 : 1 : nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        WriteNewLine(method_name, error_norm);
        
    end
    fclose(fileID);
    
end

    function WriteNewLine(method_name, error_norm)
        
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            num2str(bool_preproc),...
            num2str(noise),...
            method_name,...
            num2str(error_norm)...
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE, EX_NUM, BOOL_PREPROC, NOISE, method_name, error_norm \n');
        
    end


end

function vErrors = GetErrors(arr_hx_comp,arr_hx_exact)
% Compare each computed h{i} with actual h_{i}
%
% % Inputs
%
% arr_hx_comp : (Array of Vectors) Each vector contains coefficients of the
% polynomial h_{i}(x)
%
% arr_hx_exact : (Array of Vectors) Each vector contains coefficients of
% the polynomial h_{i}(x)

% Get number of polynomials in the array
nPolynomials_arr_hx = size(arr_hx_comp,1);

% Initialise vector to store errors
vErrors = zeros(nPolynomials_arr_hx,1);

%
for i = 1:1:nPolynomials_arr_hx
    
    % Get exact polynomial
    exact = arr_hx_exact{i}./ arr_hx_exact{i,1}(1,1);
    
    % Get computed polynomial
    comp = arr_hx_comp{i}./arr_hx_comp{i,1}(1,1);
    
    % Get Error
    vErrors(i) = norm(exact - comp) ./ norm(exact);
    
    
end
end


function PlotErrors(arr_Error, arr_DeconvolutionMethod)

% Plotting


global SETTINGS


nMethods = length(arr_DeconvolutionMethod);

if (SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf([mfilename ':' 'Deconvolution Methods Error']);
    figure('name',figure_name)
    hold on
    
    
    for i = 1 : 1 : nMethods
        
        methodName = arr_DeconvolutionMethod{i};
        
        vError = arr_Error{i,1};
        
        plot(log10(vError), '-o', 'DisplayName', methodName)
        
    end
    
    nPolynomials_hx = length(vError);
    
    legend(gca,'show');
    xlim([1 nPolynomials_hx]);
    xlabel('Factor')
    ylabel('log_{10} error')
    hold off
    
end
end
