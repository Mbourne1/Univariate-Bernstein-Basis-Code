function [] = Test_Deconvolution(ex_num,noise,bool_preproc)
% Test the different methods of deconvolving the polynomials f_{i}(x), to
% form the set of polynomials h_{i} where h_{i} = f{i}/f_{i+1}
%
% % Inputs
% 
% ex_num : Example number (String)
%
% noise : noise level
%
% bool_preproc : Bool determining whether to include preprocessing
%
%
% % Example
%
% >> Test_Deconvolution('1',1e-12,'y')

% Set settings pertaining to this test

global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;
SETTINGS.BOOL_LOG = 'n';
SETTINGS.BOOL_DENOM_SYL = 'y';
SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;
SETTINGS.SEED = 1024;



[vFactors,vMult] = Examples_Deconvolution(ex_num);

% Get highest power of any factor
highest_pwr = max(vMult);

% %
% %
% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).

% Initialise the matrix to store multiplicity of each root across its
% columns, for each f_{i}(x).
mult_mat_arr_fx = zeros(highest_pwr+1,length(vMult));

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Add the multiplicity structure to the matrix of multiplicities.
    mult_mat_arr_fx(i+1,:) = mults';
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fx{i+1,1} = prod(vFactors.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_arr_fx(i+1) = double(feval(symengine, 'degree', (arr_sym_fx{i+1})));
end

display(arr_sym_fx{1});

% Get the number of polynomials f_{i}(x)
nPolys_arr_fx = size(arr_sym_fx,1);

% Get the number of polynomials h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx - 1;


% %
% %
% Get coefficients vectors of h_{i}(x)

% Get multiplicity of the factors of h_{i}(x)
mult_mat_arr_hx = abs(diff(mult_mat_arr_fx));

% Initialise a cell array to store polynomials h_{i}(x)
arr_sym_h = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_hx

    % Get multiplicity structure of h_{i}(x)
    mults = mult_mat_arr_hx(i,:)';
   
    % Get symbolic polynomial h_{i}(x)
    arr_sym_h{i} = prod(vFactors.^(mults));
   
end



% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDeg_arr_hx = diff(vDeg_arr_fx);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDeg_arr_wx = diff([vDeg_arr_hx 0]);

% Get the multiplicities of the roots.
vMultiplicities = find(vDeg_arr_wx~=0);


arr_fx = cell(nPolys_arr_fx,1);
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_fx
    if i <= nPolys_arr_hx
        arr_fx{i,1} = sym2poly(arr_sym_fx{i})';
        arr_hx{i,1} = sym2poly(arr_sym_h{i})';
    else
        arr_fx{i,1} = 1;
    end
    
end


% %
% %
% Convert the polynomials f_{i}(x) to Bernstein Basis
for i = 1:1:length(arr_fx)
    arr_fx{i,1} = PowerToBernstein(arr_fx{i,1});
end

% Convert the polynomials h_{i}(x) to Bernstein Basis
for i = 1:1:length(arr_hx)
    arr_hx{i,1} = PowerToBernstein(arr_hx{i,1});
end



% %
% %
% %
% Add noise to the coefficients of f_{i}(x)

% Initialise a cell array to store noisy polynomials f_{i}(x)
arr_fx_noisy = cell(nPolys_arr_fx,1);

% Get noisy polynomials f_{i}(x)
for i = 1:1:nPolys_arr_fx
    arr_fx_noisy{i,1} = AddNoiseToPoly(arr_fx{i},noise);
end



%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolution Separate \n']);


arr_hx_Separate = Deconvolve_Separate(arr_fx_noisy);
vError_Separate = GetErrors(arr_hx_Separate,arr_hx);


%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolution Batch\n']);

arr_hx_Batch = Deconvolve_Batch(arr_fx_noisy);
vError_Batch = GetErrors(arr_hx_Batch,arr_hx);

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolution Batch With STLN\n']);

arr_hx_BatchSTLN = Deconvolve_Batch_With_STLN(arr_fx_noisy);
vError_BatchSTLN = GetErrors(arr_hx_BatchSTLN,arr_hx);

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf([mfilename ' : ''Deconvoltuion Batch Constrained \n']);

arr_hx_BatchConstrained = Deconvolve_Batch_Constrained(arr_fx_noisy,vMultiplicities);
vError_BatchConstrained = GetErrors(arr_hx_BatchConstrained,arr_hx);

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf([mfilename ' : ''Deconvoltuion Batch Constrained With STLN \n']);

arr_hx_BatchConstrainedSTLN = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy,vMultiplicities);
vError_BatchConstrainedSTLN = GetErrors(arr_hx_BatchConstrainedSTLN,arr_hx);
% -------------------------------------------------------------------------

% Plotting
nPolys_hx = size(arr_hx,1);
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure()
        hold on
        plot(log10(vError_Separate),'-o','DisplayName','Separate')
        plot(log10(vError_Batch),'-s','DisplayName','Batch')
        plot(log10(vError_BatchSTLN),'-s','DisplayName','Batch STLN')
        plot(log10(vError_BatchConstrained),'-s','DisplayName','Batch Constrained')
        plot(log10(vError_BatchConstrainedSTLN),'-s','DisplayName','Batch Constrained STLN')
        legend(gca,'show');
        xlim([1 nPolys_hx]);
        xlabel('Factor')
        ylabel('log_{10} error')
        hold off
    case 'n'
    otherwise
        error('err');
end

%--------------------------------------------------------------------------
% Console writing



display([mfilename ' : ' sprintf('Error Separate : %e',norm(vError_Separate))]);
display([mfilename ' : ' sprintf('Error Batch : %e', norm(vError_Batch))]);
display([mfilename ' : ' sprintf('Error Batch STLN: %e', norm(vError_BatchSTLN))]);
display([mfilename ' : ' sprintf('Error Batch Constrained : %e', norm(vError_BatchConstrained))]);
display([mfilename ' : ' sprintf('Error Batch Constrained STLN : %e', norm(vError_BatchConstrainedSTLN))]);

%--------------------------------------------------------------------------
% Writing outputs to file
A = ...
    [
    norm(vError_Separate),...
    norm(vError_Batch),...
    norm(vError_BatchSTLN),...
    norm(vError_BatchConstrained),...
    norm(vError_BatchConstrainedSTLN)
    ];

fileID = fopen('Deconvolution/Test_Deconvolution.txt','a');
fprintf(fileID,'%s %s %6.2e %6.2e %6.2e %6.2e %6.2e %6.2e\n',ex_num,bool_preproc,noise,A);
fclose(fileID);

end


function vErrors = GetErrors(arr_hx_comp,arr_hx_exact)
% Compare each computed h{i} with actual h_{i}

nPolys_hx = size(arr_hx_comp,1);

vErrors = zeros(nPolys_hx,1);

for i = 1:1:nPolys_hx
    
    exact = arr_hx_exact{i}./ arr_hx_exact{i,1}(1,1);
    comp = arr_hx_comp{i}./arr_hx_comp{i,1}(1,1);
    
    vErrors(i) = norm(exact - comp) ./ norm(exact);
    
    
end
end
