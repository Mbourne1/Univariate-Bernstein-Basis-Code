function [arr_hx] = Deconvolve_Batch_Constrained(arr_fx, vMultiplicities)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : (Array of Vectors) Array of polynomials f(x)
%
% vMult : (Vector) Multiplicities of the factors of f(x)
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x)

global SETTINGS


% Get the number of polynomials in the array
nPolys_arr_fx = size(arr_fx, 1);

if( SETTINGS.PREPROC_DECONVOLUTIONS)

    theta = GetOptimalTheta(arr_fx);
    fprintf([mfilename ' : ' sprintf('Optimal theta : %e \n',theta)])
    
else
    
    theta = 1;
    
end

% Initialise a cell-array to store preprocessed polynomials f(\omega)
arr_fw = cell(nPolys_arr_fx, 1);

% Preprocess each polynomial f_{i}(x) to get f_{i}(\omega)
for i = 1:1:nPolys_arr_fx
    
    arr_fw{i} = GetWithThetas(arr_fx{i}, theta);
    
end

% Build LHS Matrix
LHS_Matrix = BuildDTQ(arr_fw, vMultiplicities);

% Build the RHS vector
RHS_vec = BuildRHSF(arr_fw);

% Get number of polynomials in array of h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx - 1;


% Get vector x, the ls solution
x = SolveAx_b(LHS_Matrix, RHS_vec);
x_temp = x;

% Remove non-unique multiplicities.
unique_vMult = unique(vMultiplicities);

% Get the number of polynomials p_{i}(\omega)
nPolys_arr_pw = length(unique_vMult);

% Get unique polynomials p(\omega)
arr_pw = cell(nPolys_arr_pw,1);

for i = 1 : 1 : length(unique_vMult)
    
    % Get multiplicity
    factor_multiplicity = unique_vMult(i);
    
    % Get degree of polynomial p_{i}(\omega)
    factor_degree = GetDegree(arr_fx{factor_multiplicity}) - GetDegree(arr_fx{factor_multiplicity+1});
    
    % Get coefficients of polynomial p_{i}(\omega)
    temp_poly = x_temp(1:factor_degree+1);
    arr_pw{i} = temp_poly;
    
    % Remove coefficients from the vector
    x_temp(1:factor_degree+1) = [];
    
end


% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).

% Get number of entries in array of polynomials p_{i}(x)
nPolys_arr_px = size(arr_pw, 1);

% Initialise the array of h_{i}(x)
arr_hw = cell(nPolys_arr_fx - 1, 1);

% Initialise a count
count = 1;

for i = 1 : 1 : nPolys_arr_px
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1 : 1 : nReps
        
        arr_hw{count,1} = arr_pw{i};
        count = count + 1;
        
    end
    
end

% Get array of polynomials h_{i}(x)
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:nPolys_arr_hx
    
    arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
    
end

end

function LHS_Matrix = BuildDTQ(arr_fx, vMult)
% %
% %
% Build the LHS Coefficient matrix
% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]
%
%
% % Inputs
%
% arr_fx : (Array of Vectors) Array of vectors containing coefficients of
% polynomials f_{i}(x)
%
% vMult : Multiplicity of the factors of f_{0}(x)
% 
%
% % Outputs
%
% LHS_Matrix : (Matrix) Convolution matrix


% Get number of distinct polynomials h_{i}(x)
nDistinctPolys_arr_hx = length(vMult);

% Initialise a cell array to store matrices to form convolution matrix
arr_C = cell(nDistinctPolys_arr_hx);
arr_Q = cell(nDistinctPolys_arr_hx);
arr_DTQ = cell(nDistinctPolys_arr_hx);

for i = 1 : 1 : nDistinctPolys_arr_hx
    
    % Get multiplicity of previous polynonial in the array f_{i-1}(x)
    if i > 1
        old_mult = vMult(i-1);
    else % No previous polynomial so set to zero
        old_mult = 0;
    end
    
    % Get multiplicity of current f_{i}(x)
    new_mult = vMult(i);
    
    
    arr_C{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    
    % Initialise cell arrays.
    arr_Tf = cell((new_mult + 1) - (old_mult+1+1),1);
    arr_D = cell((new_mult + 1) - (old_mult+1+1),1);
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of the previous polynomial f_{i-1}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the current polynomial f_{i}(x)
        fx = arr_fx{j};
        deg_fx = GetDegree(fx);
        
        % Get the degree of the polynomial h_{i} = f_{i-1}/f_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the Cauchy like matrix T_{m_{i} - m_{i-1}}(f_{i})
        arr_Tf{j} = BuildT1(fx, deg_hx);
        
        arr_D{j} = BuildD_2Polys(deg_fx, deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        arr_C{i} = [arr_C{i} ; arr_D{j}*arr_Tf{j}];
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    arr_DTQ{i} = arr_C{i} * arr_Q{i};
    
end

%LHS_Matrix = blkdiag(arr_C{:});
LHS_Matrix = blkdiag(arr_DTQ{:});

end

function [RHS_vec] = BuildRHSF(arr_fx)
% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
%
% % Inputs
%
% arr_fx : (Array of Vectors) Vectors containing coefficients of
% polynomials f(x).
%
% % Outputs
%
% RHS_vec : (Vector)

% Get number of polynomials in the array
nPolys_arr_fx = length(arr_fx);

RHS_vec = [];

for i = 1:1:nPolys_arr_fx - 1
    
    RHS_vec = [RHS_vec ; arr_fx{i}];
    
end

end