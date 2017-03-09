function [arr_hx] = Deconvolve_Batch_With_STLN(arr_fx)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two of the deconvolutions.
%
%
% Input:
%
% arr_fx - set of input polynomials g(y) to be deconvolved. Each g_{i} has a
% different number of elements, so set_g is a cell array.
%
% Output:
%
% arr_hx (Array of Vectors) Vectors containing coefficients of polynomials
% h_{i}(x) where h_{i} = f_{i}/f_{i+1}
%


% Global Variables
global SETTINGS

% Get the number of polynomials in the array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% Get number of polynomials in the array arr_hxy of h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx - 1;


% Get the degree m_{i} of each of the polynomials f_{i}

% Initialise vector to store degrees of f_{i}(x)
vDeg_arr_fx = zeros(nPolys_arr_fx, 1);

% For each polynomial f_{i}, get its degree.
for i = 1 : 1 : nPolys_arr_fx
    
    vDeg_arr_fx(i) = GetDegree( arr_fx{i});
    
end

% Get the degrees n{i} of polynomials h_{i} = f_{i-1}/f_{i}.
vDeg_arr_hx = (vDeg_arr_fx(1:end-1) - vDeg_arr_fx(2:end));


% Define M to be the total number of coefficeints of all the polynomials
% f_{i} excluding the last f_{i}.
% f_{0},...,f_{nPolys_f-1}.
M = sum(vDeg_arr_fx+1) - (vDeg_arr_fx(end:end)+1);

% Define M1 to be the total number of all coefficients of all of the polynomials
% f_{i}
M1 = sum(vDeg_arr_fx+1);

% Define N to be the number of coefficients of all h_{i}
N = sum(vDeg_arr_hx+1);

% Preprocess the polynomials to get optimal value \theta
if(SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx);
    fprintf([mfilename ' : ' sprintf('Optimal theta : %e \n',theta)])
    
else
    
    theta = 1;
    
end

% % Get preprocessed f(\omega) from f_{i}(x)
% Initialise a cell array for preprocessed polynomials f(w)
arr_fw = cell(nPolys_arr_fx, 1);

% Preprocess each f_{i}(x) get f_{i}(w)
for i = 1:1:nPolys_arr_fx
    
    arr_fw{i} = GetWithThetas(arr_fx{i}, theta);
    
end

% % Write Deconvolutions in form [D^{-1}T(f)Q] h = RHS_f

% Get the right hand side vector of coefficients of f_{\omega}
vRHS_fw = BuildRHSF(arr_fw);

% Get the Left hand side matrix C(f1,...,fd)
DTQ = BuildDTQ(arr_fw);

% Get the solution vector h(w) in the system of equations
% DCQ * hw = RHS_vec.
v_hw = SolveAx_b(DTQ, vRHS_fw);


% Seperate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}

% GetArray of polynomials h_{i}(\omega)
arr_hw = GetArray(v_hw,vDeg_arr_hx);

% Let z be  vectors of perturbations to polynomials f_{i} such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_zw = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    % initialise polynomial z_{i} as a zero vector.
    arr_zw{i,1} = zeros(vDeg_arr_fx(i)+1,1);
    
    
end

% Build vector z, consisting of all vectors z_{i}
v_zw = cell2mat(arr_zw);

% Build the Matrix P
P = [eye(M) zeros(M,M1-M)];

% Get Vector of perturbations for RHS by multiplying perturbation vector by
% P, such that we eliminate the z_max

% Build Matrix Y, where E(z)h = Y(h)z
DYU = BuildDYU(arr_hw, vDeg_arr_fx);

% Compute the initial residual
res_vec = (vRHS_fw + (P*v_zw) - (DTQ*v_hw));

% Set the iteration counter.
ite = 1;

F = eye(N+M1);

G = [DTQ (DYU)-P];



condition(ite) = norm(res_vec)./ norm(vRHS_fw);

start_point = ...
    [
    v_hw;
    v_zw;
    ];

yy = start_point;

s = -(yy - start_point);
% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F, s, G, res_vec);
    
    % Update vector yy
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_h = y(1:N);
    delta_z = y(N+1:end);
    
    % Add structured perturbations to vector h(\omega).
    v_hw = v_hw + delta_h;
    
    % Add structured perturbations to vector z(\omega).
    v_zw = v_zw + delta_z;
    
    % Seperate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    arr_zw = GetArray(v_zw, vDeg_arr_fx);
    arr_hw = GetArray(v_hw, vDeg_arr_hx);
    
    
    %Increment s in LSE Problem
    s = -(yy - start_point);
    
    %Build iterative DYU
    DYU = BuildDYU(arr_hw, vDeg_arr_fx);
    
    % Build DCEQ
    DC_fQ = BuildDTQ(arr_fw);
    DC_zQ = BuildDTQ(arr_zw);
    
    % Build G
    G = [(DC_fQ + DC_zQ) (DYU - P)];
    
    % Update the RHS_vector
    vRHS_fw = BuildRHSF(arr_fw);
    vRHS_zw = BuildRHSF(arr_zw);
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((vRHS_fw + vRHS_zw ) - ((DC_fQ + DC_zQ)*v_hw));
    
    
    % Get the condition
    condition(ite +1) = norm(res_vec)./norm((vRHS_fw + vRHS_zw));
    
    % Increment iteration number
    ite = ite + 1;
    
end % End of loop

% Get the array of polynomials h_{i}(x) from h_{i}(w) by removing thetas.
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1:nPolys_arr_hx
    
    arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
    
end


% Print outputs to command line
fprintf([mfilename ' : ' 'Performed Deconvolutions...\n'])
fprintf([mfilename ' : ' sprintf('Iterations required for Batch Deconvolution %i\n', ite)]);

if (SETTINGS.PLOT_GRAPHS)
    fig_name = sprintf([mfilename ' : ' 'Condition at Iterations']);
    figure('name',fig_name)
    hold on
    plot(log10(condition),'-s','DisplayName','Condition')
    hold off
end

end

function Y_new = BuildDYU(arr_hx, vDeg_arr_fx)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

% Inputs.
%
% arr_hw : (Array of Vectors) Set of polynomials h_{i}(x)
%
% vDeg_arr_fx : (Vector) Vector of degrees of polynomials f_{0},...

nPolys_arr_hx = size(arr_hx,1);

for i = 1:1:nPolys_arr_hx
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hx = arr_hx{i,1};
    
    % Get degree of f_{i}
    deg_fx = vDeg_arr_fx(i+1);
    
    y{i,1} = real(BuildD0Y1U1(hx, deg_fx));
end

%Build the Coefficient Matrix C
nRows = 0;
for i = 1:length(vDeg_arr_fx)-1
    nRows = nRows + 1 + (vDeg_arr_fx(i));
end
nCols = (vDeg_arr_fx(1)+1);

xx = zeros(nRows,nCols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y_new = Y;


end

function Y1 = BuildD0Y1U1(hx, m1)
global SETTINGS

if( SETTINGS.BOOL_LOG)
    
    Y1 = BuildD0Y1U1_log(hx,m1);
    
else
    Y1 = BuildD0Y1U1_nchoosek(hx,m1);
    
end
end

function Y1 = BuildD0Y1U1_nchoosek(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% to perform the multiplication C(h_{i}) * f_{i} = f_{i+1}
%
% Inputs.
%
% hx : (Vector) Coefficients of polynomial h_{i}(x)
%
% m1 : (Int) Degree of polynomial f_{i}

% Get degree of polynomial h(x) where deg(h_{1}) = n_{1} = m_{0} - m_{1}
n1 = GetDegree(hx);

% Y1 = zeros(m0+1,m1+1);
Y1 = [];

% for each column j = 1:1:m0-m1+1
for j = 0:1:m1
    % for each row i = 1:1:m1+1
    for i = j:1:j+n1
        Y1(i+1,j+1) = ...
            hx(i-j+1) .* nchoosek(i,j) .* nchoosek(n1+m1-i,m1-j);
        
    end
end


Y1 = Y1 ./  nchoosek(n1+m1,m1);

end


function Y1 = BuildD0Y1U1_log(hx, m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
%
% % Inputs
%
% hx : (Vector) Vector of Coefficients of h(x)
%
% m1 : (Int) degree of current polynomial
%
% % Outputs
%
% Y1 : (Matrix)

% Get Degree of polynomial deg(h_{x}) = n_{1} = m_{0}-m_{1}
n1 = GetDegree(hx);

m0 = n1+m1;


% for each column i = 1:1:m0-m1+1
for k = 0:1:m1
    
    % for each row j = 1:1:m1+1
    for j = k:1:k+n1
        BinomsEval_Log = lnnchoosek(j,k) + lnnchoosek(m0-j,m1-(j-k));
        BinomsEval_Exp = 10.^BinomsEval_Log;
        Y1(j+1,k+1) = hx(j-k+1) .* BinomsEval_Exp;
    end
end

% Include the denominator
Denom_Log = lnnchoosek(n1+m1,m1);
Denom_Exp = 10.^Denom_Log;

Y1 = Y1 ./  Denom_Exp;

end





function f = BuildRHSF(fx_array)
% Build the vector f such that it contains the coefficients of each of the
% polynomials f_{i}(x) for i = 0,...,n-1
% Rhs f = [f_{0},...,f_{n-1}]
%
% % Inputs
%
% fx = (Array of vectors) array of vectors f_{0},...,f_{n}
%
%
% % Outputs
%
% f : (Vector) Coefficients of polynomials f_{0}...,f_{n-1}

% Initialise empty vector.
f = [];

% Get number of polynomials in the array
nPolys_arr_fx = length(fx_array);

% For each vector f f_{0},...,f_{n-1} in fw_array, add to right hand
% side vector
for i = 1 : 1 : nPolys_arr_fx - 1
    
    f = [f ; fx_array{i}];
    
end

end





function DCQ = BuildDTQ(arr_fx)
% Build the left hand side convolution matrix for the set of deconvolutions
%
%
% Inputs
%
% arr_fx : (Array of Vectors) Array of vectors containing coefficients of
% the polynomials f_{i}(x)

% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% For each of the polynomials f_{i}(x), excluding the final polynomial
arr_DT1Q1 = cell(nPolys_arr_fx -1);

for i = 2:1:nPolys_arr_fx
    
    % Get the polynomial f_{i} = set_f{i+1}
    fx = arr_fx{i};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fx_prev = arr_fx{i-1};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fx = GetDegree(fx);
    
    % Get the degree of polynomial f_{i-1}
    deg_fx_prev = GetDegree(fx_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fx_prev - deg_fx;
    
    % Build the Matrix T(f)
    T1 = BuildT1(fx, deg_hw);
    
    % Build the matrix D^{-1}_{}
    D = BuildD_2Polys(deg_fx, deg_hw);
    
    % Build the matrix Q_{}
    Q1 = BuildQ1(deg_hw);
    
    % Add the matrix DTQ to an array of matrices ready to construct the LHS
    % matrix
    arr_DT1Q1{i-1}  = D*T1*Q1;
    
end

%Build the Coefficient Matrix C of all matrices c
DCQ = blkdiag(arr_DT1Q1{1:length(arr_DT1Q1)});

end


function arr_hx = GetArray(v_hx, vDeg_arr_hx)
% Separate the coefficients of each polynomial h_{i}(x) into an array of
% vectors.
%
% % Inputs
%
% v_hx : (Vector) Vector containing coefficients of all polynomials
% h_{i}(x)
%
% % Outputs
%
% arr_hx : (Array of Vectors) Array containing vectors of coefficients of
% polynomials h_{i}(x)

% Get the number of polynomials in array of h_{i}(x)
nPolys_arr_hx = size(vDeg_arr_hx,1);

% Initialise an array to store the polynomials h_{i}(x)
arr_hx = cell(nPolys_arr_hx, 1);

for i = 1 : 1 : nPolys_arr_hx
    
    % Get degree of h{i}
    deg_hx = vDeg_arr_hx(i);
    
    % Get coefficients of h_{i} from the solution vector
    arr_hx{i} = v_hx(1:deg_hx+1);
    
    % Remove the coefficients from the large vector
    v_hx(1:deg_hx+1) = [];
    
end

end

