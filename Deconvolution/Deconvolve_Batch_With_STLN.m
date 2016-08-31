function [arr_hx] = Deconvolve_Batch_With_STLN(arr_fx)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two of the deconvolutions.
%
%
% Input:
%
% set_g - set of input polynomials g(y) to be deconvolved. Each g_{i} has a
% different number of elements, so set_g is a cell array.
%
% Output:
%
% h_{i} = g_{i-1}/g_{i}
%
%

% Global Variables
global SETTINGS

% Get the number of polynomials in the array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% Get number of polynomials in the array arr_hxy of h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx - 1;

% %
% %
% Get the degree m_{i} of each of the polynomials f_{i}

% Initialise vector to store degrees of f_{i}(x)
vDeg_arr_fx = zeros(nPolys_arr_fx,1);

% For each polynomial f_{i}, get its degree.
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% %
% %
% %
% Get the degrees n{i} of polynomials h_{i} = f_{i-1}/f_{i}.
vDeg_arr_hx = (vDeg_arr_fx(1:end-1) - vDeg_arr_fx(2:end));


% Define M to be the total number of coefficeints of all the polynomials
% f_{i} excluding the last f_{i}.
% f_{0},...,f_{nPolys_f-1}.
M = sum(vDeg_arr_fx+1)-(vDeg_arr_fx(end:end)+1);

% Define M1 to be the total number of all coefficients of all of the polynomials
% f_{i}
M1 = sum(vDeg_arr_fx+1);

% Define N to be the number of coefficients of all h_{i}
N = sum(vDeg_arr_hx+1);

%
% y - Preprocess
% n - Dont preprocess
SETTINGS.PREPROC_DECONVOLUTIONS;

switch SETTINGS.PREPROC_DECONVOLUTIONS
    case 'y'
        theta = GetOptimalTheta(arr_fx);
        fprintf([mfilename ' : ' sprintf('Optimal theta : %e \n',theta)])
    case 'n'
        theta = 1;
    otherwise
        error([mfilename ' : error '])
end

% %
% %
% Get preprocessed f(\omega) from f_{i}(x)

% Initialise a cell-array for f(w)
arr_fw = cell(nPolys_arr_fx,1);

% for each f_{i}(x,y) get f_{i}(w,w)
for i = 1:1:nPolys_arr_fx
    arr_fw{i} = GetWithThetas(arr_fx{i},theta);
end

% %
% %
% %
% Write Deconvolutions in form [D^{-1}T(f)Q] h = RHS_f

% %
% %
% Get the right hand side vector of coefficients of f_{\omega}
vRHS_fw = BuildRHSF(arr_fw);

% Get the Left hand side matrix C(f1,...,fd)
DT_fwQ = BuildDTQ(arr_fw);

% Get the solution vector h(w) in the system of equations
% DCQ * hw = RHS_vec.
v_hw = SolveAx_b(DT_fwQ,vRHS_fw);


% %
% %
% %
% Seperate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}
% initialise a cell array to store the coefficients of the individual
% polynomials h_{i}

% Split vec h in to an array of polynomials.

% GetArray of polynomials h_{i}(\omega)
arr_hw = GetArray(v_hw,vDeg_arr_hx);



% %
% %
% Let z be  vectors of perturbations to polynomials f_{i} such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_zw = cell(nPolys_arr_fx,1);

for i =1:1:nPolys_arr_fx
    
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
DYU = BuildDYU(arr_hw,vDeg_arr_fx);

% Compute the initial residual
res_vec = (vRHS_fw + (P*v_zw) - (DT_fwQ*v_hw));

% Set the iteration counter.
ite = 1;

F = eye(N+M1);

G = [DT_fwQ (DYU)-P];



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
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_h = y(1:N);
    delta_z = y(N+1:end);
    
    % Add structured perturbations to vector h.
    v_hw = v_hw + delta_h;
    
    % Add structured perturbations to vector z.
    v_zw = v_zw + delta_z;
    
    % Seperate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    
    arr_zw = GetArray(v_zw,vDeg_arr_fx);
    arr_hw = GetArray(v_hw,vDeg_arr_hx);
    
    % Renew Matrix Pz
    Pz = P*v_zw;
    
    %Increment s in LSE Problem
    s = -(yy-start_point);
    
    %Build iterative DYU
    DYU = BuildDYU(arr_hw,vDeg_arr_fx);
    
    % Build DCEQ
    DC_fQ = BuildDTQ(arr_fw);
    DC_zQ = BuildDTQ(arr_zw);
    
    % Build G
    G = [(DC_fQ + DC_zQ) (DYU-P)];
    
    % Update the RHS_vector
    vRHS_fw = BuildRHSF(arr_fw);
    vRHS_zw = BuildRHSF(arr_zw);
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((vRHS_fw + vRHS_zw ) - ((DC_fQ + DC_zQ)*v_hw));
    
    
    % Get the condition
    condition(ite +1) = norm(res_vec)./norm((vRHS_fw + vRHS_zw));
    
    % Increment iteration number
    ite = ite + 1;
    
end

% Get the array of polynomials h_{i}(x) from h_{i}(w) by removing thetas.
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1:nPolys_arr_hx
    arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
end


% Print outputs to command line
fprintf([mfilename ' : ' 'Performed Deconvolutions...\n'])
fprintf([mfilename ' : ' sprintf('Iterations required for Batch Deconvolution %i\n', ite)]);


fig_name = sprintf([mfilename ' : ' 'Condition at Iterations']);
figure('name',fig_name)
hold on
plot(log10(condition),'-s','DisplayName','Condition')
hold off

end

function Y_new = BuildDYU(arr_hx,vDeg_arr_fx)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

% Inputs.
%
% arr_hw : Set of polynomials h_{i}(w)
%
% vDeg_arr_fx : vector of degrees of polynomials f_{0},...

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
    
    y{i,1} = real(BuildD0Y1U1(hx,deg_fx));
end

%Build the Coefficient Matrix C
num_Rows = 0;
for i = 1:length(vDeg_arr_fx)-1
    num_Rows = num_Rows + 1 + (vDeg_arr_fx(i));
end
cols = (vDeg_arr_fx(1)+1);

xx = zeros(num_Rows,cols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y_new = Y;


end

function Y1 = BuildD0Y1U1(hx,m1)
global SETTINGS

switch SETTINGS.BOOL_LOG
    case 'y' % use logs
        Y1 = BuildD0Y1U1_log(hx,m1);
        
    case 'n' % use nchoosek
        Y1 = BuildD0Y1U1_nchoosek(hx,m1);
    otherwise
        error([mfilename ' : ' 'BOOL_LOG is either (y) or (n)'])
end
end

function Y1 = BuildD0Y1U1_nchoosek(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% to perform the multiplication C(h_{i}) * f_{i} = f_{i+1}
%
% Inputs.
%
% hx : Coefficients of polynomial h_{i}(x)
%
% m1 : Degree of polynomial f_{i}

global SETTINGS

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

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        Y1 = Y1 ./  nchoosek(n1+m1,m1);
    case 'n'
    otherwise
        error('err')
end
end


function Y1 = BuildD0Y1U1_log(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
% m0 = degree of previous polynomial
% m1 = degree of current polynomial

global SETTINGS

% Y1 = zeros(m0+1,m1+1);

% Get Degree of polynomial deg(h_{x}) = n_{1} = m_{0}-m_{1}
n1 = GetDegree(hx);

m0 = n1+m1;


Y1 = [];
% for each column i = 1:1:m0-m1+1
for k = 0:1:m1
    % for each row j = 1:1:m1+1
    for j = k:1:k+n1
        BinomsEval_Log = lnnchoosek(j,k) + lnnchoosek(m0-j,m1-(j-k));
        BinomsEval_Exp = 10.^BinomsEval_Log;
        Y1(j+1,k+1) = hx(j-k+1) .* BinomsEval_Exp;
    end
end
switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Include the denominator
        Denom_Log = lnnchoosek(n1+m1,m1);
        Denom_Exp = 10.^Denom_Log;
        
        Y1 = Y1 ./  Denom_Exp;
    case 'n'
        % Exclude the denominator
    otherwise
        error('err')
end


end




function f = BuildRHSF(fw_array)
% Build the vector f such that it contains the elements of
% Rhs f = [f_{0},...,f_{n-1}]
%
%
% fw = array of vectors f_{0},...,f_{n}
%

% Initialise empty vector.
f = [];

% for each vector f f_{0},...,f_{n-1} in fw_array, add to right hand
% side vector
for i=1:1:length(fw_array)-1
    f = [f;fw_array{i}];
end

end





function DCQ = BuildDTQ(arr_fx)
% set fw is the cell array of poly coefficiencts fw_i
%
% Inputs.

% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% For each of the polynomials f_{i}(x), excluding the final polynomial
for i = 2:1:nPolys_arr_fx
    
    % Get the polynomial f_{i} = set_f{i+1}
    fw = arr_fx{i};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fw_prev = arr_fx{i-1};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fw = GetDegree(fw);
    
    % Get the degree of polynomial f_{i-1}
    deg_fw_prev = GetDegree(fw_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fw_prev - deg_fw;
    
    % Build the Matrix T(f)
    T1 = BuildT1(fw,deg_hw);
    
    D = BuildD(deg_fw,deg_hw);
    Q1 = BuildQ1(deg_hw);
    DT1Q1{i-1}  = D*T1*Q1;
end



%Build the Coefficient Matrix C of all matrices c
DCQ = blkdiag(DT1Q1{1:length(DT1Q1)});

end


function arr_hx = GetArray(v_hw,vDeg_arr_hx)

% Get the number of polynomials in array of h_{i}(x)
nPolys_arr_hx = size(vDeg_arr_hx,1);

% Initialise an array to store the polynomials h_{i}(x)
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_hx
    
    % Get degree of h{i}
    deg_hx = vDeg_arr_hx(i);
    
    % Get coefficients of h_{i} from the solution vector
    arr_hx{i} = v_hw(1:deg_hx+1);
    
    % Remove the coefficients from the solution vector
    v_hw(1:deg_hx+1) = [];
end

end

