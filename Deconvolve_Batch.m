function arr_hw = Deconvolve_Batch(set_f)
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

global MAX_ITERATIONS_DECONVOLUTIONS
global MAX_ERROR_DECONVOLUTIONS


% let d be the number of deconvolutions = num of polynomials in set_g - 1
d = length(set_f) -1;

% Get values of m_{i}, the degrees of the polynomials f_{i}(y) for i =
% 0...d
m = zeros(1,length(set_f));
for i = 1:1:length(set_f)
    m(i) = GetDegree(set_f{i});
end

% get values of n{i}, the degree of polynomials h_{i} for i = 1....d
n = zeros(1,d);
for i = 1:1:d
    n(i) = m(i)-m(i+1);
end

%Define M to be the total number of all coefficients of the first d polynomials
%g_{0}...g_{d-1}, this is used in the RHS vector.
M = sum(m+1)-(m(end:end)+1);

% Define M1 to be the total number of all coefficients of polynomials
% g_{0}...g_{d}
M1 = sum(m+1) ;% m_{i}+1

% Define N to be the degree of sum of all h_{i}
N = sum(n+1);

% Obtain theta such that the ratio of max element to min element is
% minimised
theta = getOptimalTheta(set_f,m);


% Initialise a cell-array for f(w)
fw = cell(1,length(set_f));


for i = 1:1:length(set_f)
    fx = set_f{i};
    % Convert to scaled bernstein basis, by inclusion of theta
    % a_{i}\theta^{i}
    %fw{i} = fx*theta^(i);
    fw{i} = GetWithThetas(fx,theta);
    
end

% Write Deconvolutions in form [D^{-1}C(f)Q] h = RHS_f

RHS_f = real(BuildRHSF(fw));
DCQ = BuildDCQ(fw);

% Solve h_{0}
hw = SolveAx_b(DCQ,RHS_f);
v_h = hw;

% Seperate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}
% initialise a cell array to store the coefficients of the individual
% polynomials h_{i}

arr_hw = cell(1,length(n));
for i = 1:1:length(n)
    arr_hw{i} = hw(1:n(i)+1);
    hw(1:n(i)+1) = [];
end

hw = v_h;

% Let z be  vectors of perturbations to polynomials fi such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_z = cell(1,length(m));
for i =1:1:length(m)
    arr_z{i} = zeros(1,m(i)+1);
end

% Build vector z, consisting of all vectors z_{i}
z_o = [arr_z{1:length(arr_z)}];
z_o = z_o'    ;

% Build the Matrix P
P = [eye(M) zeros(M,M1-M)];

% Get Vector of perturbations for RHS by multiplying perturbation vector by
% P, such that we eliminate the z_max

% Build Matrix Y, where E(z)h = Y(h)z
DYU = BuildDYU(arr_hw,m);

% Compute the first residual
residual_o = (RHS_f+(P*z_o) - (DCQ*v_h));

% Set the iteration counter.
ite = 0;

F = eye(N+M1);

G = [DCQ (DYU)-P];

s = [hw ; z_o];

t = residual_o;

condition = norm(residual_o);

hw_ite = hw;

zw_ite = z_o;

start_point = ...
    [
    hw;
    z_o;
    ];

yy = start_point;

% Perform iteration to obtain perturbations

while (condition > MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F,s,G,t);
    
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_h = y(1:N);
    delta_z = y(N+1:end);
    
    % Add structured perturbations to vector h.
    hw_ite = hw_ite + delta_h;
    
    % Add structured perturbations to vector z.
    zw_ite = zw_ite + delta_z;
    
    % Seperate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    zz = zw_ite;
    
    zi_ite = cell(1,length(n)+1);
    for i = 1:1:length(n)+1
        zi_ite{i} = zz(1:m(i)+1);
        zz(1:m(i)+1) = [];
    end
    
    % Renew Matrix Pz
    Pz = P*zw_ite;
    
    %Increment s in LSE Problem
    s = -(yy-start_point);
    
    % Copy vector h_ite
    hh = hw_ite;
    % Move individual vectors hi into variable size array, emptying
    % hh
    for i = 1:1:length(n)
        hh(1:n(i)+1);
        arr_hw{i} = hh(1:n(i)+1);
        hh(1:n(i)+1) = [];
    end
    
    %Build iterative DYU
    DYU = BuildDYU(arr_hw,m);
    
    % add the structured perturbations to improved fw array.
    for i = 1:length(fw)
        new_fw{i} = fw{i} + zi_ite{i};
    end
    
    %Build DCEQ
    DCEQ = BuildDCQ(new_fw) ;
    
    % Build G
    G = [DCEQ (DYU-P)];
    
    % Calculate residual and increment t in LSE Problem
    r = ((RHS_f+Pz) - (DCEQ*hw_ite));
    t = r;
    
    condition = norm(r)./norm(RHS_f+Pz);
    
    % Increment iteration number
    ite = ite + 1;
end

% Print outputs to command line
fprintf('Performed Deconvolutions...\n')
fprintf('Iterations required for Batch Deconvolution %i\n', ite)


end

function Y_new = BuildDYU(set_hw,m)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

for i = 1:1:length(set_hw)
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hw = set_hw{i};
    
    % Get degree of f_{i}
    deg_fw = m(i+1);
    
    y{i} = real(BuildD0Y1U1(hw,deg_fw));
end

%Build the Coefficient Matrix C
num_Rows = 0;
for i = 1:length(m)-1
    num_Rows = num_Rows + 1 + (m(i));
end
cols = (m(1)+1);

xx = zeros(num_Rows,cols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y_new = Y;


end

function Y1 = BuildD0Y1U1(hx,m1)
global BOOL_LOG

switch BOOL_LOG
    case 'y' % use logs
        Y1 = BuildD0Y1U1_log(hx,m1);
        
    case 'n' % use nchoosek
        Y1 = BuildD0Y1U1_nchoosek(hx,m1);
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

global BOOL_DENOM_SYL

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

switch BOOL_DENOM_SYL
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

global BOOL_DENOM_SYL

% Y1 = zeros(m0+1,m1+1);

% Get Degree of polynomial deg(h_{x}) = n_{1} = m_{0}-m_{1}
n1 = GetDegree(hx);

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
switch BOOL_DENOM_SYL
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



function opt_theta = getOptimalTheta(set_f,v_m)
% Get optimal value of theta for the matrix. 
%
%
% Inputs
%
% set_f : Set of vectors f_{i}
%
% v_m : vector which stores each m_{i}, the degree of the polynomial f_{i}
%
%


% Let f_{i} denote the ith polynomial stored in the array set_g
%


% Get number of polynomials in set_g
nPolys = length(set_f);

%For each coefficient ai,j
% Let \lambda_{i,j} be its max value in c_i(f_i)
% Let \mu_{i,j} be its min value in c_{i}(f_i)

F_max = cell(1,nPolys-1);
F_min = cell(1,nPolys-1);

% For each polynomial f_{1},...,f_{}, note we exclude f_{0} from this,
% since f_{0} does not appear in the LHS matrix.
for i = 1:1:nPolys-1  
    
    % Get polynomial f_{i} from the set g containing all f_{i}
    prev_fw = set_f{i+1}
    fw = set_f{i+1};
    
    deg_fw = GetDegree(fw);
    deg_prev_fw = GetDegree(prev_fw);
    
    deg_hw = deg_fw - deg_prev_fw;
    
    % Each coefficient a_{j} of polynomial f_{i} appears in deg(h_{i})
    % columns, where the degree of h_{i} = m(i-1) - m(i).
    
    % Assign empty vectors for the max and minimum values of each
    % coefficient in F.
    F_max{i} = zeros(1,v_m(i+1)+1);
    F_min{i} = zeros(1,v_m(i+1)+1);
    
    % For each coefficient a_{j} in f_{i+1}
    for j = 0:1:deg_hw
        
        % Get the coefficient a_{j} of polynomial f_{i}
        aij = fw(j+1);
        
        % initialise a vector to store all the a_{i}
        x = zeros(1,deg_hw+1);
        
        % For each occurence of the coefficient ai_j in the columns of C(fi)
        for k = 0:1: deg_hw
            x(k+1) = aij .* nchoosek(j+k,k) .* nchoosek(v_m(i)-(j+k),v_m(i+1)-j) ./ nchoosek(v_m(i),v_m(i+1));
            %x(k+1) = aij .* nchoosek(deg_fw,j+k) * nchoosek(deg_hw,k) ./ nchoosek(deg_fw + deg_hw,j);
            
        end
        
        % Get max entry of each coefficient.
        F_max{i}(j+1) = max(abs(x));
        % Get min entry of each coefficient.
        F_min{i}(j+1) = min(abs(x));
        
    end
end

opt_theta = MinMaxOpt(F_max,F_min,v_m);

end

function theta = MinMaxOpt(F_max,F_min,m,nPolys)
%
% This function computes the optimal value theta for the preprocessing
% opertation as part of block deconvolution
%
% F_max   :  A vector of length m1+1 + m2+1 + ... + md+1, such that F_max(i) stores the
%            element of maximum magnitude of D(C(f))Q that contains the
%            coefficient a(i,j) of polys fi, j=0,...,m1.
%
% F_min   :  A vector of length m+1, such that F_min(i) stores the
%            element of minimum magnitude of S(f,g)Q that contains the
%            coefficient a(i) of f, i=1,...,m+1.
%

% Get the number of polynomials
nPolys = size(F_max,2);



f = [1 -1 0];

Part1 = [];

% For each Ai
for i = 1:1:nPolys
    
    % Get the max of each coefficient of polynomial fw
    fw_max = F_max{i};
    % Get Degree of the polynomial 
    deg_fw = GetDegree(fw_max)
    
    % Build the matrix A_{i}
    Ai = [ones(deg_fw+1,1) zeros(deg_fw+1,1)   -(0:1:deg_fw)'];
    
    % Append the matrix
    Part1 = [Part1 ; Ai];
end


Part2 = [];
% For each Bi
for i = 1:1:nPolys
    
    % Get the max of each coefficient in polynomial f(x)
    fw_max = F_max{i}
    
    % Get the degree of the polynomial f(x)
    deg_fw = GetDegree(fw_max);
    
    Bi = [zeros(deg_fw+1,1) -ones(deg_fw+1,1) (0:1:deg_fw)'];
    
    Part2 = [Part2 ; Bi];
end

A = -[Part1;Part2];

% Get the array of entries F_max_{i} as a vector
v_F_max = cell2mat(F_max)'
v_F_min = cell2mat(F_min)'

b = [-log10(v_F_max); log10(v_F_min)];


% Solve the linear programming problem and extract alpha and theta
% from the solution vector x.
try
    x = linprog(f,A,b);
    theta = 10^x(3);
catch
    fprintf('Error Calculating Optimal value of theta\n');
    theta = 1;
    
end

end

function DCQ = BuildDCQ(set_fw)
% set fw is the cell array of poly coefficiencts fw_i
%
% Inputs.

% create cauchy matrices c{i} for i = 1,...
c = cell(1,length(set_fw));

for i = 1:1:length(set_fw)-1
    
    % Start with set_fw{2} = f_{1}(i)
    prev_fw = set_fw{i};
    fw = set_fw{i+1};
    
    % Get the degree of f, previous f, and h
    deg_fw = GetDegree(fw);
    deg_prev_fw = GetDegree(prev_fw);
    
    deg_hw = deg_prev_fw - deg_fw;
    
    c{i} = BuildD0C1Q1(fw,deg_hw);


end

%Build the Coefficient Matrix C of all matrices c
DCQ = blkdiag(c{1:length(c)});
DCQ = real(DCQ);

end


function D0C1Q1 = BuildD0C1Q1(fx,m0_m1)
%% Build the Partition D_{i-1}C_{i}Q_{i}
%
%
% Inputs
%
% fx -
%
% m0_m1 = m_{0} - m_{1}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise global variables.
global BOOL_LOG


switch BOOL_LOG
    case 'n' 
        % use nchoosek
        D0C1Q1 = BuildD0C1Q1_nchoosek(fx,m0_m1);
        
    case 'y' 
        % use log
        D0C1Q1 = BuildD0C1Q1_log(fx,m0_m1);
    otherwise 
        error('err')
end

end


function D0C1Q1 = BuildD0C1Q1_nchoosek(fx,m0_m1)
% Build a partition D_{i-1}C_{i}Q_{i} of the DCQ matrix
%
%
% Inputs
%
% fx : Coefficients of polynomial f(x)
%
% m0_m1 : Degree of polynomial h_{1}
%

global BOOL_DENOM_SYL

% Get Degree of polynomial f(x)
m1 = GetDegree(fx);

% Define n1 the degree of h1
n1 = m0_m1;

% Preassign space for partition D_{0}C_{1}Q_{1}
D0C1Q1 = zeros(n1+m1+1,n1+1);

% For each column k = 0:1:m_{i-1} - m_{i}
for k = 0:1:n1
    % For each row j = k:1:m_{i}+k
    for j = k:1:m1+k
        D0C1Q1(j+1,k+1) = ...
            fx(j-k+1) .* ...
            nchoosek(j,k) .* nchoosek(n1+m1-j,m1-(j-k)) ;
    end
end
switch BOOL_DENOM_SYL
    case 'y'
        denom = nchoosek(n1+m1,m1);
        D0C1Q1 = D0C1Q1 ./ denom ;
    case 'n'
    otherwise
        error('err')
end
end


function D0C1Q1 = BuildD0C1Q1_log(fx,m0_m1)
% Build a partition D_{i-1}C_{i}Q_{i} of the DCQ matrix
%
% Inputs.
%
% fx -
%
% m0 - Degree of f_{i-1}
%
% m1 - Degree of f_{i}
%

global BOOL_DENOM_SYL

% Define n_{1} the degree of h_{1}
n1 = m0_m1;

% Preassign space for partition D_{0}C_{1}Q_{1}
D0C1Q1 = zeros(m0+1,n1+1);

% For each column k = 0:1:m_{i-1} - m_{i}
for k = 0:1:n1
    % For each row j = k:1:m_{i}+k
    for j = k:1:m1+k
        Numerator_log = lnnchoosek(j,k) + lnnchoosek(n1+m1-j,m1-(j-k));
        Numerator_exp = 10.^Numerator_log;
        D0C1Q1(j+1,k+1) = ...
            fx(j-k+1) .* Numerator_exp ;
    end
end

switch BOOL_DENOM_SYL
    case 'y'
        Denom_Log = lnnchoosek(n1+m1,m1);
        Denom_exp = 10.^Denom_Log;
        
        D0C1Q1 = D0C1Q1./ Denom_exp;
    case 'n'
    otherwise
        error('err')
end
end


