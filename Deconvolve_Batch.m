function hi = Deconvolve_Batch(set_g)
%  Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two deconvolutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                             Input:

% set_g - set of input polynomials g(y) to be deconvolved. Each g_{i} has a
% different number of elements, so set_g is a cell array.

% %                             Output:

% h_{i} = g_{i-1}/g_{i}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Global Variables

global max_iterations_deconvolutions
global max_error_deconvolutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% let d be the number of deconvolutions = num of polynomials in set_g - 1
d = length(set_g) -1;

% get values of m_{i}, the degrees of the polynomials f_{i}(y) for i =
% 0...d
m = zeros(1,length(set_g));
for i = 1:1:length(set_g)
    m(i) = length(set_g{i})-1;
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
theta = getOptimalTheta(set_g,m);


% for polynomial f_i
fw = cell(1,length(set_g));
for i = 1:1:length(set_g)
    fy = set_g{i};
    % Convert to scaled bernstein basis, by inclusion of theta
    % a_{i}\theta^{i}
    fw{i} = fy*theta^(i);
end

% Write Deconvolutions in form [D^{-1}C(f)Q] h = RHS_f

RHS_f = real(BuildRHSF(fw));
DCQ = BuildDCQ(fw,m);


% set bk the right hand side vector to be equal to RHS_f.
bk = RHS_f;

[~,n2] = size(DCQ);
[Q1,R] = qr(DCQ);

R1 = R(1:n2,:);
cd = Q1'*bk;
c = cd(1:n2,:);
x_ls = R1\c;

%r = norm(RHS_f-A*x_ls)
h_o = x_ls;
vec_h = h_o;


% Seperate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}
hi = cell(1,length(n));
for i = 1:1:length(n)
    hi{i} = h_o(1:n(i)+1);
    h_o(1:n(i)+1) = [];
end

h_o = vec_h;

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
DYU = BuildDYU(hi,m);

%
DCQh = DCQ*vec_h;

% First Residual
residual_o = (RHS_f+(P*z_o) - DCQh);

% Set the iteration counter.
iteration_num=0;

F = eye(N+M1);

G = [DCQ (DYU)-P];

s = [h_o ; z_o];

t = residual_o;

condition = norm(residual_o);

h_ite = h_o;

z_ite = z_o;

% Perform iteration to obtain perturbations

while condition >(max_error_deconvolutions) &&  iteration_num < max_iterations_deconvolutions
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    cond(G);
    y =LSE(F,s,G,t);
    
    % Output y gives delta h and delta z
    delta_h = y(1:N);
    delta_z = y(N+1:end);
    
    
    % Add structured perturbations to vector h.
    h_ite = h_ite + delta_h;
    
    % Add structured perturbations to vector z.
    z_ite = z_ite + delta_z;
    
    % Seperate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    zz = z_ite;
    
    zi_ite = cell(1,length(n)+1);
    for i = 1:1:length(n)+1
        zi_ite{i} = zz(1:m(i)+1);
        zz(1:m(i)+1) = [];
    end
    
    % Renew Matrix Pz
    Pz = P*z_ite;
    
    %Increment s in LSE Problem
    s = [-(h_ite-h_o); -z_ite];
    
    % Copy vector h_ite
    hh = h_ite;
    % Move individual vectors hi into variable size array, emptying
    % hh
    for i = 1:1:length(n)
        hh(1:n(i)+1);
        hi{i} = hh(1:n(i)+1);
        hh(1:n(i)+1) = [];
    end
    
    
    %Build iterative DYU
    DYU = BuildDYU(hi,m);
    
    
    % add the structured perturbations to improved fw array.
    for i = 1:length(fw)
        new_fw{i} = fw{i} + zi_ite{i};
    end
    
    %Build DCEQ
    DCEQ = BuildDCQ(new_fw,m) ;
    
    % Build G
    G = [DCEQ (DYU-P)];
    
    % Calculate residual and increment t in LSE Problem
    r = ((RHS_f+Pz) - (DCEQ*h_ite));
    t = r;
    
    condition = norm(r)./norm(RHS_f+Pz);
    
    % Increment iteration number
    iteration_num = iteration_num + 1;
end

% Print outputs to command line
fprintf('Performed Deconvolutions...\n')
fprintf('Iterations required for Batch Deconvolution %i\n', iteration_num)


end

function Y_new = BuildDYU(set_hw,m)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

for i = 2:1:length(set_hw)+1
    set_hw{i-1};
    y{i-1} = real(BuildD0Y1U1(set_hw{i-1},m(i-1),m(i)));
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

Y_new = [Y];


end

function Y1 = BuildD0Y1U1(hx,m0,m1)
global bool_log

switch bool_log
    case 'y' % use logs
        Y1 = BuildD0Y1U1_log(hx,m0,m1);
        
    case 'n' % use nchoosek
        Y1 = BuildD0Y1U1_nchoosek(hx,m0,m1);
end
end

function Y1 = BuildD0Y1U1_nchoosek(hx,m0,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% m0 = degree of previous polynomial
% m1 = degree of current polynomial

global bool_denom_syl

% Y1 = zeros(m0+1,m1+1);
Y1 = [];
% for each column j = 1:1:m0-m1+1
for j = 0:1:m1
    % for each row i = 1:1:m1+1
    for i = j:1:j+(m0-m1)
        Y1(i+1,j+1) = ...
            hx(i-j+1) .* nchoosek(i,j) .* nchoosek(m0-i,m1-j);
        
    end
end

switch bool_denom_syl
    case 'y'
        Y1 = Y1 ./  nchoosek(m0,m1);
end
end


function Y1 = BuildD0Y1U1_log(hx,m0,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
% m0 = degree of previous polynomial
% m1 = degree of current polynomial

global bool_denom_syl

% Y1 = zeros(m0+1,m1+1);

Y1 = [];
% for each column i = 1:1:m0-m1+1
for k = 0:1:m1
    % for each row j = 1:1:m1+1
    for j = k:1:k+(m0-m1)
        BinomsEval_Log = lnnchoosek(j,k) + lnnchoosek(m0-j,m1-(j-k));
        BinomsEval_Exp = 10.^BinomsEval_Log;
        Y1(j+1,k+1) = hx(j-k+1) .* BinomsEval_Exp;
    end
end
switch bool_denom_syl
    case 'y'
        % Include the denominator
        Denom_Log = lnnchoosek(m0,m1);
        Denom_Exp = 10.^Denom_Log;
        
        Y1 = Y1 ./  Denom_Exp;
    case 'n'
        % Exclude the denominator
end


end



function f = BuildRHSF(fw_array)
% Build the vector f such that it contains the elements of
% Rhs f = [f_{0},...,f_{n-1}]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fw = array of vectors f_{0},...,f_{n}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty vector.
f = [];

% for each vector f f_{0},...,f_{n-1} in fw_array, add to right hand
% side vector
for i=1:1:length(fw_array)-1
    f = [f;fw_array{i}];
end

end



function opt_theta = getOptimalTheta(set_g,m)
% Get optimal value of theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                           Inputs

% set_g - set of vectors

% m     - set of degrees, of each polynomial stored in the set_g array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let f_{i} denote the ith polynomial stored in the array set_g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get number of entries in set_g
d = length(set_g)-1;

%F_max = zeros(1,sum(m(2:end))+d);
%F_min = zeros(1,sum(m(2:end))+d);

F_max = {};
F_min = {};



%For each coefficient ai,j
% Let \lambda_{i,j} be its max value in c_i(f_i)
% Let \mu_{i,j} be its min value in c_{i}(f_i)

for i = 1:1:d % For each polynomial fi except the final f_{i}
    
    % Get polynomial f_{i} from the set g containing all f_{i}
    fi = set_g{i+1};
    
    % Assign empty vectors for the max and minimum values of each
    % coefficient in F.
    F_max{i} = zeros(1,m(i)-m(i+1));
    F_min{i} = zeros(1,m(i)-m(i+1));
    
    
    for j = 0:1:m(i+1)
        
        % For each coefficient ai_j
        aij = fi(j+1);
        x = zeros(1,m(i)-m(i+1)+1);
        
        % For each occurence of the coefficient ai_j in the columns of C(fi)
        for k = 0:1: m(i)-m(i+1)
            x(k+1) = aij .* nchoosek(j+k,k) .* nchoosek(m(i)-(j+k),m(i+1)-j) ./ nchoosek(m(i),m(i+1));
            
        end
        
        % Get max entry of each coefficient.
        F_max{i}(j+1) = max(abs(x));
        % Get min entry of each coefficient.
        F_min{i}(j+1) = min(abs(x));
        
    end
end

F_max = horzcat(F_max{:})';
F_min = horzcat(F_min{:})';

opt_theta = MinMaxOpt(F_max,F_min,m,d);

end

function theta = MinMaxOpt(F_max,F_min,m,d)

% This function computes the optimal value theta for the preprocessing
% opertation as part of block deconvolution

% F_max   :  A vector of length m1+1 + m2+1 + ... + md+1, such that F_max(i) stores the
%            element of maximum magnitude of D(C(f))Q that contains the
%            coefficient a(i,j) of polys fi, j=0,...,m1.

% F_min   :  A vector of length m+1, such that F_min(i) stores the
%            element of minimum magnitude of S(f,g)Q that contains the
%            coefficient a(i) of f, i=1,...,m+1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [1 -1 0];

Part1 = [];
% For each Ai
for i = 1:1:d
    % Build Ai
    mi = m(i+1);
    Ai = [ones(mi+1,1) zeros(mi+1,1)   -(0:1:mi)'];
    Part1 = [Part1 ; Ai];
end


Part2 = [];
% For each Bi
for i = 1:1:d
    mi = m(i+1);
    Bi = [zeros(mi+1,1) -ones(mi+1,1) (0:1:mi)'];
    Part2 = [Part2 ; Bi];
end

A = -[Part1;Part2];

b = [-log10(F_max); log10(F_min)];


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

function DCQ = BuildDCQ(set_fw,m)
% set fw is the cell array of poly coefficiencts fw_i
% m is the vector of degrees of polynomials fw_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set_fw - set of polynomials f_{i}
% m
% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create cauchy matrices c{i} for i = 1,...
c = cell(1,length(set_fw));
for i = 1:1:length(set_fw)-1
    c{i} = BuildD0C1Q1(set_fw{i+1},m(i),m(i+1));
end

%Build the Coefficient Matrix C of all matrices c
DCQ = blkdiag(c{1:length(c)});
DCQ = real(DCQ);

end


function D0C1Q1 = BuildD0C1Q1(fx,m0,m1)
%% Build the Partition D_{i-1}C_{i}Q_{i}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Inputs

% fx -

% m0 -

% m1 -

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise global variables.
global bool_log


switch bool_log
    case 'n' 
        % use nchoosek
        D0C1Q1 = BuildD0C1Q1_nchoosek(fx,m0,m1);
        
    case 'y' 
        % use log
        D0C1Q1 = BuildD0C1Q1_log(fx,m0,m1);
end

end


function D0C1Q1 = BuildD0C1Q1_nchoosek(fx,m0,m1)
% Build a partition D_{i-1}C_{i}Q_{i} of the DCQ matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% fx -

% m0 -

% m1 -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bool_denom_syl

% Define n1 the degree of h1
n1 = m0-m1;

% Preassign space for partition D_{0}C_{1}Q_{1}
D0C1Q1 = zeros(m0+1,n1+1);

% For each column k = 0:1:m_{i-1} - m_{i}
for k = 0:1:m0-m1
    % For each row j = k:1:m_{i}+k
    for j = k:1:m1+k
        D0C1Q1(j+1,k+1) = ...
            fx(j-k+1) .* ...
            nchoosek(j,k) .* nchoosek(m0-j,m1-(j-k)) ;
    end
end
switch bool_denom_syl
    case 'y'
        denom = nchoosek(m0,m1);
        D0C1Q1 = D0C1Q1 ./ denom ;
end
end


function D0C1Q1 = BuildD0C1Q1_log(fx,m0,m1)
% Build a partition D_{i-1}C_{i}Q_{i} of the DCQ matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fx -
% m0 - Degree of f_{i-1}
% m1 - Degree of f_{i}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bool_denom_syl

% Defind n1 the degree of h
n1 = m0-m1;

% Preassign space for partition D_{0}C_{1}Q_{1}
D0C1Q1 = zeros(m0+1,n1+1);

% For each column k = 0:1:m_{i-1} - m_{i}
for k = 0:1:m0-m1
    % For each row j = k:1:m_{i}+k
    for j = k:1:m1+k
        Numerator_log = lnnchoosek(j,k) + lnnchoosek(m0-j,m1-(j-k));
        Numerator_exp = 10.^Numerator_log;
        D0C1Q1(j+1,k+1) = ...
            fx(j-k+1) .* Numerator_exp ;
    end
end

switch bool_denom_syl
    case 'y'
        Denom_Log = lnnchoosek(m0,m1);
        Denom_exp = 10.^Denom_Log;
        
        D0C1Q1 = D0C1Q1./ Denom_exp;
end
end


