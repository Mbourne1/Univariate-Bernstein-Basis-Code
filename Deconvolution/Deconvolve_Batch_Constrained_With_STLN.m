function [arr_hx] = Deconvolve_Batch_Constrained_With_STLN(arr_fx,vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : Array of polynomials f(x)
%
% vMult : Multiplicities of the factors of f_{0}(x) in ascending order.
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x)

% Global Variables
global SETTINGS

% Get the number of polynomials in the array of polynomials f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% Get the number of polynomials in the array of polynomials h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx -1 ;

% % Get the degree of each of the polynomials f{i}(x)

% Initialise vector
vDeg_arr_fx = zeros(nPolys_arr_fx,1);

% For each polynomial f_{i} get the degree
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get the degree n(i) of polynomials h_{i}.
vDeg_arr_hx = (vDeg_arr_fx(1:end-1) - vDeg_arr_fx(2:end));

% Define M to be the total number of all coefficients of the first d polynomials
% f_{0}...f_{d-1},
M = sum(vDeg_arr_fx+1) - (vDeg_arr_fx(end)+1);

% Define M1 to be the total number of all coefficients of polynomials
% f_{0},...,f_{d}
nCoefficients_fx = sum(vDeg_arr_fx+1);

% Define N to be the number of coefficients of all h_{i}
nCoefficients_hx = sum(vDeg_arr_hx+1);


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
        error('err');
end

% %
% %
% Initialise a cell-array for f(w)
% %
arr_fw = cell(nPolys_arr_fx,1);

% for each f_{i} get fw_{i}
for i = 1:1:nPolys_arr_fx
    arr_fw{i} = GetWithThetas(arr_fx{i},theta);
end

% %
% %
% Build LHS Matrix C(f1,...,fd)
DT_fwQ = BuildDTQ(arr_fw,vMult);

% %
% %
% Build the RHS vector

% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec_fw = BuildRHSF(arr_fw);



v_pw = SolveAx_b(DT_fwQ,RHS_vec_fw);

nCoefficients_px = length(v_pw);

unique_vMult = unique(vMult);

nPolys_arr_px = length(unique_vMult);

vDeg_arr_px = zeros(nPolys_arr_px,1);

for i = 1:1:length(unique_vMult)
    mult = unique_vMult(i);
    deg = vDeg_arr_fx(mult) - vDeg_arr_fx(mult+1);
    vDeg_arr_px(i) = deg;
end

% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).
arr_pw = GetArray(v_pw,vDeg_arr_px);

% %
% %
% Get the polynomials p_{i}(\omega) repeated to give the set of polynomials
% h_{i}(\omega)
arr_hw = Get_hx(arr_pw,unique_vMult);

% %
% %
% Build the array of polynomails z(\omega) which are the structured
% perturbations of the array of polynomials f(x).
arr_zw = cell(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    arr_zw{i} = zeros(vDeg_arr_fx(i) + 1,1);
end

% Build vector z(\omega) consisting of all vectors in z_{i}(x)
v_zw = cell2mat(arr_zw);

% Build the matrix P
P = [eye(M) zeros(M,nCoefficients_fx-M)];

% DY_hQ
DY_hQ = BuildDYQ(arr_hw,vDeg_arr_fx);

% Set iteration number
ite = 1;

% Build the matrix F
F = eye(nCoefficients_px + nCoefficients_fx);

%
%
%
% Build the matrix G

% Build component H_h of G
H_h = DT_fwQ;

% Build component H_z of G
H_z = DY_hQ - P;

G = [H_h H_z];

% Compute the first residual
res_vec = RHS_vec_fw + (P*v_zw) - (DT_fwQ * v_pw);

% Update Matrix P*z
Pz = P*v_zw;

% Perform test
for i = 1 : 1 : nPolys_arr_fx
    vec_fw = [RHS_vec_fw; arr_fx{i}];
end
%test1 = DY_hQ * vec_fw;
%test2 = DT_fwQ * v_pw;
%test1./test2

condition(ite) = norm(res_vec)./norm(RHS_vec_fw + Pz);

start_point = ...
    [
    v_pw;
    v_zw;
    ];

%Get the iterated value
yy = start_point;

% Get
s = -(yy - start_point);



while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    % output y gives delta p and delta z
    delta_pw = y(1:nCoefficients_px);
    delta_zw = y(nCoefficients_px+1:end);
    
    % Add structured perturbations to vector p(\omega)
    v_pw = v_pw + delta_pw;
    
    % Add structured perturbations to vector z(\omega)
    v_zw = v_zw + delta_zw;
    
    % Get the updated array of polynomials p_{i}(\omega)
    arr_pw = GetArray(v_pw, vDeg_arr_px);
    
    arr_zw = GetArray(v_zw, vDeg_arr_fx);
    
    arr_hw = Get_hx(arr_pw,unique_vMult);
    
    s = -(yy-start_point);
    
    DY_hQ = BuildDYQ(arr_hw,vDeg_arr_fx);
    
    % Build the matrix C(f)
    DT_fwQ = BuildDTQ(arr_fw,vMult);
    
    % Build the matrix C(z)
    DT_zwQ = BuildDTQ(arr_zw,vMult);
    
    % Build G
    H_z = DY_hQ - P;
    H_h = DT_fwQ + DT_zwQ;
    
    G = [H_h H_z];
    
    % Update the RHS vector
    RHS_vec_fw = BuildRHSF(arr_fw);
    RHS_vec_Pz = BuildRHSF(arr_zw);
    
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((RHS_vec_fw+RHS_vec_Pz) - ((DT_fwQ + DT_zwQ)*v_pw));
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    condition(ite) = norm(res_vec)./norm(RHS_vec_fw + RHS_vec_Pz);
end

% Get array of polynomials h_{i}(x) from h_{i}(\omega)
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1:nPolys_arr_hx
    arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
end

%
%
%
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Required Number of iterations : %i \n',ite)]);
LineBreakLarge();

end

function LHS_Matrix = BuildDTQ(arr_fx,vMult)
% %
% %
% Build the LHS Coefficient matrix

% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]

% Get number of distinct polynomials in h_{i}
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i > 1
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
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
        Tf{j} = BuildT1(fx, deg_hx);
        
        D{j} = BuildD_2Polys(deg_fx, deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        Cf{i} = [Cf{i} ; D{j}*Tf{j}];
    end
    
    Q{i} = BuildQ1(deg_hx);
    
    DTQ{i} = Cf{i} * Q{i};
    
end


LHS_Matrix = blkdiag(DTQ{:});

end



function arr_zx = GetArray(v_zx,v_degree_f)
% Given the vector of perturbations of f(x) given by v_zx

% Get number of polynomials in arr_fx
nPolys_fx = length(v_degree_f);

% Initialise an array
arr_zx = cell(nPolys_fx,1);

for i = 1:1:nPolys_fx
    
    % Get the m+1 coefficients from the vector
    arr_zx{i} = v_zx(1:v_degree_f(i)+1);
    % Remove the m+1 coefficients
    v_zx(1:v_degree_f(i)+1) = [];
    
end


end

function arr_hx = Get_hx(arr_px,vUniqueMult)

% Get number of entries in the array of polynomials p_{i}(x)
nEntries_arr_px = size(arr_px,1);

% initialise count
count = 1;


for i = 1:1:nEntries_arr_px
    
    if i == 1
        nReps = vUniqueMult(i);
    else
        nReps = (vUniqueMult(i) - vUniqueMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hx{count,1} = arr_px{i};
        count = count + 1;
    end
    
end
end


function RHS_vec_fx = BuildRHSF(arr_fx)

% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% vDeg_arr_fx = zeros(nPolys_arr_fx,1);
%
% for i = 1:nPolys_arr_fx
%     vDeg_arr_fx(i) = GetDegree(arr_fx{i});
% end
%
% % Get the total number of coefficients in the first d polynomials
% % f_{0},...f_{d-1}
% nCoefficients_RHS = (sum(vDeg_arr_fx) + nPolys_arr_fx) - (vDeg_arr_fx(end)+1);

% Initialise an empty RHS Vector
RHS_vec_fx = [];

for i = 1:1:nPolys_arr_fx - 1;
    RHS_vec_fx = [RHS_vec_fx ; arr_fx{i}];
end

end


function Y_new = BuildDYQ(arr_hw,vDeg_arr_fx)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

% Inputs.
%
% arr_hw : Set of polynomials h_{i}(w)
%
% vDeg_arr_fx : vector of degrees of polynomials f_{0},...

nPolys_hw = size(arr_hw,1);

for i = 1:1:nPolys_hw
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hw = arr_hw{i,1};
    
    % Get degree of f_{i}
    deg_fw = vDeg_arr_fx(i+1);
    
    y{i,1} = real(BuildD0Y1U1(hw,deg_fw));
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
% hx : Coefficients of polynomial h_{i}(x)
%
% m1 : Degree of polynomial f_{i}



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


function Y1 = BuildD0Y1U1_log(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
% m0 = degree of previous polynomial
% m1 = degree of current polynomial


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

% Include the denominator
Denom_Log = lnnchoosek(n1+m1,m1);
Denom_Exp = 10.^Denom_Log;

Y1 = Y1 ./  Denom_Exp;


end
