function arr_hx = Deconvolve_Batch(arr_fx)
% Perform batch deconvolution

% Get the number of polynomials in the set set_f
nPolys_arr_fx = size(arr_fx,1);

% let d be the number of deconvolutions = num of polynomials in set_f - 1
nPolys_arr_hx = nPolys_arr_fx - 1;

% %
% %
% Get the degree m_{i} of each of the polynomials f_{i}

% Intialise vector to store the degree of the set of polynomaisl f_{i}
vDeg_arr_fx = zeros(1,nPolys_arr_fx);

% For each polynomial f_{i}, get its degree.
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% %
% %
% %
% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
vDeg_arr_hx = (vDeg_arr_fx(1:end-1) - vDeg_arr_fx(2:end))'; 


% %
% % Preprocessing
% %
% Obtain theta such that the ratio of max element to min element is
% minimised
%theta = GetOptimalTheta(arr_fx,vDeg_arr_fx);
theta = 1;

% % 
% %
% Initialise a cell-array for f(w)
% %
arr_fw = cell(nPolys_arr_fx,1);

% for each f_{i} get fw_{i}
for i = 1:1:length(arr_fx)
    arr_fw{i,1} = GetWithThetas(arr_fx{i},theta);
end

% % 
% %
% %
% Write Deconvolutions in form [D^{-1}C(f)Q] h = RHS_f
RHS_vec = real(BuildRHSF(arr_fw));
DCQ = BuildDCQ(arr_fw);

% Get the solution vector h(w) in the system of equations
% DCQ * hw = RHS_vec.
v_hw = SolveAx_b(DCQ,RHS_vec);


% %
% %
% %
% Seperate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}
% initialise a cell array to store the coefficients of the individual
% polynomials h_{i}

% Split vec h in to an array of polynomials.

% I
arr_hw = GetArray(v_hw,vDeg_arr_hx);
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1:nPolys_arr_hx
   arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
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


function DCQ = BuildDCQ(set_fx)
% set fw is the cell array of poly coefficiencts fw_i
%
% Inputs.

% For each of the polynomials f_{i}(x), excluding the final polynomial
for i = 1:1:length(set_fx)-1
    
    % Get the polynomial f_{i} = set_f{i+1} 
    fw = set_fx{i+1};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fw_prev = set_fx{i};
    
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
    DT1Q1{i}  = D*T1*Q1;
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
    deg_hw = vDeg_arr_hx(i);
    
    % Get coefficients of h_{i} from the solution vector
    arr_hx{i} = v_hw(1:deg_hw+1);
    
    % Remove the coefficients from the solution vector
    v_hw(1:deg_hw+1) = [];
end

end
