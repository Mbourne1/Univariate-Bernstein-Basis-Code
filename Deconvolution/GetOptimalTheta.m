function theta = GetOptimalTheta(arr_fx)
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


% Get number of polynomials in the array
nPolys_arr_fx = size(arr_fx,1);

%For each coefficient ai,j
% Let \lambda_{i,j} be its max value in c_i(f_i)
% Let \mu_{i,j} be its min value in c_{i}(f_i)

F_max = cell(nPolys_arr_fx,1);
F_min = cell(nPolys_arr_fx,1);

% Get vector of degrees of polynomials f_{i}(x)
vDeg_arr_fx = zeros(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get vector of degrees of polynomials h_{i}(x)
vDeg_arr_hx = vDeg_arr_fx(1:end-1) - vDeg_arr_fx(2:end);


% For each polynomial f_{1},...,f_{d}, note we exclude f_{0} from this,
% since f_{0} does not appear in the LHS matrix.
for i = 2:1:nPolys_arr_fx  
    
    % Get polynomial f_{i} from the set g containing all f_{i}
    fx = arr_fx{i};
    
    % Assign empty vectors for the max and minimum values of each
    % coefficient in F.
    
    m = vDeg_arr_fx(i);
    n = vDeg_arr_hx(i-1);
    
    F_max{i} = zeros(m+1,1);
    F_min{i} = zeros(m+1,1);
    
    % For each coefficient a_{j} in f_{i+1}
    for j = 0:1:m
        
        % Get the coefficient a_{j} of polynomial f_{i}
        aij = fx(j+1);
        
        % initialise a vector to store all the a_{i}
        x = zeros(n+1,1);
        
        % For each occurence of the coefficient ai_j in the columns of C_{n_{i}}(f_{i})
        for k = 0:1:n
            %x(k+1) = aij .* nchoosek(j+k,k) .* nchoosek(v_m(i)-(j+k),v_m(i+1)-j) ./ nchoosek(v_m(i),v_m(i+1));
            %x(k+1) = aij .* nchoosek(deg_fw,j+k) * nchoosek(deg_hw,k) ./ nchoosek(deg_fw + deg_hw,j);
            x(k+1) = aij .* nchoosek(m,j) * nchoosek(n,k) ./ nchoosek(m+n,j+k);
        end
        
        % Get max entry of each coefficient.
        F_max{i}(j+1) = max(abs(x));
        
        % Get min entry of each coefficient.
        F_min{i}(j+1) = min(abs(x));
        
    end
end

theta = MinMaxOpt(F_max,F_min);

end

function theta = MinMaxOpt(F_max,F_min)
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
nPolys = size(F_max,1);



f = [1 -1 0];


% For each Ai
for i = 1:1:nPolys
    
    % Get the max of each coefficient of polynomial fw
    fw_max = F_max{i,1};
    % Get Degree of the polynomial 
    deg_fw = GetDegree(fw_max);
    
    % Build the matrix A_{i}
    Ai{i,1} = [ones(deg_fw+1,1) zeros(deg_fw+1,1)   -(0:1:deg_fw)'];
    
end



% For each Bi
for i = 1:1:nPolys
    
    % Get the max of each coefficient in polynomial f(x)
    fw_max = F_max{i,1};
    
    % Get the degree of the polynomial f(x)
    deg_fw = GetDegree(fw_max);
    
    Bi{i,1} = [zeros(deg_fw+1,1) -ones(deg_fw+1,1) (0:1:deg_fw)'];
    
end

Part1 = cell2mat(Ai);

Part2 = cell2mat(Bi);

A = ...
    [
    Part1;
    Part2
    ];

% Get the array of entries F_max_{i} as a vector
v_F_max = cell2mat(F_max);
v_F_min = cell2mat(F_min);

b = [log10(v_F_max); -log10(v_F_min)];


% Solve the linear programming problem and extract alpha and theta
% from the solution vector x.
try
    x = linprog(f,-A,-b);
    theta = 10^x(3);
catch
    fprintf('Error Calculating Optimal value of theta\n');
    theta = 1;
    
end

end