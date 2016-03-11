function HYQ = BuildHYQ(dx,m,n,theta)


    HYQ_v1 = BuildHYQ1(dx,m,n,theta);
    HYQ_v2 = BuildHYQ2(dx,m,n,theta);
    
    HYQ = HYQ_v1;
    
    
end

function HYQ = BuildHYQ1(dx,m,n,theta)

% Builds the Matrix H.Y.Q where Y_{k} performs a change of variable.
% E(w).d = Y(d).w where d is vector of coefficients of GCD, and w is vector
% of perturbations to f and g

% (See Report "APF - General Case - Rearrangement of the matrix vector
% Product")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% dw -

% m -

% n -

% theta -

% alpha - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables

% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
global BOOL_LOG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch BOOL_LOG
    case 'n' % Use Nchoosek
        A = BuildHYPartition_nchoosek(dx,m,theta);
        B = BuildHYPartition_nchoosek(dx,n,theta);
        
    case 'y' % use logs
        A = BuildHYPartition_log(dx,m,theta);
        B = BuildHYPartition_log(dx,n,theta);
    otherwise 
        error('err')
end

HYQ = blkdiag(A,B);

end

function A = BuildHYPartition_nchoosek(dx,m,theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs.

% dw -

% m

% theta -

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global BOOL_DENOM_APF

t = length(dx)-1;
A = zeros(m+1,m-t+1);

% for each column j
for j = 0:1:(m-t)
    % for each row i
    for i = j:1:t+j
        %A(i+1,j+1) = dx(i-j+1).* theta^(j) .* nchoosek(i,j) .* ...
        
        % EDIT 03/06/2015 - .
        % Change to power of theta - in the theory was theta^i-j 
        A(i+1,j+1) = dx(i-j+1).* theta^(i) ...
            .* nchoosek(i,j) ...]
            .* nchoosek(m-i,t-(i-j)) ;

    end
end

switch BOOL_DENOM_APF
    case 'y'
        % include the denominator
        A = A./nchoosek(m,t);
    case 'n'
        A = A;
    otherwise 
        error(err)
end
end

function A = BuildHYPartition_log(dx,m,theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% dw -

% m

% theta -

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                         Outputs

% A : Matrix forming a partition of HY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                   Global Variables

% bool_denom_apf - (Boolean) Given the rearrangement of the Sylvester 
% matrix in the Bernstein basis, each partition of each subresultant has a 
% common divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.

global BOOL_DENOM_APF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree of polynomial dw
t = length(dx)-1;

% Initialise an empty matrix
A = zeros(m+1,m-t+1);

% for each column j of the matrix A
for j = 0:1:(m-t)
    % for each nonzero coefficient.
    for i = j:1:t+j
        
        Num_eval_log = ...
            lnnchoosek(i,j) +...
            lnnchoosek(m-i,t-(i-j));
        
        Num_eval_exp = 10.^Num_eval_log;
        
        A(i+1,j+1) = dx(i-j+1) .* theta^(i) .* Num_eval_exp;
        %A(i+1,j+1) = dx(i-j+1) .* theta^(j) .* Num_eval_exp;
    end
end


switch BOOL_DENOM_APF
    case 'y'
        Denom_eval_log = lnnchoosek(m,t);
        Denom_eval_exp = 10.^Denom_eval_log;
        A = A ./ Denom_eval_exp;
    case 'n'
    otherwise
        error(err)
end


end



function HYQ = BuildHYQ2(dx,m,n,theta)

% Get the degree of the gcd
t = length(dx)-1;

% Build Diagonal Matrix H
H = BuildH(m,n);

% Build the matrix Y = [Y1 zeros ; zeros Y2]
Y1 = zeros(m+1,m-t+1);
% For each column j
for j = 0:1:(m-t)
    % for each row i
    for i = j:1:t+j
        Y1(i+1,j+1) = dx(i-j+1) .* nchoosek(t,i-j) .* theta^i;
    end
end

Y2 = zeros(n+1,n-t+1);
for j = 0:1:(n-t)
    % for each row i
    for i = j:1:t+j
        Y2(i+1,j+1) = dx(i-j+1) .* nchoosek(t,i-j) .* theta^i;
    end
end

Y = blkdiag(Y1,Y2);

% Build Matrix Q

Q = BuildQ(m,n,t);

HYQ = H*Y*Q;
end


