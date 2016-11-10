function DPQ = BuildDP_SNTLN(m,n,alpha,theta,idx_col,k)
% BuildDP(m,n,theta,mincol,t)
%
% Build the matrix DP. Build the matrix DP such that DP * [f;g] gives the
% column of the Sylvester subresultant matrix matrix whose index is given
% by idx_col.
%
% Used in SNTLN.m
% Used in STLN.m
%
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% theta : Optimal value of \theta
%
% idx_col : Index of column c_{k} removed from S_{k}(f,g)
%
% k : Degree of GCD d(x)
%
% Outputs.
%
% DPQ : matrix DPQ

% Get the number of columns in T_{n-k}(f)
nCols_Tf = n-k+1;


% Build the matrix D^{-1}_{m+n-k}
D = BuildD(m,n-k);
    
% %
% Build the matrix P = [P1 P2]

if idx_col <= nCols_Tf % Column is in first partition T_{n-k}(f) of S_{k}
        
    % Build the matrix P
    P1 = BuildP1(m,n-k,idx_col,theta);
    
    P2 = zeros(m+n-k+1,n+1);
    
else  %  The column is from the second partiton of the Sylvester matrix poly g
    
    % Get index of column relative to second partition T_{m-k}(g)
    
    idx_col_rel = idx_col - (n-k+1);
    
    % Build the matrix P
    P1 = zeros(m+n-k+1,m+1);
    P2 = alpha .* BuildP1(n,m-k,idx_col_rel,theta);
    
end

% Build the matrix Q
Q1 = BuildQ1(m);
Q2 = BuildQ1(n);
Q = blkdiag(Q1,Q2);

DPQ = D*[P1 P2] * Q;

end


function P1 = BuildP1(m,n_k,idx_col,theta)
% Build the matrix P1, a partition of P, where P*[f;g] gives a column of
% the Sylvester subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n_k : Degree of polynomial v(x)
%
% idx_col : Index of column of T_{n-k}(f) which is being constructed.

% Split P1 into three parts [P1_a ; P1 b : P1_c]

P1 = zeros(m+n_k+1,m+1);

%
theta_mat = diag(theta.^(0:1:m));

%
temp_mat = theta_mat .* nchoosek(n_k,idx_col-1);


P1(idx_col:idx_col+m,:) = temp_mat;

end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function G = LeftDG_Build(m,n,t,idx_col,theta)
% Build the matrix G, where G forms part of DP.
%
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Degree of GCD
%
% idx_col : Index of column of S_{k} being constructed by P*[f;g]
%
% theta : Optimal value of \theta_{1}
%
% Outputs.
%
% G : Matrix G


% Global Variables
global SETTINGS



switch SETTINGS.BOOL_LOG
    
    case 'y' % use log
        G = LeftDG_Build_log(m,n,t,idx_col,theta);
        
    case 'n' % use nchoosek
        G = LeftDG_Build_nchoosek(m,n,t,idx_col,theta);
        
    otherwise
        error('SETTINGS.BOOL_LOG must be either y or n')
end
end

function G = RightDG_Build(n,m,t,mincol,theta)
%
%
%
% Inputs.
%
% m : Degree of polynomial f
%
% n : Degree of polynomial g
%
% t : Degree of GCD
%
% mincol : index of column of S_{t}
%
% theta : Optimal value of theta


% Global Variables.
global SETTINGS

switch SETTINGS.BOOL_LOG
    case 'y' % use logs
        G = RightDG_Build_log(n,m,t,mincol,theta);
    case 'n' % use nchoosek
        warning('off','all');
        G = RightDG_Build_nchoosek(n,m,t,mincol,theta);
        warning('on','all');
    otherwise
        error('SETTINGS.BOOL_LOG must be either y or n')
end
end

function [G] = LeftDG_Build_nchoosek(m,n,t,num_mincolumn,theta)
%
%
% Inputs
%
% m : Degree of polynomial f
%
% n : Degree of polynomial g
%
% t : Degree of GCD
%
% mincol :
%
% theta :


% Global Variables
global SETTINGS

% Initalise vector G
G = zeros(1,m+1);

% Get column to be removed.
j = num_mincolumn-1;

for i = 0:1:m
    G(i+1) =  (theta.^(i)) .* nchoosek(i+j,j) .* nchoosek(m+n-t-(i+j),n-t-j);
end



% Create diagonal matrix G from values in vector G.
G = diag(G);


switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Included denominator in coefficient matrix.
        G = G./ nchoosek(m+n-t,n-t);
    case 'n'
        % Exclude
    otherwise
        error(' SETTINGS.BOOL_DENOM_SYL must be either y or n')
end

end

function [G] = LeftDG_Build_log(m,n,t,idxMinCol,theta)
% Get G \in\mathbb{R}^{(m+1)\times(m+1)}, where G is a diagonal matrix.
%
% Inputs
%
% m - Degree of polynomial f(x)
%
% n - Degree of polynomial g(x)
%
% t - Degree of GCD d(x)
%
% idxMinCol - Index of column removed from Sylvester Matrix
%
% theta - Optimal value of \theta

% Global Variables
global SETTINGS

% Initialise empty vector G.
G = zeros(1,m+1);

j = idxMinCol - 1;

for i = 0:1:m
    
    % Evaluate binomial coefficients in the numerator in logs
    Numerator_Eval_log = ...
        lnnchoosek(i+j,j) + ...
        lnnchoosek(m+n-t-(i+j),n-t-j);
    
    % Convert to standard form
    Numerator_Eval_exp = 10.^Numerator_Eval_log;
    
    %
    G(i+1) =  (theta.^(i)) .* Numerator_Eval_exp;
    
end

% Create diagonal matrix G from values in vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Included denominator in coefficient matrix.
        
        % Evaluate binomial coefficient in denominator in logs
        Denom_Eval_log = lnnchoosek(m+n-t,n-t);
        
        % Conver to standard form
        Denom_Eval_exp = 10.^Denom_Eval_log;
        
        % Divide G by the denominator.
        G = G./Denom_Eval_exp;
    case 'n'
        % Exclude denominator from coefficient matrix
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
        
end




end


function [G] = RightDG_Build_nchoosek(n,m,t,idxMinCol,theta)
%
%
% Inputs
%
% m - Degree of polynomial f(x)
%
% n - Degree of polynomial g(x)
%
% t - Degree of GCD
%
% idxMinCol -
%
% theta -



global SETTINGS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty vector G.
G = zeros(1,n+1);

j = idxMinCol-n+t-2;

for i=0:1:n
    
    G(i+1) = (theta.^(i)) .* nchoosek(i+j,j) .* nchoosek(m+n-t-(i+j),m-t-j) ;
    
end

% Create diagonal matrix G from vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Denominator is included in coefficient matrix.
        G = G./ nchoosek(m+n-t,m-t);
    case 'n'
        % Denominator is excluded from the coefficient matrix
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
        
end



end




function [G] = RightDG_Build_log(n,m,t,idxMinCol,theta)
%
%
% Inputs
%
% m - Degree of polynomial f
%
% n - Degree of polynomial g
%
% t - Degree of GCD
%
% mincol -
%
% theta -

% Global Variables
global SETTINGS

% Initialise the Right Matrix DG
G = zeros(1,n+1);


j = idxMinCol-n+t-2;


for i=0:1:n
    Numerator_Eval_log =...
        lnnchoosek(i+j,i) + ...
        lnnchoosek(m+n-t-(i+j),m-t-j);
    
    Numerator_Eval_exp = 10.^Numerator_Eval_log;
    
    G(i+1) = theta.^(i) .* Numerator_Eval_exp ;
    
end

% Create diagonal matrix G from vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Include the denominator
        
        Denom_log = lnnchoosek(m+n-t,n);
        
        Denom_exp = 10.^Denom_log;
        
        G = G./ Denom_exp;
    case 'n'
        % Exclude the denominator
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
end

end
