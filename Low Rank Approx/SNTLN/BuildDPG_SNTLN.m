function DPG = BuildDPG_SNTLN(m, n, k, alpha, theta, idx_col)
% BuildDPQ(m,n,theta,mincol,t)
%
% Build the matrix DP. Build the matrix DP such that DP * [f;g] gives the
% column of the Sylvester subresultant matrix matrix whose index is given
% by idx_col.
%
% Used in SNTLN.m
%
%
% Inputs
%
% [m, n] : Degree of polynomial f(x) and g(x)
%
% alpha : Optimal value of \alpha
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
D = BuildD_2Polys(m,n-k);

% Build the matrices P_{1} and P_{2}

if idx_col <= nCols_Tf % Column is in first partition T_{n-k}(f) of S_{k}
    
    % Build the matrix P
    P1 = BuildP1(m,n-k,idx_col);
    
    P2 = zeros(m+n-k+1,n+1);
    
else  %  The column is from the second partiton of the Sylvester matrix poly g
    
    % Get index of column relative to second partition T_{m-k}(g)
    
    idx_col_rel = idx_col - (n-k+1);
    
    % Build the matrix P_{1}
    P1 = zeros(m+n-k+1,m+1);
    
    % Build the matrix P_{2}
    P2 =  BuildP1(n,m-k,idx_col_rel);
    
end

% Build the matrices Q
Q1 = BuildQ1(m);
Q2 = BuildQ1(n);

% Get thetas associated with polynomial f(x) and g(x)
th_f = diag(theta.^(0:1:m));
th_g = diag(theta.^(0:1:n));

DPG = D*[P1*Q1*th_f alpha.*P2*Q2*th_g];

end

