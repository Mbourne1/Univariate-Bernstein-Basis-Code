function DPQ = BuildDPG_STLN(m,n,k,idx_col)
% BuildDPQ(m,n,k,idx_col)
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

% Build DPG
DPQ = D*[P1*Q1 P2*Q2];

end