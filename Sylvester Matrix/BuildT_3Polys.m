function T = BuildT_3Polys(fx, gx, hx, k)
%  Build the Toeplitz matrix T = [T1 T2], consisting of coefficients of 
% f(x) and g(x).
%
%
% % Input
%
% [fx, gx, hx] : Vector of coefficients of f(x), g(x) and h(x) in the 
% Bernstein basis. a_{i}
%
% k : index of subresultant S_{k} to be formed. (Also degree of GCD)
%
% % Output
%
% T : the partitioned matrix T = [T(f) T(g)].
%

% Get the degree of polynomail f(x)
m = GetDegree(fx);

% Get the degree of polynomial g(x)
n = GetDegree(gx);

% Get the degree of polynomial h(x)
o = GetDegree(hx);

% % Build the block diagonal matrix of T_{n-k}(f) and T_{o-k}(f)
% Build T1
T1 = BuildT1(fx,n-k);

% Build T2
T2 = BuildT1(fx,o-k);

diag_sec = blkdiag(T1,T2);

% % Build the column section of the matrix 

% Build T_{m-k}(g)
T3 = BuildT1(gx,m-k);

% Build T_{m-k}(h)
T4 = BuildT1(hx,m-k);

% Build
col_sec = [T3 ; T4];

% Concatenate the partitions.
T = [diag_sec col_sec];


end
