function T = BuildT(fx,gx,t)
%  Build the Toeplitz matrix T = [T1 T2], consisting of coefficients of 
% f(x) and g(x).
%
%
% % Input
%
%
% fx: vector of coefficients of f(x) in the standard bernstein basis. a_{i}
%
% gx: vector of coefficients of g(x) in the standard Bernstein basis. b_{i}
%
% t : index of subresultant S_{t} to be formed. (Also degree of GCD)
%
% % Output
%
%
% T : the partitioned matrix T = [T(f) T(g)].
%
%

%%
% Get degree of polynomail f
m = GetDegree(fx);

% Get degree of polynomial g.
n = GetDegree(gx);

% Build Toeplitz matrix of f, the first partiton.
T1 = BuildT1(fx,n-t);

% Build Toeplitz matrix of g, the second partition.
T2 = BuildT1(gx,m-t);

% Concatenate the partitions.
T = [T1 T2];
end
