function Q = BuildQ(m,n,t)
% Build the diagonal matrix Q corresponding to the binomial coefficients
% of coprime polynomials u and v.
%
%
%
% Inputs
%
% m : Degree of polynomial f.
%
% n : Degree of polynomial g.
%
% t : Index of subresultant S_{t} to be formed.
%
%
%
% Outputs.
%
% Q : The diagonal matrix of binomial coefficients corresponding to coprime
%       polynomials u and v.
%
%

% Build first partition of Q corresponding to the binomial coefficients of
% v(y). \binom{n-k}{i}
Q1 = BuildQ1(n-t);

% Build second partition of Q corresponding to the binomial coefficients of
% u(y). \binom{m-k}{i}
Q2 = BuildQ1(m-t);

% Join the two partitions as a diagonal matrix.
Q = blkdiag(Q1,Q2);
end
