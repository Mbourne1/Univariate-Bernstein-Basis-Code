function Q1 = BuildQ1(n,t)
% Build the partition of the matrix Q, corresponding to the binomial
% coefficients of v(y). (Note this function works for both v(y) and u(y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Inputs
%
% n : degree of polynomial g(y)
%
% t : index of subresultant S_{t} to be built.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Produce a vector of elements of Q1, then diagonalise it to form a matrix.

% Initialise empty vector Q1.
Q1 = zeros(1,n-t+1);

% for each value of Q, \binom{n-k}{i}
for i = 0:1:n-t
    Q1(i+1) = nchoosek(n-t,i);
end

% Diagonalise vector Q1 to form matrix Q1.
Q1 = diag(Q1);

end