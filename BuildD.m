function D = BuildD(m,n,k)
% Build Matrix D^{-1} the diagonal matrix of binomial coefficients. Used to
% construct the Sylvester Subresultant matrices D_{k}^{-1}T_{k}(f,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                Input
%
% m : Degree of polynomial f.
%
% n : Degree of polynomial g.
%
% k : Index of subresultant.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Produce a vector of elements of D, then diagonalise it to form a matrix.

% Initialise empty vector D.
D = zeros(1,m+n-k+1);

% for each element in D, assign value \frac{1}{\binom{m+n-k}{i}}
for i = 0:1:m+n-k
    D(i+1) = nchoosek(m+n-k,i);
end

% Diagonalise the vector to form square matrix D.
D = diag(1./D);

end