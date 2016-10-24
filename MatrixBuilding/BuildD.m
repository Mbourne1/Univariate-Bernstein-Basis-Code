function D = BuildD(m,n_k)
% Build the matrix D^{-1} which is the diagonal matrix of binomial 
% coefficients. Used to construct the Sylvester Subresultant matrices in 
% format D_{k}^{-1}*T_{k}(f,g)*Q.
% 
%
% Input
%
% m : Degree of polynomial f.
%
% n_k : Degree of polynomial v.
%
% n : Degree of polynomial g.
%
% k : Index of subresultant.


% Produce a vector of elements of D, then diagonalise it to form a matrix.
D = diag(1./GetBinomials(m+n_k));

end