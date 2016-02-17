function [] = DegreeElevateMatrixMethod(fx,r)
% Given the coefficients of the polynomial f(x) in Bernstein basis, and 
% the number of degree elevations r. Perform r fold degree elevation.

% let m be the degree of polynomial f
m = length(fx)-1;

for i = 0:1:m+r
    D(i+1) = nchoosek(m+r,i)
end

D = diag(1./D);

E = zeros(m+r+1,m+1)
% For each column of the matrix E
for j = 0:1:m
    % for each row of the matrix E
    for i = j:1:j+r
        E(i+1,j+1) = nchoosek(r,i-j) * nchoosek(m,j)
    end
end

f = fx';

D*E*f