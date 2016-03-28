function [m1,m2] = GetDegree_Bivariate(fxy)
% GetDegree(fxy)
%
% Get the degree of the polynomial f(x,y), whose coefficients are given as 
% a matrix.

% Get number of rows and number of columns in matrix of coefficients of f(x,y)
[r,c] = size(fxy);

% Get degree m1 of f(x,y) with respect to x
m1 = r - 1;

% Get degree m2 of f(x,y) with respect to y
m2 = c - 1;

end