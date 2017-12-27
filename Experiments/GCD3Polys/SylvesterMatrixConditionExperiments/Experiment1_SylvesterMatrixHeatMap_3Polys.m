function [] = Experiment1_SylvesterMatrixHeatMap_3Polys(m, n, o, k)
% This experiment considers heat maps of the coefficient multipliers in the
% subresultant matrices
%
% Consider three alternative orderings of polynomials f(x), g(x) and h(x)
% in the k-th sylvester subresutlant matrix
% 
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
% 
% n : (Int) Degree of polynomial g(x)
%
% o : (Int) Degree of polynomial h(x) 
%
% k : (Int) Index of k-th subresultant matrix
%
%
% % Examples
%
% Experiment1_SylvesterMatrixHeatMap_3Polys(28, 18, 19, 1)

close all;
clc;


%
bool_reorder_polynomials = false;

% Get heat map of S(f,g,h)
SylvesterMatrixHeatMap_3Polys(m, n, o, k, bool_reorder_polynomials);

% Get heat map of S(g,f,h)
SylvesterMatrixHeatMap_3Polys(n, m, o, k, bool_reorder_polynomials);

% Get heat map of S(h,f,g)
SylvesterMatrixHeatMap_3Polys(o, m, n, k, bool_reorder_polynomials);

end