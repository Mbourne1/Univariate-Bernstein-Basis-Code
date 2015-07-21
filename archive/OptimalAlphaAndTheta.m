function [ alpha,theta ] = OptimalAlphaAndTheta(fx_n,gx_n,k,BOOL_Q,BOOL_DENOM,BOOL_LOG)
%OPTIMALALPHAANDTHETA Obtain Optimal Values of alpha and theta
% fx_n - coefficients in Bernstein basis $a_{i}$
% gx_n - coefficients in Bernstein basis $b_{i}$
% k - integer such that the subresultant $S_{k}$ is formed.
% BOOL_Q - whether to include the binomial coefficients of x in the
% Sylvester matrix or not, in the system Ax = b, A is the sylvester matrix,
% x is the solution vector [v;-u]
% BOOL_DENOM - Whether the common denominator in the Sylvester Matrix is
% included.


% For the polynomial fx, calculate the vectors F_max and F_min, such
% that F_max(i) and F_min(i) store the entries of maximum and
% minimum magnitude in s_test that contain the coefficient a(i) of F.
% The vectors G_max and G_min are defined identically, but for G.

% Get the degrees of polynomials f and g
    m = length(fx_n) - 1;
    n = length(gx_n) - 1;

% Build a subresultant such that alpha and theta are 1
    S_test          =   Subresultant(fx_n,gx_n,k,BOOL_Q,BOOL_DENOM,BOOL_LOG);
    
% For each coefficient ai of F, obtain the max and min such that F_max =
% [max a0, max a1,...] and similarly for F_min, G_max, G_min

    [F_max,F_min,G_max,G_min] = sel_maxmin2(S_test,m,n,k);
  

% Calculate the optimal value of alpha and theta for the kth
% subresultant matrix.
    try
        [theta,alpha] = optimal(F_max,F_min,G_max,G_min);
    
               
    
    catch
        fprintf('Failed to obtain optimal values of alpha and theta.\n');
        fprintf('alpha and theta assigned value 1');
        theta = 1;
        alpha = 1; 
    end
end

