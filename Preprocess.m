function [lambda, mu, alpha, theta] = Preprocess(fx,gx,k)
% Preprocess(fx,gx,k)
% Get the optimal values of lamdba, mu, alpha and theta.
%
% Inputs.
%
%
% fx : Vector of coefficients of f(x)
%
% gx : Vector of coefficients of g(x)
%
% k : Degree of common divisor d_{k}(x)
%
% Outputs.
%
% lamda : Geometric mean of entries of f(x) in k-th subresultant matrix
%
% mu : Geometric mean of entries of g(x) in k-th subresultant matrix
%
% alpha : Optimal value of alpha
%
% theta : Optimal value of theta


% Global variables
global BOOL_PREPROC

% Get Degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

switch BOOL_PREPROC
    case 'y'
        
        % Get Unprocessed partitions (Including Geometric Mean)
        C_f_unproc = BuildT1(fx,n-k);
        C_g_unproc = BuildT1(gx,m-k);
        
        % Reason for performing BuildToeplitz before taking geometric
        % mean is that the MATLAB function geomean requires the
        % unprocessed matrix to calculate the GM.
        
        [lambda, mu] = GetGeometricMean(fx,gx,k);
        
        % Divide the Sylvester Matrix partitions by Geometric mean.
        C_f_unproc_gm = C_f_unproc ./ lambda;
        C_g_unproc_gm = C_g_unproc ./ mu;
        
        % Build subresultant S_{k}, and add to array of Sk
        Sylvester_unproc = [C_f_unproc_gm C_g_unproc_gm];
        
        % For each coefficient ai of F, obtain the max and min such that F_max =
        % [max a0, max a1,...] and similarly for F_min, G_max, G_min
        [F_max,F_min,G_max,G_min] = GetMaxMin(Sylvester_unproc,m,n,k);
        
        % Calculate the optimal value of alpha and theta for the kth
        % subresultant matrix.
        [alpha,theta] = OptimalAlphaTheta(F_max,F_min,G_max,G_min);
        
    case 'n'
        alpha = 1;
        theta = 1;
        mu = 1;
        lambda = 1;
        
    otherwise
        error('err : Preprocess()');
end

end