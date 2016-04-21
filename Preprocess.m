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
global BOOL_ALPHA_THETA

% Get Degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Get the mean of the entries.
lambda = GetMean(fx,n-k);
mu = GetMean(gx,m-k);

% 
fx_n = fx./ lambda;
gx_n = gx./ mu;


switch BOOL_ALPHA_THETA
    case 'y'
  
        % For each coefficient ai of F, obtain the max and min such that F_max =
        % [max a0, max a1,...] and similarly for F_min, G_max, G_min
        
        [v_F_max,v_F_min] = GetMaxMin(fx_n,n-k);
        [v_G_max,v_G_min] = GetMaxMin(gx_n,m-k);
        
        
        f_max = max(v_F_max);
        f_min = min(v_F_min);
        g_max = max(v_G_max);
        g_min = min(v_G_min);
        
        
        PrintToFile(f_max,f_min,g_max,g_min,m,n,k, 'without thetas',1,1)
        
        % Calculate the optimal value of alpha and theta for the kth
        % subresultant matrix.
        [alpha,theta] = OptimalAlphaTheta(v_F_max,v_F_min,v_G_max,v_G_min);
        
        % Having calculated optimal values of alpha and theta, get the max
        % and min entries in C(f) and C(g)
        fx_n = fx ./ lambda;
        gx_n = gx ./ mu;
        
        fw = GetWithThetas(fx_n,theta);
        gw = GetWithThetas(gx_n,theta);
        
        % Get max and min values;
        [v_F_max,v_F_min] = GetMaxMin(fw,n-k);
        [v_G_max,v_G_min] = GetMaxMin(alpha.*gw,m-k);
        
        f_max = max(v_F_max);
        f_min = min(v_F_min);
        g_max = max(v_G_max);
        g_min = min(v_G_min);
        
        PrintToFile(f_max,f_min,g_max,g_min,m,n,k,'With Thetas',alpha,theta)
        
    case 'n'
        alpha = 1;
        theta = 1;
        
    otherwise
        error('err : Preprocess()');
end

end

function [] = PrintToFile(F_max,F_min,G_max,G_min,m,n,k,comment,alpha,theta)

global MEAN_METHOD

fullFileName = 'o_Preprocessing.txt';


if exist('Preprocessing.txt', 'file')
    fileID = fopen('Preprocessing.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %5d \t %d \t %d \t %d \t %s \t %5d \t %5d\n',...
        m,n,k,MEAN_METHOD,F_max,F_min,G_max,G_min,comment,alpha,theta);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end

end