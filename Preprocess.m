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
global SETTINGS

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Get the mean of the entries.
lambda = GetMean(fx,n-k);

mu = GetMean(gx,m-k);

% Normalize f(x) and g(x) by geometric means
fx_n = fx./ lambda;
gx_n = gx./ mu;


switch SETTINGS.BOOL_ALPHA_THETA
    case 'y'
  
        % For each coefficient ai of F, obtain the max and min such that F_max =
        % [max a0, max a1,...] and similarly for F_min, G_max, G_min
        
        [v_F_max,v_F_min] = GetMaxMin(fx_n,n-k);
        [v_G_max,v_G_min] = GetMaxMin(gx_n,m-k);
        
        
        f_max = max(v_F_max);
        f_min = min(v_F_min);
        g_max = max(v_G_max);
        g_min = min(v_G_min);
        
        
        PrintToFile(f_max,f_min,g_max,g_min,m,n,k,1,1)
        
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
        
        PrintToFile(f_max,f_min,g_max,g_min,m,n,k,alpha,theta)
        
        % Testing to see if preprocessing lowers condition number
        fprintf([mfilename ' : ' sprintf('Condition S(f(x),g(x)) : %2.4e \n', cond(BuildDTQ(fx,gx,k)))]);
        fprintf([mfilename ' : ' sprintf('Conditon S(f(w),alpha g(w)) : %2.4e \n', cond(BuildDTQ(fw,alpha.*gw,k)))]);
        
    case 'n'
        alpha = 1;
        theta = 1;
        
    otherwise
        error('err : Preprocess()');
end

end

function [] = PrintToFile(F_max,F_min,G_max,G_min,m,n,k,alpha,theta)

global SETTINGS

fullFileName = 'Results_Preprocessing.txt';


if exist('Results_Preprocessing.txt', 'file')
    fileID = fopen('Results_Preprocessing.txt','a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s,\t %s,\t %s, \t %s \n',...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(k),...
        SETTINGS.MEAN_METHOD,...
        F_max,...
        F_min,...
        G_max,...
        G_min,...
        num2str(alpha),...
        num2str(theta)...
        );
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end

end