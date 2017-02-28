function [GM_fx, GM_gx, GM_hx, alpha, theta] = Preprocess_3Polys(fx, gx, hx, k)
% Preprocess_2Polys(fx,gx,k)
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

% Get degree of polynomials f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Get the mean of the entries of f(x) in T_{n-k}(f) and T_{o-k}(f)


GM_fx = GetMean(fx,n-k);
GM_fx2 = GetMean(fx, o-k);
GM_gx = GetMean(gx,m-k);
GM_hx = GetMean(hx,m-k);

% Normalize f(x) and g(x) by geometric means
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;
hx_n = hx./ GM_hx;


if(SETTINGS.BOOL_ALPHA_THETA)
    
    % For each coefficient ai of F, obtain the max and min such that F_max =
    % [max a0, max a1,...] and similarly for F_min, G_max, G_min
    
    % Get maximum and minimum entries of each a_{i} in the first
    % partition T_{n-k}(f)
    [v_F_max1,v_F_min1] = GetMaxMin(fx_n,n-k);
    
    
    % Get maximum and minimum entries of each a_{i} in the second
    % partition T_{o-k}(f)
    [v_F_max2,v_F_min2] = GetMaxMin(fx_n,o-k);
    
    % Get maximum and minimum entries of each b_{i} in the third
    % partition T_{m-k}(g)
    [v_G_max,v_G_min] = GetMaxMin(gx_n,m-k);
    
    % Get maximum and minimum entries of each c_{i} in the third
    % partition T_{m-k}(h)
    [v_H_max,v_H_min] = GetMaxMin(hx_n,m-k);
    
    
    print(v_F_max1,v_F_min1,v_G_max,v_G_min,m,n,k);
    
    % Calculate the optimal value of alpha and theta for the kth
    % subresultant matrix.
    [alpha,theta] = OptimalAlphaTheta(v_F_max1, v_F_min1, v_G_max, v_G_min);
    
    % Having calculated optimal values of alpha and theta, get the max
    % and min entries in C(f) and C(g)
    fx_n = fx ./ GM_fx;
    gx_n = gx ./ GM_gx;
    
    fw = GetWithThetas(fx_n,theta);
    gw = GetWithThetas(gx_n,theta);
    
    % Get max and min values;
    [v_F_max1,v_F_min1] = GetMaxMin(fw,n-k);
    [v_G_max,v_G_min] = GetMaxMin(alpha.*gw,m-k);
    
    print(v_F_max1,v_F_min1,v_G_max,v_G_min,m,n,k);
    
    % Testing to see if preprocessing lowers condition number
    %fprintf([mfilename ' : ' sprintf('Condition S(f(x),g(x)) : %2.4e \n', cond(BuildDTQ(fx,gx,k)))]);
    %fprintf([mfilename ' : ' sprintf('Conditon S(f(w),alpha g(w)) : %2.4e \n', cond(BuildDTQ(fw,alpha.*gw,k)))]);
    
    
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

function [] = print(v_F_max,v_F_min,v_G_max,v_G_min,m,n,k)
% Get maximum entry of all entries of T_{n-k}(f)
f_max = max(v_F_max);

% Get minimum entry of all entries of T_{n-k}(f)
f_min = min(v_F_min);

% Get maximum entry of all entries of T_{m-k}(g)
g_max = max(v_G_max);

% Get minimum entry of all entries of T_{m-k}(g)
g_min = min(v_G_min);

% Print max and minimum entries
PrintToFile(f_max,f_min,g_max,g_min,m,n,k,1,1)
end