function [lambda, mu] = GetGeometricMean(fx,gx,k)
% Calculate Geometric means of input polynomials f(x,y) and g(x,y)

global GEOMETRIC_MEAN_METHOD


% Get degree of polynomial f(x)
[r,~] = size(fx);
m = r-1;

% Get degree of polynomial g(x)
[r,~] = size(gx);
n = r-1;


switch GEOMETRIC_MEAN_METHOD
    case 'matlab'
        
        C_f_unproc = BuildT1(fx,1,n,k);
        C_g_unproc = BuildT1(gx,1,m,k);
        
        % Use Matlab method for calculating GM
        lambda = geomean(abs(C_f_unproc(C_f_unproc~=0)));
        mu = geomean(abs(C_g_unproc(C_g_unproc~=0)));
        
    case 'mymethod'
        % Get Geometric means my method
        lambda = GeometricMean(fx,n,k);
        mu = GeometricMean(gx,m,k);
        
    case 'none'
        % Do not calculate GeometricMean
        lambda(k) = 1;
        mu(k) = 1;
        
    case 'mean'
        lambda(k) = mean(fx);
        mu(k) = mean(gx);
        
    otherwise
        error('err: geometricMeanMethod must be either matlab or mymethod')
end
end