function [lambda] = GetMean(fx, n_k)
% Calculate Geometric means of the polynomial f(x) in the matrix
% C_{n-k}(f(x))
%
% % Inputs
%
% fx : (Vector) The vector of coefficients of the polynomial f(x)

global SETTINGS


switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        
        lambda = GetGeometricMeanMatlabMethod(fx, n_k);
       
        
    case 'Geometric Mean My Method'
        
        lambda = GetGeometricMean(fx, n_k);
        
        
    case 'Arithmetic Mean'
        
        lambda = GetArithmeticMean(fx, n_k);
        
        
    case 'None'
        
        lambda = 1;
             
    otherwise
        error('err: MEAN_METHOD must be either *Geometric Mean Matlab Method or *Geometric Man My Method or *None')
end
end