function [lambda] = GetMean(fx,n_k)
% Calculate Geometric means of input polynomials f(x,y) and g(x,y)


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