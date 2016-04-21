function [lambda] = GetMean(fx,n_k)
% Calculate Geometric means of input polynomials f(x,y) and g(x,y)


global MEAN_METHOD


switch MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        C_f_unproc = BuildDT1Q1(fx,n_k);
        lambda = geomean(abs(C_f_unproc(C_f_unproc~=0)));
        
    case 'Geometric Mean My Method'
        lambda = GetGeometricMean(fx,n_k);
        
    case 'None'
        lambda = 1;
             
    otherwise
        error('err: MEAN_METHOD must be either matlab or mymethod')
end
end