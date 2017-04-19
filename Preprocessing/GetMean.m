function [lambda] = GetMean(fx,n_k)
% Calculate Geometric means of input polynomials f(x,y) and g(x,y)


global SETTINGS


switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        % Build the partition of the Sylvester matrix
        C_f_unproc = BuildSubresultant_Partition_2Polys(fx,n_k);
        
        % Get geometric mean of non-zero entries
        lambda = geomean(abs(C_f_unproc(C_f_unproc~=0)));
        
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