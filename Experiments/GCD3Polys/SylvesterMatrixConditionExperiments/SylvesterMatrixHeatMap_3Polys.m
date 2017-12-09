function [] = SylvesterMatrixHeatMap_3Polys(m, n, o, k)
%
% % Inputs
%
% m : (Int) Degree of f(x)
%
% n : (Int) Degree of g(x)
%
% o : (Int) Degree of h(x)
%
% k : (Int) index of Sylvester subresultant matrix
%
% % Outputs
%
%
% % Example
%
% >> SylvesterMatrixHeatMap_3Polys(5, 15, 7, 3, true)


fx = ones(m + 1, 1);
gx = ones(n + 1, 1);
hx = ones(o + 1, 1);


%if bool_reorder_polys
%    [fx, gx, hx, m, n, o] = ReorderPolys(fx, gx, hx, m, n, o);
%end


arrSubresultantFomat = {'T', 'DT', 'TQ', 'DTQ'};
SYLVESTER_EQUATIONS = '2';


nFormats = length(arrSubresultantFomat);


for i = 1:1:nFormats
    
    subresultantFomat = arrSubresultantFomat{i};
    
    
    switch subresultantFomat
        
        case 'T'
            
            switch SYLVESTER_EQUATIONS
                
                case '2'
                    T = BuildT_3Polys(fx, gx, hx, k);
                    
                case '3'
                    T = BuildT_3Polys_New(fx, gx, hx, k);
                    
                otherwise
                    error('err')
                    
            end
            
            
            Sk = T;
            
        case 'DT'
            
            switch SYLVESTER_EQUATIONS
                
                case '2'
                    
                    D = BuildD_3Polys(m, n - k, o - k);
                    T = BuildT_3Polys(fx, gx, hx, k);
                    
                case '3'
                    D = BuildD_3Polys_New(m, n, o, k);
                    T = BuildT_3Polys_New(fx, gx, hx, k);
                otherwise
                    error('err')
            end
            
            Sk = D*T;
            
        case 'TQ'
            
            switch     SYLVESTER_EQUATIONS
                case '2'
                    T = BuildT_3Polys(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m, n, o, k);
                case '3'
                    T = BuildT_3Polys_New(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m,n,o,k);
                otherwise
                    error('err')
                    
                    
            end
            
            Sk = T*Q;
            
            
        case 'DTQ'
            
            switch SYLVESTER_EQUATIONS
                case '2'
                    D = BuildD_3Polys(m, n-k, o-k);
                    T = BuildT_3Polys(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m, n, o, k);
                    
                case '3'
                    D = BuildD_3Polys_New(m,n,o,k);
                    T = BuildT_3Polys_New(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m,n,o,k);
                otherwise
                    error('err')
            end
            
            Sk = D*T*Q;
            
            
            
        case 'DTQ Denominator Removed'
            
            error('not completed')
            
            
        otherwise
            error('%s is not a valid format', subresultantFomat)
            
    end
    
    figure_name = sprintf('%s : Heat Map', subresultantFomat);
    figure('Name',figure_name)
    
    % Note - neccesary to take absolute values since one partition is
    % negative
    Sk_rounded = log10(abs(Sk));
    
    hm = heatmap((Sk_rounded));
    
end


end


function [fx,gx,hx,m,n,o] = ReorderPolys(fx,gx,hx,m,n,o)

[~, index] = min([m n o]);


switch index
    
    case 1 % fx is minimum so dont change order
        
        
    case 2 % gx is iminimum, swap f and g
        
        [fx,gx] = PolySwap(fx,gx);
        [m, n] = DegreeSwap(m,n);
        
    case 3 % hx is minimum, swap f and h
        
        [fx, hx] = PolySwap(fx, hx);
        [m, o] = DegreeSwap(m, o);
        
    otherwise
        
        
end


end


function [fx, gx] = PolySwap(fx, gx)


fx_temp = fx;
fx = gx;
gx = fx_temp;



end


function [m,n] = DegreeSwap(m,n)


m_temp = m;

m = n;
n = m_temp;

end
