function [] = Test_BuildSylvester(ex_num)
% Testing the new Method of building Sylvester matrices given a previous
% Sylvester subresultant Matrix has been constructed


addpath 'Examples'

format short

global SETTINGS
SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';

% Get coefficients of f(x,y) and g(x,y)

[fx, gx] = Examples_GCD_FromCoefficients(ex_num);

% Get degree of polynomial f(x,y)
m = GetDegree(fx);

% Get degree of polynomial g(x,y)
n = GetDegree(gx);


arrDTQ = cell(min(m,n)+1,1);
arrDTQ_newMethod = cell(min(m,n) + 1,1);

% Build each subresultant matrix in the SSMS
for k = 1:1:min(m, n)
    
    % Build the matrix D^{-1}_{m+n-k}
    D = BuildD(m, n - k);
    
    % Build the matrix T_{n-k}(f)
    Tf = BuildT1(fx, n - k);
    
    % Build the matrix T_{m-k}(g)
    Tg = BuildT1(gx, m - k);
    
    % Build the matrix Q = [Q_{n-k} 0 ; 0 Q_{m-k}]
    Q = BuildQ(m, n, k);
    
    % Build the kth Sylvester subresultant matrix
    arrDTQ{k} = D*[Tf Tg]*Q;
    
    arrV{k} = BuildV(m, n, k);
    arrW{k} = BuildW(m, n, k);
    
    if k > 1
        
        arrDTQ_newMethod{k} = arrV{k}*arrDTQ{k-1}*arrW{k}
                
        display(norm(arrDTQ_newMethod{k} - arrDTQ{k}));
                        
        [qv{k},rv{k}] = qr(arrV{k});    
        [Qs{k-1},Rs{k-1}] = qr(arrDTQ{k-1});
        [qw{k},rw{k}] = qr(arrW{k}); 

        [q2,r2] = qr(rv{k} * Qs{k-1})

        
        
        % qw only changes the column ordering, so it is possible to ignore
        % this step. in which case q4 = Identity matrix and r4 = r3
        
        r5 = r2*Rs{k-1}*rw{k}
        q5 = q2 
        
        q5*r5
        
        
        display(arrW{k});
        display(qw{k});
        display(rw{k});
        
        
    end
    
end


end