function [] = Test_BuildSylvester(ex_num)

addpath 'Examples'

format short

global SETTINGS
SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';

% Get coefficients of f(x,y) and g(x,y)

[f,g] = Examples_GCD_FromCoefficients(ex_num);

% Get degree of polynomial f(x,y)
m = GetDegree(f);

% Get degree of polynomial g(x,y)
n = GetDegree(g);

% Build each subresultant matrix in the SSMS
for k = 1:1:min(m,n)
    
    % Build the matrix D^{-1}_{m+n-k}
    D = BuildD(m,n-k);
    
    % Build the matrix T_{n-k}(f)
    Tf = BuildT1(f,n-k);
    
    % Build the matrix T_{m-k}(g)
    Tg = BuildT1(g,m-k);
    
    % Build the matrix Q = [Q_{n-k} 0 ; 0 Q_{m-k}]
    Q = BuildQ(m,n,k);
    
    % Build the kth Sylvester subresultant matrix
    DTQ{k} = D*[Tf Tg]*Q;
    
    V{k} = BuildV(m,n,k);
    W{k} = BuildW(m,n,k);
    
    if k > 1
        
        DTQ_newMethod{k} = V{k}*DTQ{k-1}*W{k}
                
        display(norm(DTQ_newMethod{k} - DTQ{k}));
                        
        [qv{k},rv{k}] = qr(V{k});    
        [Qs{k-1},Rs{k-1}] = qr(DTQ{k-1});
        [qw{k},rw{k}] = qr(W{k}); 

        [q2,r2] = qr(rv{k} * Qs{k-1})

        
        
        % qw only changes the column ordering, so it is possible to ignore
        % this step. in which case q4 = Identity matrix and r4 = r3
        
        r5 = r2*Rs{k-1}*rw{k}
        q5 = q2 
        
        q5*r5
        
        
        display(W{k});
        display(qw{k});
        display(rw{k});
        
        
    end
    
end


end