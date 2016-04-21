function [v_max_ai,v_min_ai] = GetMaxMin(fx,n_k)
% Get the maximum and minimum occurence of a_{i} in C_{n-k}(f) for all
% a_{i}, where a_{i} are the m+1 coefficients of the polynomial f(x)

global BOOL_DENOM_SYL

m = GetDegree(fx);
fx = abs(fx);


v_max_ai = zeros(m+1,1);
v_min_ai = zeros(m+1,1);

% For each coefficient
for i = 0:1:m
    v_ai_value = zeros(n_k+1,1);
    for j = 0:1:n_k
 

        cell_val = fx(i+1) * nchoosek(i+j,i) * nchoosek(m+n_k-i-j,m-i) ;
        
        switch BOOL_DENOM_SYL
            case 'y'
                denom = nchoosek(m+n_k,m);
            case 'n'
                denom = 1;
        end
        
        cell_val = cell_val ./ denom;
        
        v_ai_value(j+1) = cell_val;
        
    
    end
    
    v_max_ai(i+1) = max(v_ai_value);
    v_min_ai(i+1) = min(v_ai_value);
    
end



