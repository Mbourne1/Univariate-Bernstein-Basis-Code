function [vMax_ai,vMin_ai] = GetMaxMin(fx,n_k)
% Get the maximum and minimum occurence of a_{i} in C_{n-k}(f) for all
% a_{i}, where a_{i} are the m+1 coefficients of the polynomial f(x)

global SETTINGS

% Get the degree of f(x)
m = GetDegree(fx);

% Get absolute values of f(x).
fx = abs(fx);

% Initialise vectors to store maximum and minimum of each a_{i} in
% T_{n-k}(f)
vMax_ai = zeros(m+1,1);
vMin_ai = zeros(m+1,1);


% First determine the denominator of each entry
switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        denom = nchoosek(m+n_k,m);
    case 'n'
        denom = 1;
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL is either y or n')
end

% For each coefficient a_{i} of f(x)
for i = 0:1:m
    
    % Initialise a vector to store a_{i} from each column
    v_ai_value = zeros(n_k+1,1);

    % For each column 0,...,n-k of T_{n-k}(f), get entry which contains
    % a_{i,j}
    for j = 0:1:n_k
        
        
        switch SETTINGS.SYLVESTER_BUILD_METHOD
            case 'Standard'
                
                
                cell_val = fx(i+1) * nchoosek(m,i) * nchoosek(n_k,j) ./ nchoosek(m+n_k,i+j);
            case 'Rearranged'
                
                cell_val = fx(i+1) * nchoosek(i+j,i) * nchoosek(m+n_k-i-j,m-i) ./ denom ;
                
            otherwise 
                error([mfilename ' : ''err']);
        end
        
        
        % Store entry in vector.
        v_ai_value(j+1) = cell_val;
        
        
    end
    
    % Get maximum entry of a_{i} in T_{n-k}(f), and store in vector.
    vMax_ai(i+1) = max(v_ai_value);
    
    % Get minimum entry of a_{i} in T_{n-k}(f), and store in vector.
    vMin_ai(i+1) = min(v_ai_value);
    
end



