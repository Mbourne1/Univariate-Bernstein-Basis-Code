function DTQ = BuildDTQ_ElementWise(fx,gx,t)
% Build the Matrix D_{k}T_{k}(f,\alpha g)Q
%
%                           Inputs
%
%
% f : Coefficients of polynomial f(\omega,\theta).
%
% g : Coefficients of polynomial g(\omega,\theta).
%
% t : Degree of GCD d(x)


% Global Variables
% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
global SETTINGS

% Get degree of polynomial f
m = GetDegree(fx);

% Get degree of polynomial g
n = GetDegree(gx);

switch SETTINGS.BOOL_LOG
    case 'y' 
        % Use log method
        DT1Q1 = BuildDT1Q1_log(fx,n-t);
        DT2Q2 = BuildDT1Q1_log(gx,m-t);
        
    case 'n' 
        DT1Q1 = BuildDT1Q1_nchoosek(fx,n-t);
        DT2Q2 = BuildDT1Q1_nchoosek(gx,m-t);
    otherwise
        error('SETTINGS.BOOL_LOG must be either y or n')
end


DTQ = [DT1Q1 DT2Q2];

end

function DT1Q1 = BuildDT1Q1_nchoosek(fx,n_t)
% Build D^{-1}T_{n-k}(f)Q_{n-k} partition using matlabs nchoosek function.
%
% Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% n_t : Degree of polynomial v(x,y) = n - t

% Global Variables
global SETTINGS

% Get Degree of input polynomial
m = GetDegree(fx);

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DT1Q1 = zeros(m+n_t+1,n_t+1);

% for each column k in the matrix D^{-1}T_{n-k}(f)Q_{n-k}.
for j = 0:1:n_t
    % for each coefficient a_{i} in the polynomial f(x)
    for i = j:1:m+j
        DT1Q1(i+1,j+1) = ...
            fx(i-j+1) .*...
            nchoosek(m+n_t-i,m-(i-j)) .* ...
            nchoosek(i,j);
    end
end

% Include/Exclude the denominator common to each element of
% D^{-1}T_{n-k}(f)Q_{n-k} in rearranged format.
switch SETTINGS.BOOL_DENOM_SYL
    case 'y' 
        % Common Denominator is included in the coefficient matrix.
        DT1Q1 = DT1Q1 ./ nchoosek(m+n_t,n_t);
    case 'n' 
        % Common Denominator is excluded in the coefficient matrix
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
end
end

function DT1Q1 = BuildDT1Q1_log(f,n_t)
% Build D^{-1}T_{n-k}(f)Q_{n-k} partition using logs function.
%
% Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% n_t : Degree of polynomial v(x,y) = n - t

% Global Variables
global SETTINGS

% Get degree of polynomial f.
m = GetDegree(f);

% Initialise the partition of D^{-1}T_{n-k}(f)Q_{n-k} \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DT1Q1 = zeros(m+n_t+1,n_t+1);

% for each column in partition of D^{-1}T_{n-k}(f)Q_{n-k}
for j = 0:1:n_t
    % for each coefficient f_{i-j} in the polynomial f
    for i = j:1:m+j
        
        % Evaluate binomial coefficients in the numerator in terms of logs
        Numerator_eval_log = ...
            lnnchoosek(m+n_t-i,m-(i-j)) +...
            lnnchoosek(i,j);
        
        % Convert to normal numeric form
        Numerator_eval_exp = 10.^Numerator_eval_log;
        
        % Enter the coefficient in the Sylvester matrix.
        DT1Q1(i+1,j+1) = f(i-j+1) .* Numerator_eval_exp;
        
    end
end

% Include/Exclude the denominator common to each element of
% D^{-1}T_{n-k}(f)Q_{n-k} in rearranged format.
switch SETTINGS.BOOL_DENOM_SYL
    case 'y' % If denominator is included in the coefficient matrix.
        
        % Evaluate the binomial coefficient in the denominator in terms of
        % logs
        Denom_eval_log = lnnchoosek(m+n_t,n_t);
        
        % Convert to normal numeric form
        Denom_eval_exp = 10.^Denom_eval_log;
        
        % Divide the partition by the common denominator.
        DT1Q1 = DT1Q1 ./ Denom_eval_exp ;
    case 'n'
        % Denominator is excluded
    otherwise 
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
end



end