function DTQ = BuildDTQ_ElementWise(fw,gw,alpha,t)
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


%

%                       Global Variables

% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
global BOOL_LOG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomial f
m = size(fw,1) - 1;

% Get degree of polynomial g
n = size(gw,1) - 1;

switch BOOL_LOG
    case 'y' 
        % Use log method
        DT1Q1 = BuildDT1Q1_log(fw,n-t);
        DT2Q2 = BuildDT1Q1_log(gw,m-t);
        
    case 'n' 
        % Use nchoosek method
        warning('off','all');
        DT1Q1 = BuildDT1Q1_nchoosek(fw,n-t);
        DT2Q2 = BuildDT1Q1_nchoosek(gw,m-t);
        warning('on','all');
    otherwise
        error('bool_log must be either y or n')
end


DTQ = [DT1Q1 alpha.*DT2Q2];

end

function DT1Q1 = BuildDT1Q1_nchoosek(fw,n_t)
% Build DTQ partition by, using matlabs nchoosek function.
%
%                       Inputs
%
%
% f : Coefficients of polynomial f(\omega,\theta)
%
% n : Degree of other polynomial g(\omega,\theta).
%
% t : Degree of GCD.

%%
%                   Global Variables

% BOOL_DENOM_SYL - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global BOOL_DENOM_SYL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get Degree of input polynomial
m = size(fw,1) - 1;

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DT1Q1 = zeros(m+n_t+1,n_t+1);

%fw = fx .* theta.^(0:1:m);

% for each column k in the partition of DTQ.
for j = 0:1:n_t
    % for each coefficient in the polynomial f
    for i = j:1:m+j
        DT1Q1(i+1,j+1) = ...
            fw(i-j+1) .*...
            nchoosek(m+n_t-i,m-(i-j)) .* ...
            nchoosek(i,j);
    end
end

switch BOOL_DENOM_SYL
    case 'y' 
        % Common Denominator is included in the coefficient matrix.
        DT1Q1 = DT1Q1 ./ nchoosek(m+n_t,n_t);
    case 'n' 
        % Common Denominator is excluded in the coefficient matrix
    otherwise
        error('bool_denom_syl must be either y or n')
end
end

function DT1Q1 = BuildDT1Q1_log(f,theta, n,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                        Inputs

% f : - Coefficients of polynomial f(\omega,\theta).

% n :- degree of other polynomial.

% t :- Degree of GCD.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                      Global Variables

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global BOOL_DENOM_SYL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomial f.
m = length(f)-1;

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DT1Q1 = zeros(m+n-t+1,n-t+1);

% for each column in partition of DTQ
for j = 0:1:n-t
    % for each coefficient f_{i-j} in the polynomial f
    for i = j:1:m+j
        
        % Evaluate binomial coefficients in the numerator in terms of logs
        Numerator_eval_log = ...
            lnnchoosek(m+n-t-i,m-(i-j)) +...
            lnnchoosek(i,j);
        
        % Convert to normal numeric form
        Numerator_eval_exp = 10.^Numerator_eval_log;
        
        % Enter the coefficient in the Sylvester matrix.
        DT1Q1(i+1,j+1) = fw(i-j+1) .* Numerator_eval_exp;
        
    end
end

switch BOOL_DENOM_SYL
    case 'y' % If denominator is included in the coefficient matrix.
        
        % Evaluate the binomial coefficient in the denominator in terms of
        % logs
        Denom_eval_log = lnnchoosek(m+n-t,n-t);
        
        % Convert to normal numeric form
        Denom_eval_exp = 10.^Denom_eval_log;
        
        % Divide the partition by the common denominator.
        DT1Q1 = DT1Q1 ./ Denom_eval_exp ;
    case 'n'
        % Denominator is excluded
    otherwise 
        error('err')
end



end