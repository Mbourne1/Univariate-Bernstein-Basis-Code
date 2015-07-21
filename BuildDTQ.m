
function DTQ = BuildDTQ(fx,gx,alpha, theta, t)
% Build the Matrix D_{k}T_{k}(f,\alpha g)Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% f : - Coefficients of polynomial f(\omega,\theta).

% alpha_g :- Coefficients of \alpha.*g(\omega,\theta).

% t:- Degree of GCD.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables

% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
global bool_log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomial f
m = length(fx) -1;

% Get degree of polynomial g
n = length(gx) -1;

switch bool_log
    case 1 % Use log method
        DTQ1 = BuildDTQ_Partition_log(fx,theta,n,t);
        DTQ2 = BuildDTQ_Partition_log(gx,theta,m,t);
    case 0 % Use nchoosek method
        warning('off','all');
        DTQ1 = BuildDTQ_Partition_nchoosek(fx,theta,n,t);
        DTQ2 = BuildDTQ_Partition_nchoosek(gx,theta,m,t);
        warning('on','all');
end


DTQ = [DTQ1 alpha.*DTQ2];

end

function DTQ1 = BuildDTQ_Partition_nchoosek(f,theta,n,t)
% Build DTQ partition by, using matlabs nchoosek function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs

% f : - Coefficients of polynomial f(\omega,\theta).

% n :- degree of other polynomial.

% t :- Degree of GCD.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                   Global Variables

% bool_denom_syl - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get Degree of input polynomial
m = length(f)-1;

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DTQ1 = zeros(m+n-t+1,n-t+1);

% for each column k in the partition of DTQ.
for j = 0:1:n-t
    % for each coefficient in the polynomial f
    for i = j:1:m+j
        DTQ1(i+1,j+1) = ...
            f(i-j+1) .* theta^(i-j) .*...
            nchoosek(m+n-t-i,m-(i-j)) .* ...
            nchoosek(i,j);
    end
end


% for each column k in the partition of DTQ.
for j = 0:1:n-t
    % for each coefficient in the polynomial f
    for i = j:1:m+j
        DTQ1(i+1,j+1) = ...
            f(i-j+1) .* theta^(i-j) .*...
            nchoosek(m+n-t-i,m-(i-j)) .* ...
            nchoosek(i,j);
    end
end

switch bool_denom_syl
    case 1 % Common Denominator is included in the coefficient matrix.
        DTQ1 = DTQ1 ./ nchoosek(m+n-t,n-t);
end
end

function DTQ1 = BuildDTQ_Partition_log(f,theta, n,t)
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
global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomial f.
m = length(f)-1;

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DTQ1 = zeros(m+n-t+1,n-t+1);

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
        DTQ1(i+1,j+1) = f(i-j+1) .* theta^(i-j) .* Numerator_eval_exp;
        
    end
end

switch bool_denom_syl
    case 1 % If denominator is included in the coefficient matrix.
        
        % Evaluate the binomial coefficient in the denominator in terms of
        % logs
        Denom_eval_log = lnnchoosek(m+n-t,n-t);
        
        % Convert to normal numeric form
        Denom_eval_exp = 10.^Denom_eval_log;
        
        % Divide the partition by the common denominator.
        DTQ1 = DTQ1 ./ Denom_eval_exp ;
        
end



end