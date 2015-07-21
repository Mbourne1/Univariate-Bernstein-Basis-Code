
function H1C1G = BuildHCGPart(uw,t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.

% uw - input polynomial

% t - degree of GCD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables


% bool_log - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.
global bool_log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch bool_log
    case 1 % Use logs
        H1C1G = BuildHCGPart_log(uw,t);
    case 0 % Use nchoosek
        H1C1G = BuildHCGPart_nchoosek(uw,t);
end

end

function H1C1G = BuildHCGPart_nchoosek(uw,t)
% Build Partition of the HCG matrix using nchoosek

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs.

% uw :  Coefficients of polynomial u in scaled bernstein basis, to be put
%       into matrix form.

% t  :- Degree of GCD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global bool_denom_apf

m_minus_t = length(uw)-1;

m = m_minus_t + t;

H1C1G = zeros(m+1,t+1);

% for each column 0:1:t
for j = 0:1:t
    %for each row
    for i = j:1:(m_minus_t)+j
        H1C1G(i+1,j+1) = ...
            uw(i-j+1)...
            .* nchoosek(i,j) ...
            .* nchoosek(m-i,t-j);
        
    end
end

switch bool_denom_apf
    case 1 % Common Denominator included
        H1C1G  = H1C1G ./ nchoosek(m,t);
end
end

function H1C1G = BuildHCGPart_log(uw,t)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% uw :  Input polynomial

% t :   Degree of GCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables.

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
% (NOTE THIS IS ALWAYS ONE)
global bool_denom_apf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree of polynomial uw, deg(u) = m-t.
m_minus_t = length(uw)-1;

% Get m - the degree of polynomial f.
m = m_minus_t + t;

H1C1G = zeros(m,t);
% for each column 0:1:t
for j = 0:1:t
    %for each row
    for i = j:1:(m_minus_t)+j
        
        Numerator_eval_log = lnnchoosek(i,j) + lnnchoosek(m-i,t-j);
        
        Num_eval_exp = 10.^Numerator_eval_log;
        
        H1C1G(i+1,j+1) = uw(i-j+1) .* Num_eval_exp;
        
        
    end
end

% If include common denominator of the partition of HCG
switch bool_denom_apf
    case 1 % Include the common denominator
        
        Denom_Eval_log = lnnchoosek(m,t);
        
        Denom_Eval_exp  = 10.^Denom_Eval_log;
        
        H1C1G = H1C1G ./Denom_Eval_exp;
end

end