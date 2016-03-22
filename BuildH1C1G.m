function H1C1G = BuildH1C1G(uw,t)
% Build the matrix H1C1G
%
% Inputs.
%
% uw : input polynomial
%
% t : degree of GCD.


% Global Variables
global BOOL_LOG

switch BOOL_LOG
    case 'y' 
        % Use logs
        H1C1G = BuildH1C1G_log(uw,t);
    case 'n' 
        % Use nchoosek
        H1C1G = BuildH1C1G_nchoosek(uw,t);
    otherwise
        error('err')
end

end

function H1C1G = BuildH1C1G_nchoosek(uw,t)
% Build Partition of the HCG matrix using nchoosek
%
%
% Inputs.
%
%
% uw :  Coefficients of polynomial u in scaled bernstein basis, to be put
%       into matrix form.
%
% t  :- Degree of GCD
%
%


global BOOL_DENOM_APF
%

% Get degree of polynomial u(w)
m_minus_t = size(uw,1) - 1;

% Get degree of polynomial f(w)
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

switch BOOL_DENOM_APF
    case 'y' 
        % Include Common Denominator in the matrix
        H1C1G  = H1C1G ./ nchoosek(m,t);
    case 'n'
        % do nothing
    otherwise 
        error('err')
        
end
end

function H1C1G = BuildH1C1G_log(uw,t)
% Build the partition H1C1G where HCG = [H1C1G | H2C2G]
%
%                           Inputs.
%
%
% uw :  Input polynomial
%
% t :   Degree of GCD
%

%%
%                           Global Variables.

global BOOL_DENOM_APF

%%

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
switch BOOL_DENOM_APF
    case 'y' % Include the common denominator
        
        Denom_Eval_log = lnnchoosek(m,t);
        
        Denom_Eval_exp  = 10.^Denom_Eval_log;
        
        H1C1G = H1C1G ./Denom_Eval_exp;
end

end