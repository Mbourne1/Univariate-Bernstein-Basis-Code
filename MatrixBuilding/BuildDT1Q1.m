function DT1Q1 = BuildDT1Q1(fx,n_k)
% Build Toeplitz matrix for Sylvester Matrix Partitions
%
%
% Inputs.
%
% fx  :  Coefficients of polynomial f given in Bernstein basis n  : degree
%          of polynomial g k  : the index of the subresultant being built
%
% n_k :  Degree of polynomial which f is to be multiplied by.

global SETTINGS


switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'Standard'
        
        % Get the degree of polynomial f(x)
        m = GetDegree(fx);
        
        % Build matrices D^{-1}
        D = BuildD(m,n_k);
        
        % Build the matrix T1
        T1 = BuildT1(fx,n_k);
        
        % Build the matrix Q1
        Q1 = BuildQ1(n_k);
        
        % Get the matrix D^{-1} * T_{n-k}(f) * Q_{n-k}
        DT1Q1 = D*T1*Q1;
        
        
        
    case 'Rearranged'
        DT1Q1 = BuildDT1Q1_Rearranged(fx,n_k);
    otherwise
        error('Error : SETTINGS.BUILD_METHOD is either (Standard) or (Rearranged)')
end

end

function DT1Q1 = BuildDT1Q1_Rearranged(fx,n_k)
% Build the matrix D^{-1}*T_{n-k}(f) * Q_{n-k} in its rearranged format.
%
% Inputs 
%
% fx : Coefficients of polynomial f(x,y)
%
% n_k : Degree of polynomial v_{k}(x,y) 

global SETTINGS
% If Q is included, use the rearrangment such that each Toeplitz
% matrix has a common divisor in each element.

switch SETTINGS.BOOL_LOG
    case 'n'
        % Build Toeplitz Matrix using nchoosek
        DT1Q1 = BuildDT1Q1_Rearranged_nchoosek(fx,n_k);
    case 'y'
        % Build Toeplitz Matrix using log version of nchoosek.
        DT1Q1 = BuildDT1Q1_Rearranged_log(fx,n_k);
    otherwise
        error('error : bool_log must be either (y) or (n)')
end
end


function DT1Q1 = BuildDT1Q1_Rearranged_log(fx,n_k)
% Build Toeplitz matrix D^{-1}T_{n-k}(f)Q_{n-k}
%
% Inputs.
%
%
% fx :  Coefficients of polynomial f(x) in Bernstein basis.
%
% n_k : Degree of polynomial v(x).


% Global Variables
global SETTINGS

% Get degree of polynomial f(w)
m = GetDegree(fx);

% Build an empty matrix
DT1Q1 = zeros(m+n_k+1,n_k+1);


% for each column
for j=0:1:n_k
    % for each entry in the vector f.
    for i = j:1:j+m
        
        % Get the two binomial coefficients in the numerator in terms of
        % logs.
        Numerator_log = lnnchoosek(i,j) + lnnchoosek(m+n_k-i,n_k-j);
        
        % Convert to normal form.
        Numerator_exp = 10.^Numerator_log;
        
        % Multiply binomial evaluation by coefficient f{i-j} from
        % polynomial f.
        DT1Q1(i+1,j+1) =  fx(i-j+1).* Numerator_exp  ;
        
    end
end

% Include/Exclude the denominator which is common to all elements of
% D^{-1}*T_{n-k}(f)*Q_{n-k}.

switch SETTINGS.BOOL_DENOM_SYL
    case 'y' % Include denominator
        
        % Get log of the binomial coefficient.
        Denom_log = lnnchoosek(m+n_k,n_k);
        
        % Get exponent of the log of nchoosek
        Denom_exp = 10.^Denom_log;
        
        % Divide each element of T by the common denominator.
        DT1Q1 = DT1Q1./Denom_exp;
    case 'n' % Exclude denominator
        
    otherwise
        error('err')
end
end

function T = BuildDT1Q1_Rearranged_nchoosek(fx,n_k)
% Build the matrix D^{-1}*T_{n-k}(f)Q_{n-k}
%
% Inputs.
%
% fx :  Coefficients of polynomial f in bernstein basis.
%
% n_k :   Degree of polynomial v(x)
%


% Global variables
global SETTINGS


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Build an empty matrix
T = zeros(m+n_k+1,n_k+1);

% for each column j of the Sylvester matrix DTQ
for j=0:1:n_k
    % for each row from i = j,...
    for i = j:1:(j+m)
        T(i+1,j+1) = ...
            fx(i-j+1) ...
            .* nchoosek(i,j) ...
            .* nchoosek(m+n_k-i,n_k-j);
    end
end


switch SETTINGS.BOOL_DENOM_SYL
    case 'y' % Include denominator
        scalar = nchoosek(m+n_k,n_k);
        T = T./scalar;
    case 'n'
    otherwise
        error('error : bool_denom_syl must be either (y) or (n)')
end

end









