function DT1Q1 = BuildDT1Q1(fx,n_k)
% Build Toeplitz matrix for Sylvester Matrix Partitions
%
%
% Inputs.
%
% fx    :  Coefficients of polynomial f given in Bernstein basis n  : degree
%          of polynomial g k  : the index of the subresultant being built
%
% theta :  Optimal value of theta
%
% n     :  Degree of polynomial g(x)
%
% k     :  Index of subresultants S_{k}
%
%

global SYLVESTER_BUILD_METHOD


switch SYLVESTER_BUILD_METHOD
    case 'Standard'
        
        % Get the degree of polynomial f(x)
        m = GetDegree(fx);
        
        % Build matrices D^{-1}
        D = BuildD(m,n_k);
        
        % Build the matrix T1
        T1 = BuildT1(fx,n_k);
        
        % Build the matrix Q1
        Q1 = BuildQ1(n_k);
        DT1Q1 = D*T1*Q1;
        
    case 'Rearranged'
        DT1Q1 = BuildDT1Q1_Rearranged(fx,n_k);
    otherwise
        error('Error : Build Method is either (Standard) or (Rearranged)')
end

end

function DT1Q1 = BuildDT1Q1_Rearranged(fx,n_k)
%
%
%
%
%

global BOOL_LOG
% If Q is included, use the rearrangment such that each Toeplitz
% matrix has a common divisor in each element.
switch BOOL_LOG
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
% Build Toeplitz matrix D^{-1}T(f,g)Q
%
% Inputs.
%
%
% fx :  Coefficients of polynomial f(x) in bernstein basis.
%
% n :   Degree of polynomial g.
%
% k :   Index of subresultant to be built.
%
%

%                       Global Variables.


global BOOL_DENOM_SYL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

switch BOOL_DENOM_SYL
    case 'y' % Include denominator
        
        % Get log of the binomial coefficient.
        Denom_log = lnnchoosek(m+n_k,n_k);
        
        % Get exponent of the log of nchoosek
        Denom_exp = 10.^Denom_log;
        
        % Divide each element of T by the common denominator.
        DT1Q1 = DT1Q1./Denom_exp;
    case 'n'
    otherwise
        error('err')
end
end

function T = BuildDT1Q1_Rearranged_nchoosek(fw,n_k)
% Build Toeplitz matrix D^{-1}T(f,g)Q
%
%
%                           Inputs.
%
% fx :  Coefficients of polynomial f in bernstein basis.
%
% theta :   Optimal value of theta where \theta\omega = y
%
% n :   Degree of polynomial g
%
% k :   Index of the subresultant S_{k}
%
%

%                       Global Variables.
%   BOOL_DENOM - The rearrangment of the Sylvester Matrix in the Bernstein
%   basis reveals common denominators.
%     1 :- Include Common Denominators in the Sylvester Matrix 0 :- Exclude
%     Common Denominators from the Sylvester Matrix.
global BOOL_DENOM_SYL


% Get degree of polynomial f
m = GetDegree(fw);

% Build an empty matrix
T = zeros(m+n_k+1,n_k+1);

% for each column j of the Sylvester matrix DTQ
for j=0:1:n_k
    % for each row from i = j,...
    for i = j:1:(j+m)
        T(i+1,j+1) = ...
            fw(i-j+1) ...
            .* nchoosek(i,j) ...
            .* nchoosek(m+n_k-i,n_k-j);
    end
end


switch BOOL_DENOM_SYL
    case 'y' % Include denominator
        scalar = nchoosek(m+n_k,n_k);
        T = T./scalar;
    case 'n'
    otherwise
        error('error : bool_denom_syl must be either (y) or (n)')
end

end









