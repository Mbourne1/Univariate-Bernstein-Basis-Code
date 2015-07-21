function DC1Q1 = BuildToeplitz(fx,theta,n,k)
% Build Toeplitz matrix for Sylvester Matrix Partitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs


% fx :  Coefficients of polynomial f given in Bernstein basis n  : degree
%       of polynomial g k  : the index of the subresultant being built

% theta :   Optimal value of theta

% n :   Degree of polynomial g(x)

% k :   Index of subresultants S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.

% BOOL_Q (Boolean)
%   1 :- Q included in the Sylvester Matrix S(f,g) = D^{-1}T(f,g)Q 0 :- Q
%   excluded from Sylvester Matrix S(f,g) = D^{-1}T(f,g)
global bool_q

% BOOL_LOG -
%   1 :- Use logs to calculate coefficients 0 :- Use nchoosek to calculate
%   coefficients.
global bool_log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch bool_q
    case 0
        % If Q is not included, use the standard Sylvester form
        % D^{-1}S(f,g).
        switch bool_log
            case 0
                % Build Toeplitz Matrix using nchoosek
                DC1Q1 = BuildDT_nchoosek(fx,theta,n,k);
            case 1
                % Build Toeplitz Matrix using log version of nchoosek.
                DC1Q1 = BuildDT_log(fx,theta,n,k);
                
        end
    case 1
        % If Q is included, use the rearrangment such that each Toeplitz
        % matrix has a common divisor in each element.
        switch bool_log
            case 0
                % Build Toeplitz Matrix using nchoosek
                DC1Q1 = BuildDTQ_Rearranged_nchoosek(fx,theta,n,k);
            case 1
                % Build Toeplitz Matrix using log version of nchoosek.
                DC1Q1 = BuildDTQ_Rearranged_log(fx,theta,n,k);
        end
        
end


end


function DC1Q1 = BuildToeplitz_Naive(fx,theta,n,k)

m = length(fx)-1;


% Build matrix D
D = BuildD(m,n,k);

% Build matrix T
% for each column j
for j = 0:1:n-k
    % for each row i
    for i = j:1:m+j
        T1(i+1,j+1) = fx(i-j+1) .* theta^(i-j) .* nchoosek(m,i-j);
    end
end
        
% Build matrix Q1
Q1 = BuildQ1(n,t);

end




function DC1Q1 = BuildDTQ_Rearranged_log(fx,theta,n,k)
% Build Toeplitz matrix D^{-1}T(f,g)Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs. 


% fx :  Coefficients of polynomial f in bernstein basis. 

% n :   Degree of polynomial g. 

% k :   Index of subresultant to be built. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.


global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree of polynomial f
m = length(fx)-1;

% Build an empty matrix
DC1Q1 = zeros(m+n-k+1,n-k+1);


% for each column
for j=0:1:n-k
    % for each entry in the vector f.
    for i = j:1:j+m
        
        % Get the two binomial coefficients in the numerator in terms of
        % logs.
        Numerator_log = lnnchoosek(i,j) + lnnchoosek(m+n-k-i,n-k-j);
        
        % Convert to normal form.
        Numerator_exp = 10.^Numerator_log;
        
        % Multiply binomial evaluation by coefficient f{i-j} from
        % polynomial f.
        DC1Q1(i+1,j+1) =  fx(i-j+1) .*theta^(i-j) .* Numerator_exp  ;
        
    end
end

switch bool_denom_syl
    case 1 % Include denominator
        
        % Get log of the binomial coefficient.
        Denom_log = lnnchoosek(m+n-k,n-k);
        
        % Get exponent of the log of nchoosek
        Denom_exp = 10.^Denom_log;
        
        % Divide each element of T by the common denominator.
        DC1Q1 = DC1Q1./Denom_exp;
end
end

function T = BuildDTQ_Rearranged_nchoosek(fx,theta,n,k)
% Build Toeplitz matrix D^{-1}T(f,g)Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.

% fx :  Coefficients of polynomial f in bernstein basis. 

% theta :   Optimal value of theta where \theta\omega = y

% n :   Degree of polynomial g

% k :   Index of the subresultant S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.

%   BOOL_DENOM - The rearrangment of the Sylvester Matrix in the Bernstein
%   basis reveals common denominators.
%     1 :- Include Common Denominators in the Sylvester Matrix 0 :- Exclude
%     Common Denominators from the Sylvester Matrix.
global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomial f
m = length(fx)-1;

% Build an empty matrix
T = zeros(m+n-k+1,n-k+1);

% for each column j of the Sylvester matrix DTQ
for j=0:1:n-k
    % for each row from i = j,...
    for i = j:1:(j+m)
        T(i+1,j+1) = fx(i-j+1) .* theta^(i-j) * nchoosek(i,j) * nchoosek(m+n-k-i,n-k-j);
    end
end
switch bool_denom_syl
    case 1 % Include denominator
        scalar = nchoosek(m+n-k,n-k);
        T = T./scalar;
end

end


function T = BuildDT_nchoosek(fx,theta,n,k)
% Build Toeplitz matrix of D{-1}T(f,g),
% this is used when we consider without Q.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs.

% fx

% theta

% n

% k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Get Degree of polynomial f
m = length(fx)-1;
% Buid an empty matrix
T = zeros(m+n-k+1,n-k+1);

% for each column j
for j = 0:1:n-k
    %for each row i
    for i = j:1:(j+m)
        T(i+1,j+1) = fx(i-j+1) .*theta^(i-j) .* nchoosek(m,i-j) ./ nchoosek(m+n-k,i);
    end
end

end

function T = BuildDT_log(fx,theta,n,k)
% Build Toeplitz matrix of D^{-1}T(f,g) using logs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs.

% fx

% theta

% n

% k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Degree of polynomial f
m = length(fx) - 1;

% Build an empty matrix
T = zeros(m+n-k+1,n-k+1);

% for each column j
for j = 0:1:n-k
    %for each row i
    for i = j:1:(j+m)
        
        numerator_binom = lnnchoosek(m,i-j);
        denominator_binom = lnnchoosek(m+n-k,i);
        binom = numerator_binom - denominator_binom;
        
        T(i+1,j+1) = fx(i-j+1) .*theta^(i-j) * (10^binom);
    end
end

end


function nCk = lnnchoosek(n,k)
% Perform nchoosek in logs.

nCk = log10(factorial(n)) - log10(factorial(n-k)) - log10(factorial(k));
end






