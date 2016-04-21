function DT1 = BuildDT1(fw,n_k)

global BOOL_LOG

switch BOOL_LOG
    case 'n'
        % Build Toeplitz Matrix using nchoosek
        DT1 = BuildDT1_nchoosek(fw,n_k);
    case 'y'
        % Build Toeplitz Matrix using log version of nchoosek.
        DT1 = BuildDT1_log(fw,n_k);
    otherwise
        error('error : bool_log must be either (y) or (n)')
end

end


function [DT1] = BuildDT1_nchoosek(fw,n_k)
% Build Toeplitz matrix of D{-1}T(f,g), this is used when we consider
% without Q.
%
%                       Inputs.
%
%
% fx : Input Polynomial f(x)
%
% theta : Optimal value of \theta for change of variable
%
% n : Degree of polynomial g(x)
%
% k : Degree of common divisor d(x)
%
%



% Get Degree of polynomial f
m = GetDegree(fw);

% Buid an empty matrix
DT1 = zeros(m+n-k+1,n-k+1);

% Get fw_bi
fw_bi = fw .* GetBinomials(m);

% for each column j
for j = 0:1:n_k
    %for each row i
    for i = j:1:(j+m)
        DT1(i+1,j+1) = ...
            fw_bi ./ nchoosek(m+n_k,i);
    end
end





end



function T = BuildDT1_log(fw,n_k)
% BuildDT1_log(fw,n_k)
%
% Build Toeplitz matrix of D^{-1}T(f,g) using logs.
%
%
%
%                       Inputs.
%
% fx
%
% theta
%
% n
%
% k
%
%

% Get Degree of polynomial f
m = GetDegree(fx);

% Build an empty matrix
T = zeros(m+n-k+1,n-k+1);

% for each column j
for j = 0:1:n_k
    %for each row i
    for i = j:1:(j+m)
        
        numerator_binom = lnnchoosek(m,i-j);
        denominator_binom = lnnchoosek(m+n-k,i);
        binom = numerator_binom - denominator_binom;
        
        T(i+1,j+1) = fw(i-j+1) .* (10^binom);
    end
end

end
