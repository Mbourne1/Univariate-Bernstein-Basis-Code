function T1 = BuildT1(fw,n_t)
% Build a Toeplitz Matrix of coefficients of f(x).
% T1 \in \mathbb{R}^{(m+n-k+1)\times(n-k+1)}
%
%
%                           Inputs.
%
%
% fx : coefficients of polynomial f
%
% n : degree of polynomial g
%
% t :  index of subresultant S_{t} to be formed. (Also degree of GCD)
%
%
%%

% Get degree of polynomail f
[nRows,~] = size(fw);
m = nRows-1;

% Initialise empty matrix T1, for storing Toeplitz T_{k}(f)
T1 = zeros(m+n_t+1,n_t+1);

% Get f(w) with binomial coefficients;
fw_bi = fw.* GetBinomials(m);

% for each column of T1
for j = 0:1:n_t
    T1(j+1:m+j+1,j+1) = fw_bi;
end


end