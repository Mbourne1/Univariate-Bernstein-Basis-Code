function C1a = BuildToeplitz_fromPrev(m,n,k,C_0a)
% Build Toeplitz Matrix (Partition of Sylvester Matrix) in the 'Build Up' method.
% m : degree of polynomial f n : degree of polynomial g k : index of
% subresultant S_{k} to be built. C_0a : the preceeding Partition, from
% which C1a will be built. BOOL_DENOM : BOOL_LOG : Unused.

global BOOL_DENOM_SYL


% Where C_0 is the previous Cauchy Matrix,
knew = k+1;
% Build Matrix A
A = [zeros((m+n-k),1) diag(1./(1:1:(m+n-k))) ];

% Build Matrix B
Ba = [...
    zeros(1,n-k);...
    diag(1:1:(n-k))
    ];

% Build C1
C1a = A * C_0a * Ba;

%
switch BOOL_DENOM_SYL
    case 1
        % if Denominator is included in building toeplitz, then update
        % denominator for next S_k
        ratio = nchoosek(m+n-(k),n-(k)) ./ nchoosek(m+n-(knew),n-(knew));
        C1a = C1a .* ratio;
    case 0
end

warning('on')

end
