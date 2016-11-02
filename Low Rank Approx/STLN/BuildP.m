function P = BuildP(idx_Col,m,n,t)
% Build the matrix P_{t}, where h_{t} = P_{t}z

if idx_Col <= n-t+1
    % First Partition
    j = idx_Col;
    
    % Get binomials corresponding to f(x)
    bi_m = GetBinomials(m);
    
    bi_denom = zeros(m+1,1);
    for k = 0:1:m
        bi_denom(k+1) = nchoosek(m+n-t,(j-1)+k);
    end
    
    
    G = bi_m .* nchoosek(n-t,j-1) ./ bi_denom;
    
    P = ...
        [
        zeros(j-1,m+1)      zeros(j-1,n+1);
        diag(G)        zeros(m+1,n+1);
        zeros(n-t-j+1,m+1)  zeros(n-t-j+1,n+1);
        ];
else
    % Second Partition
    j = idx_Col - (n-t+1);
    
    
    bi_n = GetBinomials(n);
    
    bi_denom = zeros(n+1,1);
    for k = 0:1:n
        bi_denom(k+1) = nchoosek(m+n-t,j+k-1);
    end
    
    G = (bi_n .* nchoosek(m-t,j-1)) ./ bi_denom;
    
    P = ...
        [
        zeros(j-1,m+1)      zeros(j-1,n+1);
        zeros(n+1,m+1)      diag(G) ;
        zeros(m-t-j+1,m+1)  zeros(m-t-j+1,n+1);
        ];
    
end

end