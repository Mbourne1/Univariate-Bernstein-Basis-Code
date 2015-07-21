function C1 = BuildCauchy(m,n,k)

C1 = zeros(m+n-k+1,n-k+1);

% for each column
for j = 0:1:(n-k)
    % for each coefficient
    for i = j:1:j+m
        C1(i+1,j+1) = nchoosek(m+n-k-i,n-k-j) .* nchoosek(i,j) ./ ...
            nchoosek(m+n-k,n-k);
    end
end

end