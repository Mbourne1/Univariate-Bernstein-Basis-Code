function [GM] = GetGeometricMeanFromPrevious(fx,GM_prev,m,n_k)



p1 = (m+n_k+1) ./ (n_k+1);
p1 = p1 .^((n_k+2)./(n_k+1));

p2 = GM_prev .^((n_k+2)./(n_k+1));

p3 = nchoosek(m+n_k,m) .^ (1./(n_k+1));

prod_a = 1;
prod_b = 1;
for i = 0:1:m
    prod_a = prod_a * abs(fx(i+1));
    prod_b = prod_b * nchoosek(n_k+1+i,i);
end

p4 = (1./(prod_a * prod_b * prod_b)) .^(1./((n_k+1)*(m+1)));

GM = p1 * p2 * p3 * p4;

end
