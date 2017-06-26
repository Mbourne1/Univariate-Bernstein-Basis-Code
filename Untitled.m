
m = 5;
n = 8;
k = 2;


prod_a = 1;
prod_b = 1;
prod_c = 1;

for i = 0 : 1 : m 
    for j = 0:1:n-k

        prod_a = prod_a * nchoosek(m,i);
        prod_b = prod_b * nchoosek(n-k,j);
        prod_c = prod_c * nchoosek(m + n - k , i + j);
        
    end
end


display(prod_a)
display(prod_b)
display(prod_c)

prod_a * prod_b