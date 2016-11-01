function V = BuildV(m,n,k)
% Build the matrix V

V1 = BuildPartitionV(m,n-k);

V = [V1];
end

function V1 = BuildPartitionV(m,n_k)

vec = m+n_k+1:-1:1;
vec = (m+n_k+1) ./ vec; 

mat = diag(vec);
mat = [mat zeros(m+n_k+1,1)];

V1 = mat;
    

end

