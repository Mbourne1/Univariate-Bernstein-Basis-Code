function H1 = BuildH1(m)

% Initialize a zero vector
H1 = zeros(m+1,1);
for i = 0:1:m
    H1(i+1) = nchoosek(m,i);
end

H1 = diag(1./H1);

end