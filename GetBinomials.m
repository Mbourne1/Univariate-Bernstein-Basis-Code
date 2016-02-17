function Bi_m = GetBinomials(m)

Bi_m = zeros(m+1,1);

for i=1:1:m+1
    Bi_m(i) = nchoosek(m,i-1);
end
end