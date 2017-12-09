

m = 5;


for m = 1 : 1 : 10

    
    fx = rand(m+1,1);

    vLambda(m) = GetLambda(m,fx);
    
    vMu(m) = GetMu(m,fx);
    
    
end

figure()
hold on
plot(vLambda)
plot(vMu)
hold off

display(vLambda)
display(vMu)
display(vLambda./vMu)

function lambda = GetLambda(m, fx)

prod = 1;
for i = 0:1:m 
    prod = prod .* fx(i+1) * nchoosek(m, i) ;
end



lambda = prod.^ (1./ (m + 1)) ;
end

function mu = GetMu(m, fx)

prod = 1;
for i = 0:1:m-1
    prod = m * (fx(i+2) - fx(i+1)) * nchoosek(m-1,i);
end

mu = prod.^ (1./ m);

end