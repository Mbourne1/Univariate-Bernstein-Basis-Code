m = 10;
n = 8;
k = 3;




sum1 = 0;

for i = 0 : 1 : m - 1
    
    for j = 0 : 1 : m - k
   
        sum1 = sum1 +  nchoosek(i + j, i);
        
    end
end


sum2 = 0;

for i = 0 : 1 : m
    
    for j = 0 : 1 : m - k - 1
   
        sum2 = sum2 +  nchoosek(i + j, i);
        
    end
end

display(sum1)
display(sum2)

sum1 - sum2


sum3 = 0;

for i = 0 : 1 : m - 1
   
    sum3 = sum3  + nchoosek(m - k + i, i);
    
end

sum4 = 0;

for j = 0 : 1 : m - k - 1
   
    sum4 = sum4 + nchoosek(m + j, j);

end

display(sum3)
display(sum4)

sum3 - sum4