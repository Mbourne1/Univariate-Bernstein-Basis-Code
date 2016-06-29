function [] = test

myArray = zeros(10,1);

parfor i = 1:1:10
    myArray(i) = 2*i;
end

myArray


end