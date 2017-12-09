function vDegree = GetDegree_Array(arr_fx)
% Get the degree of an array of polynomials f_{i}(x)


% Get number of polynomials in the array
nPolys = length(arr_fx);

% Initialise a vector to store the degrees of each polynomial
vDegree = zeros(nPolys,1);


for i = 1 : 1: nPolys
    
    fx = arr_fx{i};
    
    vDegree(i) = GetDegree(fx);
    
end



end